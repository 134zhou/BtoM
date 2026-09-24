# 反演方法（B2M）

代码：`notebooks/B2M.ipynb`（旧数据集）、`notebooks/B2M_newprob.py` / `.ipynb`（新小探头数据集）。
两版的正演、正则化、求解完全一致，只有"数据装配"（对称面与镜像）不同。

---

## 1. 正问题：每个体素一个点偶极子

离散后未知量是每个体素的磁化强度 $\mathbf m_j$，它对测点 $\mathbf r_i$ 的贡献是点偶极场：

$$\mathbf B(\mathbf r_i)=\frac{\mu_0}{4\pi}\sum_{j} V_j\,
\frac{3(\mathbf m_j\cdot\hat{\mathbf d}_{ij})\hat{\mathbf d}_{ij}-\mathbf m_j}{|\mathbf d_{ij}|^{3}},
\qquad \mathbf d_{ij}=\mathbf r_i-\mathbf r_j$$

写成矩阵形式就是 $\mathbf B=\mathbf A\mathbf m$，$\mathbf A$ 是 $3M\times3N$ 的稠密矩阵。
本项目**不显式构造** $\mathbf A$：用 `scipy.sparse.linalg.LinearOperator` 只提供

- `matvec(v)` ＝ $\mathbf A\mathbf v$（正演，求场）
- `rmatvec(w)` ＝ $\mathbf A^{\mathsf T}\mathbf w$（伴随，反投影）

两个核都用 `numba.njit(parallel=True)` 重写并 `prange` 并行。这是整个程序最耗时的地方，
也是"点偶极近似"最直接的体现。

!!! note "单位"
    几何用 **mm**、体积用 **mm³**，但偶极公式里的 $\mu_0/4\pi=10^{-7}$ 是 SI 常数，
    所以做正演前要把坐标换成 **m**、体积换成 **m³**（乘 $10^{-9}$），得到的 $\mathbf B$ 是 **T**。
    CSV 里的磁场是 **µT**，读入后乘 `B_unit_conversion = 1e-6`。

---

## 2. 网格：只在磁铁里放未知量

用 `discretize.TreeMesh` 建八叉树自适应网格（生产配置里 `voxel_size_base = 1.5` mm）：

1. 先把**磁铁外接盒**按基础体素铺一层，再 `refine_bounding_box` 把盒内整体加密；
2. 再按几何条件把**磁铁实体所在的单元**加密到最细档；
3. 未知量只取**单元中心落在磁铁实体内**的那些单元 —— 外面是空气，$\mathbf M\equiv0$，不设未知量。

旧数据集里还用过"按 $y$ 分界、左右两半不同细化档"的网格（`Y_split = 144`），
用来检验解对网格的依赖性 —— 这正是后面「换未知量」那组实验的原型。

```python
# notebooks/B2M.ipynb（节选）
mesh = discretize.TreeMesh([...], origin=..., diagonal_balance=False)
mesh.refine_bounding_box(BBox, level=absolute_coarse_level, finalize=False)
mesh.refine(refine_func, finalize=False)     # 盒内加密
mesh.finalize(); mesh.number()
voxel_points = mesh.cell_centers[sel_magnet]   # 只保留磁铁内的单元
voxel_vol    = mesh.cell_volumes[sel_magnet]
```

---

## 3. 目标函数：数据项 + 边缘保持正则化

$$\min_{\mathbf m}\;\;
\underbrace{\bigl\|\mathbf A\mathbf m-\mathbf B\bigr\|^{2}}_{E_{\text{data}}}
\;+\;\lambda\,\varepsilon^{2}\sum_f w_f
\left(\sqrt{1+\left(\frac{\|\nabla \mathbf m\|_f}{\varepsilon}\right)^{2}}-1\right)$$

| 项 | 写法 | 作用 |
|---|---|---|
| 数据项 | $\|\mathbf A\mathbf m-\mathbf B\|^2$ | 拟合测量值 |
| 罚函数 | $\varepsilon^2\bigl(\sqrt{1+x^2/\varepsilon^2}-1\bigr)$ | 伪 Huber / Charbonnier：$x\ll\varepsilon$ 时是二次的（平滑），$x\gg\varepsilon$ 时是线性的（保边缘） |
| $w_f$ | 面对偶体积 $S_fh_f$ | **按体素体积加权的关键**，见下一节 |
| 约束 | $\|\mathbf m\|\le$ `max_M` $=1.5\times10^{6}$ A/m | 盒式约束，量级取材料饱和磁化强度 |

$\nabla\mathbf m$ 是 TreeMesh 的**面梯度** `mesh.cell_gradient`（每对相邻单元一个面）。
梯度由解析式给出（`jac=True`），交给 **L-BFGS-B** 求解：

```python
loss_data = np.sum(r**2)
loss_reg  = lambda_reg * eps**2 * np.sum(
    V_face * (np.sqrt(1.0 + (grad_norm/eps)**2) - 1.0))
grad_data = 2.0 * A_op.rmatvec(r)
grad_reg  = lambda_reg * (G_sub.T @ (V_face[:, None] * grad * weights[:, None])).ravel()
```

生产参数：`lambda_reg = 1e-10`、`huber_epsilon = 4e3`、`max_iter = 50`。

---

## 4. 为什么 $w_f$ 必须按体积加权

$w_f=S_fh_f$ 是 $\int\rho(\|\nabla \mathbf M\|)\,\mathrm dV$ 的**正确求积权重**：

$$\sum_f w_f\,\rho(\|\nabla \mathbf m\|_f)\;\xrightarrow[\;h\to0\;]{}\;\int\rho(\|\nabla \mathbf M\|)\,\mathrm dV
\quad\text{（误差 }O(h^{2})\text{）}$$

而不加权的 $\sum_f\rho(\|\nabla\mathbf m\|_f)$ 随网格加密**发散**（$O(h^{-3})$），
也就是说它本身就是一个网格相关的量 —— 同一个磁铁切细一点，正则化项的数值就变了。

受控合成实验（面梯度型泛函 vs 解析积分，网格逐级加密）：

| 每方向基础网格数 | 体积加权误差 | 不加权误差 |
|---|---|---|
| 16 | 5.01% | 338.8% |
| 32 | 1.68% | 3467% |
| 64 | **0.46%** | **28628%** |

这两列不是调参调出来的：$w_f$ 就是体积积分在每个面上的求积权重，所以加权那列随加密收敛到
解析值；不加权那列等价于把权重全取成 1，而网格越细内部面越多，它自然按 $O(h^{-3})$ 发散 ——
也就是说，不加权的正则化项**本身就是一个网格相关的量**。

---

## 5. 未知量选 $m_i$ 还是 $\sqrt{v_i}\,m_i$？

这是本项目最新的一条结论，也是消掉"粗网格偏大、细网格偏小"的**根治手段**。

求解器从 $\mathbf m=0$ 出发，隐式挑的是"未知量平方和最小"的那组解，于是：

| 未知量 $x_i$ | 隐式最小化的量 | 等效给 $m$ 的权重 | 结果 |
|---|---|---|---|
| $m_i$ | $\sum m_i^2$ | $v^{0}$ | **偏粗格子**（细格子被压小） |
| $v_im_i$ | $\sum (v_im_i)^2$ | $v^{2}$ | 偏到另一侧（细大粗小） |
| $\sqrt{v_i}\,m_i$ | $\sum v_im_i^2=\int m^2\,\mathrm dV$ | $v^{1}$ | **物理范数，无偏** ✔ |

验法只用一条判据（$\lambda=0$，完全没有正则化）：真值取整块均匀磁化、磁场无噪声，
造两张**互为粗细互换**的网格（A：左半细；B：右半细），比较左右半区的体积加权平均 $|m|$。
物理上左右必须相等：

| 设置 | 变量 | 网格A(左细) 左/右 | 网格B(右细) 左/右 | 互换不变性误差 |
|---|---|---|---|---|
| LSQR@50 | $m$ | 0.63 | 1.58 | 58.0% |
| LSQR@50 | $\sqrt v\,m$ | **1.00** | **1.00** | **0.1%** |
| LSQR@200 | $m$ | 0.70 | 1.43 | 43.4% |
| LSQR@200 | $\sqrt v\,m$ | **0.99** | **1.01** | **0.7%** |
| L-BFGS-B@50 | $m$ | 0.53 | 1.89 | 89.4% |
| L-BFGS-B@50 | $\sqrt v\,m$ | **1.00** | **1.00** | **0.0%** |

$\sqrt v\,m$ 解出的区域均值是 $99953\sim100407$（真值 $10^5$）—— **绝对值也对**。

!!! warning "尚未接入生产"
    受控实验里已经有 `variable="m" | "sqrt_v_m"` 的开关（`VARIABLES`、`variable_scale()`、
    `to_m()`），但**生产 notebook 的未知量仍然是 $m_i$**。
    换变量的实现要注意一个坑：**只有正则化项的梯度要乘 $\partial m/\partial x$，数据项不乘**
    （数据项梯度就是 $2\mathbf A^{\mathsf T}\mathbf r$）。两项都乘会让 L-BFGS-B 第 0 次迭代
    就 `ABNORMAL`、解恒为 0；`b2m_core.check_gradient` 用有限差分守着这条。

---

## 6. $\lambda$ 怎么定

$\lambda$ 的**绝对值没有跨配置意义**（它取决于 $w_f$ 的量纲与归一化），所以分两步：

1. **力平衡标定**：$\lambda_{bal}=\dfrac{\|\nabla E_{\text{data}}(\mathbf M_{\text{ref}})\|}{\|\nabla E_{\text{reg}}(\mathbf M_{\text{ref}})\|_{\lambda=1}}$，
   参考模型取"材料饱和磁化强度量级 × 正弦扰动"即可，**不需要真值**；
2. **按噪声定档**：在 $\lambda=\lambda_{bal}\cdot10^{k}$ 里取满足"残差 ≈ 噪声 $\sigma$"的
   **最大** $\lambda$（discrepancy principle）。

受控实验：这样选出来的 $\lambda$ 与上帝视角最优完全一致（$k=+1$ vs $+1$，误差同为 0.29%）。
可用窗口是**单边**的 $k\in[0,+1]$，**宁大勿小**。

---

## 7. 已知问题与 TODO

| 事项 | 现状 |
|---|---|
| 生产正则化用的是 `volume_raw`（$w_f$ 以 mm³ 计），$\lambda=10^{-10}$ | 能跑，但 $\lambda$ 会随网格加密漂移 ≈5 倍。推荐改成 `w_f/W` 归一化（`volume_norm`），或至少记住这个换算 |
| 未知量仍是 $m_i$ | 建议换 $\sqrt{v_i}\,m_i$，见 §5 |
| 新数据集的磁铁几何是**占位值** | 必须用卡尺实测厚度 $T$ 并替换 `MAGNET_*`，然后置 `GEOMETRY_IS_PLACEHOLDER = False` |
| 点偶极近似 | 测点离磁体表面太近（< 4~6 个单元尺度）时误差不可忽略；受控实验里踩过"真值处残差 = 90σ"的坑 |
