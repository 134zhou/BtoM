# BtoM —— 反演磁铁内部的磁化强度分布

**B**-field **to** **M**agnetization。用磁铁**外部**测到的磁场数据 $\mathbf B$，
反推出磁铁**内部**的磁化强度分布 $\mathbf M$ —— 无损、不解剖。

正问题是线性的（每个体素当点偶极子，偶极矩 $\mathbf M_jV_j$）：

$$\mathbf B(\mathbf r_i)=\frac{\mu_0}{4\pi}\sum_j V_j\,
\frac{3(\mathbf M_j\cdot\hat{\mathbf d}_{ij})\hat{\mathbf d}_{ij}-\mathbf M_j}{|\mathbf d_{ij}|^{3}},
\qquad \mathbf d_{ij}=\mathbf r_i-\mathbf r_j$$

反问题 $\mathbf A\mathbf M=\mathbf B$ 的难点不是"算不出来"，而是**解不唯一**：
能拟合同一批 $B$ 的 $M$ 有无穷多组。本仓库关心的就是"在无穷多组解里怎么挑出物理上
对的那一组"，以及挑解规则必须**与网格、与单位无关**。

📖 **在线文档：<https://134zhou.github.io/BtoM/>**

---

## 三条主要结论

### 1. 正则化必须按体素体积加权 —— 它才是"正确的离散化"

$w_f=S_fh_f$（面对偶体积，`mesh.get_face_inner_product().diagonal()`）是
$\int\rho(\lvert\nabla \mathbf M\rvert)\,\mathrm dV$ 的正确求积权重，误差 $O(h^2)$；
不加权的写法随网格加密**发散**（$O(h^{-3})$）：

| 每方向基础网格数 | 体积加权误差 | 不加权误差 |
|---|---|---|
| 16 | 5.01% | 338.8% |
| 32 | 1.68% | 3467% |
| 64 | **0.46%** | **28628%** |

### 2. "加了体积项结果就非物理"的真因是 λ 的标度被换掉了

加权重后，**同一物理强度对应的 λ 要乘 $W/n_f$**（平均面对偶体积，正比于 $h^3$）：

| 网格 | $W/n_f$（λ 需要乘的倍数） |
|---|---|
| 本实验 quick（4/8 mm） | ≈108 |
| 本实验 full（2/4 mm） | ≈14 |
| 生产网格（1.0~1.4 mm） | **≈5** |
| 体积写成 m³ | **×10⁹** |

λ 没跟着换算 → 等效正则化强度掉进欠正则化区 → 磁化强度堆到大单元上，
表现为**细网格偏小、粗网格偏大**。可用窗口是单边的 $k\in[0,+1]$，**宁大勿小**。

### 3. 根治办法：换未知量，求解 $\sqrt{v_i}\,m_i$ 而不是 $m_i$

求解器从 0 出发时隐式最小化 $\sum x_i^2$：

| 未知量 | 隐式最小化 | 等效权重 | 结果 |
|---|---|---|---|
| $m_i$ | $\sum m_i^2$ | $v^{0}$ | 偏粗格子 |
| $\sqrt{v_i}\,m_i$ | $\sum v_im_i^2=\int m^2\mathrm dV$ | $v^{1}$ | **物理范数，无偏** ✔ |

实测（$\lambda=0$，完全无正则化；真值均匀磁化，左右两半区必须相等）：

| 变量 | 网格A(左细) 左/右 | 网格B(右细) 左/右 | 互换不变性误差 |
|---|---|---|---|
| $m$（LSQR@50 / L-BFGS-B@50） | 0.63 / 0.53 | 1.58 / 1.89 | **58% / 89%** |
| $\sqrt v\,m$ | **1.00 / 1.00** | **1.00 / 1.00** | **0.1% / 0.0%** |

---

## 目录结构

```
BtoM/
├── src/                       公用代码
│   ├── symmetry.py            对称面自动检测 + 镜像重合处理（库 + CLI + 自测）
│   ├── csv2vtk.py             CSV -> 结构化 VTK 点云
│   └── build_notebooks.py     由 .py（# %% 标记）生成 .ipynb
├── notebooks/                 全部 notebook 与对应的 .py 源
│   ├── B2M.ipynb              生产版：specialU 数据集
│   ├── B2M_noreg.py/.ipynb    去掉正则化的对照
│   ├── B2M_newprob.py/.ipynb  新版小探头数据集 + 对称面自动检测
│   └── plot_2d_sections.py/.ipynb  文章配图
├── meshtest/                  ⭐ 对照实验与数学推导（29 条判据全过）
│   ├── THEORY.md              数学推导（含"大白话"开篇）
│   ├── README.md              实验说明与结论
│   ├── results.md             自动生成的判据汇总
│   ├── b2m_core.py truth.py run_experiments.py
│   └── figures/               （不入库，跑实验后生成）
├── data/
│   ├── raw/                   ✅ 原始测量数据（入库）
│   └── derived/               ❌ 反演产物 .vtu/.vtk（不入库，可重跑）
├── figures/                   ❌ 论文/站点配图源文件（不入库，手工挑选后放进站点）
├── archive/                   历史文件（旧 notebook、早期试验、参数快照）
├── docs/                      ⭐【构建产物】GitHub Pages 发布的目录（由 mkdocs build 生成）
└── web/                       站点工程
    └── docs/                  站点源码：所有分页 markdown 平铺在这里
```

> ⚠️ 两个 `docs/` 别搞混：`web/docs/` 是 **markdown 源**（人写），
> 仓库根的 `docs/` 是 **HTML 产物**（`mkdocs build` 生成，要提交，别手改）。
> Pages 走 **Settings → Source: Deploy from a branch → `main` / `/docs`**，
> 不需要 `gh-pages` 分支。

所有脚本/notebook 开头都有一段**仓库路径引导**（向上找含 `.git` 的仓库根），
所以**在仓库里的任何位置运行都能找到同一批文件**。

---

## 快速开始

```bash
PY=E:/Python/Miniforge/envs/discretize/python.exe     # 算法环境，依赖已装好

# 反演（新版小探头数据集，含对称面自动检测）
$PY notebooks/B2M_newprob.py

# 反演（旧数据集，生产版 notebook）
$PY -m jupyter nbconvert --to notebook --execute notebooks/B2M.ipynb --inplace

# 对称面检测 + 报告 + 出图
$PY src/symmetry.py data/raw/UshapeNormal_New_prob

# 对照实验（29 条判据，快速约 1 分钟）
cd meshtest && $PY run_experiments.py --quick
```

文档站的构建与发布见 [`web/README.md`](web/README.md)。

---

## 数据集

| 目录 | 测点数 | 说明 |
|---|---|---|
| `data/raw/specialU/` | 3696 + 2816 + 720 | 早期"特殊 U"磁铁 |
| `data/raw/UshapeNormal/` | 4200 + 7056 | U 形磁铁，圆弧区 4 mm、磁腿区 2 mm 间距 |
| `data/raw/UshapeNormal_New_prob/` | **48784**（6 个文件） | **小探头**重测，测点跨越磁铁对称面 |
| `data/raw/legacy/` | 975 + 975 | 更早的两组测量，仅存档 |

对称面自动检测结果（新数据集）：$h_{\text{raw}}=7.698$ mm → 吸附到 **$h=8.000$ mm**，
重合 **18405/48784 = 37.7%** 的镜像副本被丢弃，最终 **79163** 点。

---

## 状态与限制

- ⚠️ **`UshapeNormal_New_prob` 的磁铁几何是占位值**（从扫描范围反推的粗估），
  厚度 $T$ 必须用卡尺实测。该数据集的反演绝对幅值暂时不要引用。
- ⚠️ **$\sqrt{v_i}\,m_i$ 换变量只在 `meshtest/` 里验证过，尚未接进生产 notebook**。
- 真实磁铁没有真值，所以方法学结论全部建立在 `meshtest/` 的**合成真值对照实验**上
  （独立 1 mm 细网格正演，避免 inverse crime）。
- 未做：各向异性 TV、L0/L1 稀疏先验、真实 U 形几何 + 真实点云上的 λ 标定值。

## 踩过的坑

完整清单在 [`meshtest/README.md`](meshtest/README.md) 的"踩过的坑"一节（8 条）。最容易再犯的：

1. **几何必须与网格对齐**，否则按"单元中心在盒内"选出的磁体体积会偏大（踩过 +39%）。
2. **λ 不能用绝对数值扫**，必须先用 $\lambda_{bal}$ 标定，否则整段落在过正则化区。
3. **换变量时只有正则化项要除 $\sqrt v$**；数据项梯度就是 $2\mathbf A^{\mathsf T}\mathbf r$。
   两项都除会让 L-BFGS-B 第 0 次迭代就 `ABNORMAL`、解恒为 0。
4. **$\sqrt v$ 要用 SI**（$v$ 以 m³ 计），混用 mm³ 会让变量尺度到 $10^5$、优化器失效。
5. **测点离磁体要足够远**（≥4~6 个单元尺度），否则点偶极近似本身就不可用。

## 许可

暂无。如需引用或复用请先联系仓库作者。
