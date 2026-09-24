# BtoM：反演磁铁内部的磁化强度分布

**B**-field **to** **M**agnetization。用磁铁**外部**测到的磁场数据 $\mathbf B$，
反推出磁铁**内部**的磁化强度分布 $\mathbf M$ —— 不解剖、不切样品的无损表征。

正问题是线性的：

$$\mathbf B(\mathbf r_i)=\frac{\mu_0}{4\pi}\sum_j V_j\,
\frac{3(\mathbf M_j\cdot\hat{\mathbf d}_{ij})\hat{\mathbf d}_{ij}-\mathbf M_j}{|\mathbf d_{ij}|^{3}},
\qquad \mathbf d_{ij}=\mathbf r_i-\mathbf r_j$$

每个体素当成一个点偶极子，偶极矩是 $\mathbf M_jV_j$。反问题就是解 $\mathbf A\mathbf M=\mathbf B$。

难点不在"算不算得出来"，而在**解不唯一**：能拟合同一批 $B$ 的 $M$ 有无穷多组。
本项目关心的就是"在无穷多组解里怎么挑出物理上对的那一组"，以及挑解规则必须
**与网格、与单位无关**。

<div class="grid cards" markdown>

-   **方法一：反演方法（B2M）**

    ---

    点偶极正演的离散（$B=A M$）、`discretize.TreeMesh` 自适应网格、
    边缘保持正则化、以及用 $\sqrt{v_i}\,m_i$ 当未知量来消掉粗细网格的区别。

    [阅读](methods/b2m.md)

-   **方法二：对称面检测与镜像**

    ---

    小探头测穿了对称面之后，手工给的镜像面高度会把镜像点叠到实测点上。
    用数据自身定出对称面 $h$，并让**实测值永远优先**。

    [阅读](methods/symmetry.md)

-   **理论推导**

    ---

    从正问题的离散出发，证明"按体素体积加权"才是 $\int\rho(|\nabla M|)\mathrm dV$
    的正确离散化；以及为什么 $\lambda$ 取值不当会表现为"细网格偏小、粗网格偏大"。

    [阅读](theory/index.md)

-   **对照实验**

    ---

    E0–E7 七组可复现实验、29 条判据：$O(h^{2})$ 收敛性、$\lambda$ 换算表、
    欠正则化幅值指纹、网格无关性、$\lambda$ 自动选择、粗细网格互换不变性。

    [阅读](experiments/index.md)

-   **数据与结果**

    ---

    三套测量数据集、几何参数、以及反演出来的磁化强度分布截面图。

    [阅读](results/new-prob.md)

-   **复现指南**

    ---

    目录结构、Python 环境、每个 notebook / 脚本 / 实验的命令行。

    [阅读](reproduce.md)

</div>

## 一句话结论

| 问题 | 结论 |
|---|---|
| 正则化该不该按体素体积加权？ | **该**。$w_f=S_fh_f$（面对偶体积）是 $\int\rho(\|\nabla M\|)\mathrm dV$ 的正确求积权重，误差 $O(h^2)$；不加权那版随加密**发散**（实测 338% → 3467% → 28628%） |
| 那为什么"加了体积项反而结果非物理"？ | 因为 $\lambda$ 的**标度**被换掉了。加权重后同一物理强度对应的 $\lambda$ 要乘 $W/n_f$（平均面对偶体积，生产网格上 ≈5；若把体积换成 m³ 则是 $10^9$）。$\lambda$ 没跟着换算就掉进欠正则化区，解就堆到大格子上 |
| 怎么根治"粗格子偏大、细格子偏小"？ | 换未知量：求解 $\sqrt{v_i}\,m_i$ 而不是 $m_i$。这样求解器隐式最小化的是 $\sum v_im_i^2=\int m^2\mathrm dV$（物理范数）。实测互换粗细网格的误差从 **58%~89% 降到 ≤0.7%** |
| $\lambda$ 怎么定？ | 先力平衡标定 $\lambda_{bal}$，再按"残差 ≈ 噪声 $\sigma$"取最大的那个 $\lambda$。实测与上帝视角最优 $\lambda$ 完全一致 |

!!! warning "两个必须知道的限制"
    1. **新数据集（`UshapeNormal_New_prob`）的磁铁几何目前是占位值**（从扫描范围反推的粗估），
       厚度 $T$ 必须用卡尺实测。用占位几何算出的结果是**临时**的，本站的定性结论不受影响，
       但绝对幅值不要引用。
    2. **$\sqrt{v_i}\,m_i$ 换变量目前只在 `meshtest/` 里验证过，尚未接进生产 notebook**。
       生产 notebook 的未知量仍是 $m_i$。

## 记号约定

| 符号 | 含义 | 代码里 |
|---|---|---|
| $\mathbf B$ | 测点处的磁通密度（T） | CSV 里是 µT，乘 `B_unit_conversion = 1e-6` |
| $\mathbf M$ / $m_i$ | 体素磁化强度（A/m） | 未知量 |
| $v_i$ / $V_c$ | 体素体积 | `mesh.cell_volumes` |
| $w_f$ | 面对偶体积 $S_fh_f$ | `mesh.get_face_inner_product().diagonal()` |
| $\lambda$ | 正则化强度 | `lambda_reg`（生产取 `1e-10`） |
| $\varepsilon$ | 正则化"边缘尺度" | `huber_epsilon`（生产取 `4e3` A/m） |
| $h$ | 磁铁厚度方向的对称面高度（mm） | `SYMMETRY_H`，由 `symmetry.py` 自动检测 |

!!! note "关于本站"
    站点用 **MkDocs + Material for MkDocs** 构建，公式由 **MathJax 3** 在浏览器端渲染。
    `理论推导` / `对照实验` 两栏是仓库里 `meshtest/*.md` 的**自动同步副本**
    （见 `web/tools/sync_docs.py`），事实来源始终是仓库中的源文件。
