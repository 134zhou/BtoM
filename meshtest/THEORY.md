# 正则化体积加权的数学依据

**从连续正问题到离散反问题：为什么必须按面对偶体积加权，为什么权重取错会导致欠正则化**

本文档给出 `meshtest` 数值结论的数学依据。所有定理都给出精确陈述、假设与证明（或证明骨架＋附录依据），
并在每条后标注对应的数值证据编号（证据来源见附录 C，复现命令见附录 D）。

---

## 开篇：用大白话说清楚（只用 v、m、B、A、λ，高中数学足够）

**一句话**：$B$ 只认得每个格子的"**体积 × 磁化强度**"，认不出"磁化强度"本身；格子怎么切是我们自己选的，
$B$ 里根本没有这个信息。数据管不了的那些自由度，由"谁来挑一个解"决定——$\lambda$ 太小时，
挑解的是求解器的默认习惯（让 $\sum m_i^2$ 最小），而这个习惯**偏爱大格子**。

### 引子 1 · $B$ 只认识"体积 × 磁化强度"

把磁铁切成若干格子，第 $i$ 个格子的体积是 $v_i$、磁化强度是 $m_i$。第 $i$ 个格子对磁场的贡献
正比于 $v_i\,m_i$（偶极子公式里体积就是积分权重），所有格子叠加：

$$
B\;=\;\sum_i \bigl(v_i\,m_i\bigr)\times(\text{几何因子}).
$$

几何因子只跟"格子在哪里、测点在哪里"有关，跟格子多大无关。所以 $B$ 真正认识的量是**磁矩** $v_i m_i$，
而不是 $m_i$ 本身。

### 引子 2 · 同一个磁铁，换一套格子，$B$ 一模一样

同一个磁铁，切成粗格子还是细格子，测出来的 $B$ 是同一个 $B$（这是物理，不是近似）。
所以**网格只是我们自己记账的方式**：同样一个磁矩，既可以记成"大格子里的一个磁化强度"，
也可以记成"几个小格子里的磁化强度"。$B$ 对此一无所知。

**推论**：反演里 $m_i$ 只以 $v_i m_i$ 的形式出现，所以"解 $A(v_i m_i)=B$"**解不出唯一的 $m$**，
只能解出磁矩；把磁矩换算回磁化强度要除以 $v_i$，而 $v_i$ 是我们自己定的。

### 引子 3 · 为什么有一大堆解都能拟合 $B$

测点在磁铁外面，离得越远磁场越"糊"。两个挨着的格子（一粗一细）在同一个测点上产生的磁场形状几乎一样
（差别是"格子尺寸 ÷ 距离"的量级）。所以只要把磁矩在这两个格子之间挪一挪、合计不变，$B$ 几乎不变。
**也就是说：数据在这些方向上几乎没有约束力。**

### 引子 4 · 那这些自由度是谁填的

反演从 $m_i=0$ 开始迭代。第一步往哪个格子放多少，取决于"改这个格子的 $m$ 对 $B$ 的影响有多大"，
而由 0.1 这个影响正比于 $v_i$——**大格子总是先被填、填得多**。
更准确地说：从 0 出发的梯度类算法（L-BFGS、LSQR 都是）在"能拟合 $B$ 的所有解"里，
挑的是 $\sum_i m_i^2$ **最小**的那个。$\lambda$ 足够大时，正则化项会把这个默认习惯压过去；
$\lambda$ 太小（或者权重取错）时，默认习惯说了算。

### 引子 5 · 为什么"最小 $\sum m_i^2$"就偏大格子：两个格子手算

粗格子体积 $8v$、细格子体积 $v$。真值是均匀磁化：两个格子都是 $m$，总磁矩 $=8v\cdot m+v\cdot m=9vm$。

| 方案 | 粗格子 $m$ | 细格子 $m$ | $\sum m_i^2$ | 体积加权 $\sum v_i m_i^2$ |
|---|---|---|---|---|
| 甲：真值 | $m$ | $m$ | $m^2+m^2=2m^2$ | $8v m^2+v m^2=9vm^2$ |
| 乙：全塞进粗格子 | $9vm/(8v)=1.125m$ | $0$ | $(1.125m)^2=1.27m^2$ | $8v(1.125m)^2=10.1vm^2$ |

- 按 $\sum m_i^2$ 挑：$1.27<2$，选**乙**——粗格子偏大（$1.125m$）、细格子偏小（$0$），
  这就是"粗大细小"的来历。
- 按 $\sum v_i m_i^2$ 挑：$9<10.1$，选**甲**（真值）。

### 引子 6 · 一般结论（任意多个格子）

在"总磁矩固定"这个约束下：

- 让 $\sum_i m_i^2$ 最小（求导等于 0）：$2m_i=$ 常数 $\times v_i$，即 $m_i$ **正比于** $v_i$。
  体积是 8 倍的粗格子，磁化强度就是细格子的 **8 倍**（$8=2^3$，线尺度 2 倍）。
- 让 $\sum_i v_i m_i^2$ 最小：$2v_i m_i=$ 常数 $\times v_i$，即 $m_i$ **处处相等**——均匀磁铁该有的样子。

注意 $\sum_i m_i^2$ 本身**依赖切法**（把一格切成两格，这个和就变），所以它不是一个物理量；
$\sum_i v_i m_i^2$ 才是 $\int m^2\,\mathrm dV$ 的离散版，与切法无关。**用前者挑解，就会挑出网格的痕迹。**

### 引子 6.1 · 换变量能把这件事修掉（E7 实测）

既然偏置来自"求解器隐式最小化 $\sum x_i^2$"，那就换变量：令 $x_i=\sqrt{v_i}\,m_i$，于是
$$\textstyle\sum_i x_i^2=\sum_i v_i m_i^2=\int_\Omega m^2\,\mathrm dV ,$$
求解器隐式用的正好是物理范数。验法只用一条：真值取整块均匀磁化、磁场无噪声（物理上左右两半区
均值必须相等），造两张**互为粗细互换**的网格，看左右区均值。实测（$\lambda=0$，完全没有正则化）：

| 变量 | 网格A(左细) 左/右 | 网格B(右细) 左/右 | 互换不变性误差 |
|---|---|---|---|
| $m_i$（LSQR@50 / L-BFGS-B@50） | 0.63 / 0.53 | 1.58 / 1.89 | **58% / 89%** |
| $\sqrt{v_i}\,m_i$ | **1.00 / 1.00** | **1.00 / 1.00** | **0.1% / 0.0%** |

即：换变量后粗/细网格的区别消失，且对调粗细结果不变。（$\lambda=0$ 时结果仍受"迭代到哪一步"
影响，那是另一个问题，不在这里讨论。）

### 引子 7 · 与实测对照

| 情形 | 粗格子/细格子（纯数据拟合，$\lambda$ 极小） |
|---|---|
| 测点放到 200 mm 远（几何因子近似常数，上面的推断前提成立） | **7.997**（理论就是体积比 8） |
| 真实几何（测点 24 mm），几何能部分区分格子 | 1.716（方向不变，量级被削弱） |
| 改用体积加权规则挑解 | **1.000** |

### 引子 8 · 因此必须做两件事

1. **正则化必须按 $v_i$ 加权**：否则它自己就是网格相关的（同一个磁铁切细一点，正则化项的数值就变），
   而且它的偏向与引子 5 的偏向**同方向叠加**，把"粗大细小"放大。
2. **$\lambda$ 不能太小**：要让按 $v_i$ 加权的正则化项去决定那堆同样能拟合 $B$ 的解，
   而不是让求解器的默认习惯去决定。经验做法：取"残差刚好降到噪声水平"的那个 $\lambda$，
   或取力平衡标定值的 $1\sim10$ 倍。

**收尾一句**：不是"反着解不出来"，而是"反着解出来的解不唯一，挑哪一个"这件事——
$\lambda$ 太小时，挑解的是一个偏爱大格子的坏习惯。

---



## §0 记号与问题的数学设定

### 0.1 连续正问题

磁铁占据区域 $\Omega\subset\mathbb R^3$，内部磁化强度为矢量场 $\mathbf M:\Omega\to\mathbb R^3$。
在自由空间中，磁化强度产生的磁场（偶极子积分）为

$$
\mathbf B(\mathbf r)\;=\;\frac{\mu_0}{4\pi}\int_\Omega
\frac{3\bigl(\mathbf M(\mathbf r')\cdot\hat{\mathbf u}\bigr)\hat{\mathbf u}-|\mathbf u|^{2}\mathbf M(\mathbf r')}{|\mathbf u|^{5}}\,\mathrm dV',
\qquad \mathbf u=\mathbf r-\mathbf r',\quad \hat{\mathbf u}=\mathbf u/|\mathbf u| .
\tag{0.1}
$$

等价的核写法：$\mathbf B(\mathbf r)=\frac{\mu_0}{4\pi}\int_\Omega \mathsf K(\mathbf r,\mathbf r')\mathbf M(\mathbf r')\mathrm dV'$，
其中 $3\times3$ 对称核

$$
\mathsf K(\mathbf r,\mathbf r')=\frac{3\hat{\mathbf u}\hat{\mathbf u}^{\mathsf T}-\mathsf I}{|\mathbf u|^{3}} .
\tag{0.2}
$$

### 0.2 连续反问题

给定测点 $\mathbf r_i$（$i=1,\dots,N_{\text{meas}}$）与其上测得的磁场 $\mathbf b_i$，反演 $\mathbf M$：

$$
\min_{\mathbf M}\;\; \underbrace{\sum_{i=1}^{M}\bigl|\mathbf B[\mathbf M](\mathbf r_i)-\mathbf b_i\bigr|^{2}}_{=:E_{\text{data}}[\mathbf M]}
\;+\;\lambda\underbrace{\int_\Omega \rho\bigl(|\nabla \mathbf M|\bigr)\,\mathrm dV}_{=:R[\mathbf M]} .
\tag{0.3}
$$

罚函数取 Charbonnier/Huber 型（与代码一致）：

$$
\rho(x)=\epsilon^{2}\Bigl(\sqrt{1+(x/\epsilon)^{2}}-1\Bigr),\qquad
\rho'(x)=\frac{x}{\sqrt{1+(x/\epsilon)^{2}}},
\tag{0.4}
$$

满足 $\rho(0)=0$、$\rho$ 凸、$\rho$ Lipschitz（常数 $L=\epsilon$）、$\rho'(x)\le \epsilon$。
两个极限：

$$
\rho(x)\;\simeq\;\tfrac{x^{2}}{2}\quad (x\ll\epsilon),
\qquad
\rho(x)\;\simeq\;\epsilon\,x-\tfrac{\epsilon^{2}}{2}\quad (x\gg\epsilon).
\tag{0.5}
$$

即 $\epsilon$ 控制"二次区（处处光滑）$\to$ L1/TV 区（边缘保持、但带收缩偏置）"的过渡。
记号 $\rho'\le\epsilon$ 在 §4.2 是幅值收缩的**定量**来源。

> 这里要提前指出一处容易被忽略的事实：代码里 $\rho$ 是**逐面**作用的（见 §2.3 定理 2.3 与注 2.4），
> 因此离散泛函的连续极限在 TV 区是**方向可分离（各向异性）**的泛函 $R_{\mathcal G}$，
> 而不是教科书上的各向同性 $\int\rho(|\nabla\mathbf M|)\,\mathrm dV$。两者仅在 $\rho$ 为二次时重合。
> 这个差别在本算例上可达 $50\%$ 量级（证据 C5），论文里必须写明用的是哪一个。

### 0.3 网格与离散未知量

采用与代码一致的**非均匀树网格**（`discretize.TreeMesh`），把 $\Omega$ 剖分为单元 $c\in\mathcal C$：

- 单元 $c$：体积 $V_c$、特征线尺度 $h_c$（立方单元时 $V_c=h_c^{3}$）；
- **内部面** $f\in\mathcal F$：面积 $S_f$，两侧相邻单元中心距 $h_f$；
  **面对偶体积**（本文的核心对象）
  $$
  \boxed{\;w_f\;=\;S_f\,h_f\;}
  \tag{0.6}
  $$
- 单元常值展开：$\mathbf M_h(\mathbf r)=\sum_{c\in\mathcal C}\mathbf m_c\,\mathbf 1_c(\mathbf r)$，
  未知量 $\mathbf m=(\mathbf m_c)_{c\in\mathcal C}\in\mathbb R^{3N}$，$N=|\mathcal C|$；
- 弱形式梯度算子 $G\in\mathbb R^{n_F\times N}$（`mesh.cell_gradient`，标量）；
- 面质量矩阵 $W=\mathrm{diag}(w_f)\in\mathbb R^{n_F\times n_F}$（`get_face_inner_product` 的对角元）；
- 对有 $n$ 个相邻单元的面，$(G\mathbf m)_f$ 由这些单元值线性组合而成（(2.1) 给出两单元情形）。

离散反问题：

$$
\min_{\mathbf m\in\mathbb R^{3N}}\;\; \bigl\|A\mathbf m-\mathbf b\bigr\|_2^{2}
\;+\;\lambda\,R_h(\mathbf m),
\qquad
R_h(\mathbf m)=\sum_{f\in\mathcal F_{\text{int}}} w_f\,
\rho\Bigl(\bigl|(G\mathbf m)_f\bigr|_2\Bigr),
\tag{0.7}
$$

其中 $|(G\mathbf m)_f|_2=\bigl(\sum_{i=1}^{3}\bigl((Gm_i)_f\bigr)^{2}\bigr)^{1/2}$ 是对 $\mathbf M$ 三个分量取欧氏范数，
$\mathcal F_{\text{int}}$ 为内部面集合（与代码一致：只统计两侧单元都在磁体内的面）。

**本文要回答的两个问题**

- **(Q1)** 为什么 $R_h$ 的权重必须取 $w_f$（而不是 1，也不是 $V_c$、$V_c^{2/3}$ 之类）？
- **(Q2)** 为什么 $\lambda$ 或权重取得不当，结果会"欠正则化"，并表现为**网格密的区域偏小、稀疏的区域偏大**？

---

## §1 正问题的离散：单元体积是求积权重（后续论证的杠杆）

### 引理 1.1（列因子分解）

把 (0.1) 用 (0.6) 的单元常值展开做中值求积（取单元中心 $\mathbf x_c$）：

$$
\mathbf B_h(\mathbf r_i)=\frac{\mu_0}{4\pi}\sum_{c\in\mathcal C}V_c\,
\mathsf K(\mathbf r_i,\mathbf x_c)\,\mathbf m_c .
\tag{1.1}
$$

因此灵敏度矩阵可以**精确地**写成

$$
\boxed{\;A=\bar A\,D,\qquad D=\mathrm{diag}\bigl(V_c\mathsf I_3\bigr)_{c\in\mathcal C}\;}
\tag{1.2}
$$

其中 $\bar A$ 只含几何量（测点与单元中心的相对位置），**不含任何体积**。

*证明*：直接比较 (1.1) 与 $A\mathbf m$ 的分量形式即可；$D$ 是把每个单元的 $3$ 个分量都乘以 $V_c$ 的分块对角矩阵。$\square$

**意义**：单元体积在正问题里就是积分权重（$V_c$ 是求积系数），它出现在 $A$ 的**列**上。
这一条会在 §4.1 变成"最小 $\ell^{2}$ 范数解偏向大单元"的代数根源。

### 推论 1.2（两个范数给出相反的幅值分配）

考虑数据只约束总量（低分辨极限）的退化情形 $\bar A D\mathbf m=\mathbf b\iff \sum_c V_c\,\mathbf m_c=\mathbf C$，
比较两种"最小范数"准则：

$$
\text{(i)}\;\min\|\mathbf m\|_2^{2}\ \text{s.t.}\ \sum_cV_c\mathbf m_c=\mathbf C;
\qquad
\text{(ii)}\;\min\|\mathbf m\|_{V}^{2}\ \text{s.t.}\ \sum_cV_c\mathbf m_c=\mathbf C,
\tag{1.3}
$$

其中**体积度量**定义为

$$
\|\mathbf m\|_V^{2}:=\sum_{c}V_c\,|\mathbf m_c|^{2}
\qquad\Bigl(\text{即 }\int_\Omega|\mathbf M_h|^{2}\,\mathrm dV\text{ 的离散对应}\Bigr).
\tag{1.4}
$$

用拉格朗日乘子（$\mathcal L=\sum_c a_c|\mathbf m_c|^{2}-\mu(\sum_cV_c\mathbf m_c-\mathbf C)$，$a_c\in\{1,V_c\}$）：

$$
\text{(i)}\;2\mathbf m_c=\mu V_c\ \Rightarrow\ \mathbf m_c\propto V_c
\ \Rightarrow\
\frac{|\mathbf m|_{\text{粗}}}{|\mathbf m|_{\text{细}}}=\frac{V_{\text{粗}}}{V_{\text{细}}}=\Bigl(\frac{h_{\text{粗}}}{h_{\text{细}}}\Bigr)^{3};
\tag{1.5}
$$

$$
\text{(ii)}\;2V_c\mathbf m_c=\mu V_c\ \Rightarrow\ \mathbf m_c=\text{const}
\ \Rightarrow\
\frac{|\mathbf m|_{\text{粗}}}{|\mathbf m|_{\text{细}}}=1 .
\tag{1.6}
$$

**即：$\ell^{2}$ 准则把磁化强度堆到大单元上（比值 $=$ 体积比 $=8$，对 2:1 的加密档），
体积度量准则则给出均匀分配。** 数值核验：远场受控算例实测 $7.997$（理论 $8.0$）；
体积度量解析解 $1.000$（证据 C6、附录 D-②）。

---

## §2 为什么体积加权是对的

### 2.1 面梯度与面对偶体积

**引理 2.1（离散梯度的显式形式）**
对树网格，若面 $f$ 的两侧单元中心为 $\mathbf x_{c^+},\mathbf x_{c^-}$，法向单位矢 $\hat{\mathbf n}_f$（由 $c^-$ 指向 $c^+$），则

$$
\bigl(G\mathbf m\bigr)_f=\frac{\mathbf m_{c^+}-\mathbf m_{c^-}}{h_f},
\qquad h_f=\bigl|(\mathbf x_{c^+}-\mathbf x_{c^-})\cdot\hat{\mathbf n}_f\bigr|,
\tag{2.1}
$$

且 $W$ 的对角元恰为 $w_f=S_fh_f$。等价地，$G$ 的元素是 $\pm 1/h_f$，行和为 $0$。

*数值核验*：在非均匀（含悬挂面）树网格上逐面打印 $G$ 的系数与 $w_f$，得系数 $=\pm 1/h_f$、
由 $w_f/h_f$ 反推的面积精确等于面面积（证据 C1）。

于是对任意 $\mathbf m$：

$$
\bigl(G\mathbf m\bigr)^{\mathsf T}W\bigl(G\mathbf m\bigr)
=\sum_{f}w_f\bigl|(G\mathbf m)_f\bigr|_2^{2}
=\sum_{f}S_f\,\frac{\bigl|\Delta\mathbf m_f\bigr|^{2}}{h_f},
\qquad \Delta\mathbf m_f:=\mathbf m_{c^+}-\mathbf m_{c^-}.
\tag{2.2}
$$

### 2.2 二次情形：$w_f$ 是唯一使泛函收敛到 $\int|\nabla\mathbf M|^{2}$ 的权重

**引理 2.2（弱形式恒等式与收敛）**
设 $\mathbf m$ 是光滑 $\mathbf M$ 的单元中心取样，$\rho(x)=x^{2}/2$。则

$$
\sum_{f}w_f\bigl|(G\mathbf m)_f\bigr|_2^{2}
=\int_\Omega\bigl|\nabla\mathbf M\bigr|_F^{2}\,\mathrm dV\;+\;O(h^{2}),
\tag{2.3}
$$

其中 $|\nabla\mathbf M|_F^{2}=\sum_{i,j}(\partial_jM_i)^{2}$。

*证明*：由 (2.1) 与中值定理，
$\Delta\mathbf m_f=h_f\,\bigl(\hat{\mathbf n}_f\cdot\nabla\bigr)\mathbf M(\mathbf x_f)+O(h^{2})$，
$\mathbf x_f$ 为面心。代入 (2.2)：

$$
\sum_f S_f\frac{|\Delta\mathbf m_f|^{2}}{h_f}
=\sum_f S_fh_f\bigl|\bigl(\hat{\mathbf n}_f\cdot\nabla\bigr)\mathbf M(\mathbf x_f)\bigr|^{2}+O(h^{2}).
\tag{2.4}
$$

按面法向把 $\mathcal F$ 分成 $x,y,z$ 三族。对每族，$\sum_f S_fh_f\varphi(\mathbf x_f)$ 是体积分
$\int\varphi\,\mathrm dV$ 的**中值求积**（每个单元在自己那一族里恰好贡献一个 $S_fh_f\approx V_c$ 的
求积块，见 2.4），误差 $O(h^{2})$。三族相加即得

$$
\int\sum_{d=x,y,z}\bigl|(\hat{\mathbf e}_d\cdot\nabla)\mathbf M\bigr|^{2}\mathrm dV
=\int|\nabla\mathbf M|_F^{2}\mathrm dV . \qquad\square
$$

**数值核验（E1）**：非均匀网格上 $n_{\text{base}}=16/32/64$ 时相对误差
$5.01\%\to1.68\%\to0.46\%$，观测收敛阶 $1.57\sim1.72$（证据 C2）。

**反面（关键）**：若把权重换成 $1$，则同一被积函数变成
$\sum_f|\Delta\mathbf m_f|^{2}/h_f^{2}=\sum_f h_f^{-2}\cdot h_f^{-1}\cdot(S_fh_f)|\partial\mathbf M|^{2}$，
即每个求积块被额外乘上 $h_f^{-3}$：

$$
\sum_{f}\bigl|(G\mathbf m)_f\bigr|^{2}
\;\approx\;\sum_f \frac{1}{h_f^{3}}\,S_fh_f\,\bigl|\partial\mathbf M\bigr|^{2}
\;\longrightarrow\;\int_\Omega\frac{\bigl|\nabla\mathbf M\bigr|^{2}}{h(\mathbf x)^{3}}\,\mathrm dV .
\tag{2.5}
$$

当 $h\to0$ 时该泛函**发散**（$\propto h^{-3}$），因此它**不是任何连续泛函的离散化**。
数值上不加权的能量误差为 $338.8\%\to3467\%\to28628\%$（证据 C2），与 (2.5) 一致。

### 2.3 一般罚函数：极限泛函是"方向可分离"的那个

**定理 2.3（一致性）**
设 $\rho\in C^{1}([0,\infty))$、凸、$\rho(0)=0$、$|\rho'|\le L$。记 $\mathbf m$ 为光滑 $\mathbf M$ 的
单元中心取样。则

$$
R_h(\mathbf m)=\sum_{f}w_f\,\rho\bigl(|(G\mathbf m)_f|_2\bigr)
\;=\;R_{\mathcal G}[\mathbf M]\;+\;O(h^{\alpha}),\qquad \alpha\ge1,
\tag{2.6}
$$

$$
\boxed{\;R_{\mathcal G}[\mathbf M]:=\int_\Omega\sum_{d=x,y,z}
\rho\Bigl(\bigl|\partial_d\mathbf M\bigr|_2\Bigr)\mathrm dV ,
\qquad \bigl|\partial_d\mathbf M\bigr|_2^{2}=\sum_{i=1}^{3}\bigl(\partial_dM_i\bigr)^{2}.}
\tag{2.7}
$$

*证明骨架*：
(1) 把每个内部面归属给其**一侧单元**（例如 $+x,+y,+z$ 三个方向的面），则每个单元恰好拥有 3 个面，
$w_f=S_fh_f\approx V_c$，于是

$$
R_h=\sum_{c}\sum_{d=1}^{3}w_{f(c,d)}\,\rho\bigl(|(G\mathbf m)_{f(c,d)}|_2\bigr).
$$

(2) 由 (2.1) 与中值定理，$(G\mathbf m)_{f(c,d)}=(\hat{\mathbf e}_d\cdot\nabla)\mathbf M+O(h)$，
再用 $\rho$ 的 Lipschitz 性 $|\rho(a)-\rho(b)|\le L|a-b|$ 把误差从二次情形传到 $\rho$ 情形。

(3) 对每个 $(c,d)$ 以中值求积替换求和（同引理 2.2），得
$\sum_cV_c\sum_d\rho(|\partial_d\mathbf M|_2)+O(h^{\alpha})=R_{\mathcal G}[\mathbf M]+O(h^{\alpha})$。$\square$

**注 2.4（各向异性 vs 各向同性，必须写明）**
$R_{\mathcal G}$ 是把 $\rho$ 分别作用在**三个坐标方向**的导数上再求和（方向可分离、TV 区各向异性），
而教科书上的各向同性泛函是

$$
R_{\text{iso}}[\mathbf M]=\int_\Omega\rho\bigl(|\nabla\mathbf M|_F\bigr)\mathrm dV .
\tag{2.8}
$$

两者关系：由 $\sum_d a_d\ge\sqrt{\sum_d a_d^{2}}$（$a_d\ge0$）有 $R_{\mathcal G}\ge R_{\text{iso}}$，
一般相差一个 $O(1)$ 因子（在 TV 区最大可达 $\sqrt3$ 量级）。

- 当 $\rho(x)=x^{2}/2$（二次区，即 $|\nabla\mathbf M|\ll\epsilon$）：
  $\sum_d\frac12|\partial_d\mathbf M|^{2}=\frac12|\nabla\mathbf M|_F^{2}$，两者**重合**。
- 在深 L1 区两者**不重合**。数值核验（证据 C5）：取光滑紧支矢量场，
  $n_{\text{base}}=16/32/64$ 时离散值与 $R_{\mathcal G}$ 的相对差为 $5.30\%\to2.93\%\to2.05\%$（收敛），
  而与 $R_{\text{iso}}$ 的相对差为 $44.7\%\to48.3\%\to49.6\%$（不收敛，比值 $1.53$）；
  改用二次区（$\epsilon$ 很大）时两个参照重合（比值 $1.0000$），离散值以 $5.01\%\to1.68\%\to0.46\%$ 收敛。

**工程含义**：代码现在的写法等价于**方向可分离 TV**，这对"盒式/树网格 + 有限体积"是标准且合理的选择，
但它不是各向同性 TV。若论文要声称 $\int\rho(|\nabla\mathbf M|)\mathrm dV$，需要改成单元装配
$\sum_cV_c\rho(|\nabla\mathbf M(\mathbf x_c)|_F)$（例如先把面梯度平均回单元中心），
或明确声明所解的是 $R_{\mathcal G}$。（本文其余部分对 $R_{\mathcal G}$ 与 $R_{\text{iso}}$ 的结论完全平行，
因为一致性论证只用到"每个单元 3 个求积块、每块权重 $\approx V_c$"这一点。）

### 2.4 $W$ 的几何意义与"按体积公平"

**推论 2.5** 每个内部面被唯一归属给一个单元、每单元 3 个面、$w_f=S_fh_f\approx V_c$，故

$$
\sum_{f\in\mathcal F_{\text{int}}}w_f\;=\;3V_\Omega\bigl(1+O(h)\bigr),
\tag{2.9}
$$

并且**权重按体积分配**：某档网格单元所占磁体体积的比例，等于该档 $w_f$ 之和占 $\sum w_f$ 的比例。

*数值核验*：生产型网格（`voxel_size_base` $=6/3/1.5$ mm）实测
$\sum w_f/V_{\text{磁体}}=2.31\to2.69\to2.84\to3$；细档占磁体体积 $18\%$，
占 $\sum w_f$ 亦为 $18\%$（证据 C3、C4）。这条排除了"悬挂面把权重量算错"的怀疑。

---

## §3 不加权写法到底是什么：$O(h^{-3})$ 且**空间失真**

### 命题 3.1（局部权重）

$$
\tilde R_h(\mathbf m):=\sum_{f}\rho\bigl(|(G\mathbf m)_f|_2\bigr)
=\sum_{f}\frac{1}{h_f^{3}}\,\underbrace{w_f}_{=S_fh_f\approx h_f^{3}}\,\rho(\cdot)
\;\longrightarrow\;\int_\Omega\frac{\sum_d\rho\bigl(|\partial_d\mathbf M|_2\bigr)}{h(\mathbf x)^{3}}\,\mathrm dV .
\tag{3.1}
$$

*证明*：把 $w_f$ 拆成 $h_f^{3}\cdot(S_f/h_f^{2})$，立方单元 $S_f=h_f^{2}$ 故 $S_f/h_f^{2}=1$；
再用定理 2.3 的中值求积。$\square$

### 推论 3.2（$h^{-3}$ 律）

在均匀网格族上，$\tilde R_h=\dfrac{1}{\langle w\rangle}R_h\,(1+O(h))$，其中

$$
\langle w\rangle:=\frac{W}{n_f}=\frac{\sum_fw_f}{|\mathcal F_{\text{int}}|}=O(h^{3}).
\tag{3.2}
$$

**因此 $\tilde R_h/R_h=O(h^{-3})$：网格每加密一倍（线性尺度 $2:1$），比值涨 $2^{3}=8$ 倍。**

*数值核验*（证据 C2）：体积加权与不加权能量之比
$4.62\to36.28\to288.61$，即每加密 $7.86\times$、$7.95\times$（理论 $8$）。

### 推论 3.3（$\lambda$ 换算定理）

**严格版本（相差常数因子的权重族）**：若两套权重满足 $\tilde w_f=\kappa\,w_f$（同一个常数 $\kappa$，
例如 $w_f$、$w_f/W$、$w_f/10^{9}$ 三者之间），则 $\tilde R_h=\kappa R_h$ **对一切 $\mathbf m$ 精确成立**，
于是同一物理强度对应

$$
\boxed{\;\tilde\lambda\;=\;\frac{1}{\kappa}\,\lambda\;}
\qquad\text{特别地}\quad
w_f\to\frac{w_f}{W}:\ \tilde\lambda=W\lambda ;
\qquad
w_f\to\frac{w_f}{10^{9}}:\ \tilde\lambda=10^{9}\lambda .
\tag{3.3a}
$$

**平均版本（$1$ 与 $w_f$ 之间）**：不加权写法 $\tilde R_h=\sum_f\rho_f$ 与体积加权 $R_h=\sum_fw_f\rho_f$
一般**不成比例**——比例是 $1/\langle w\rangle_\rho$，其中 $\langle w\rangle_\rho$ 是以 $\rho_f$ 为权的
$w_f$ 加权平均。只有"所有 $w_f$ 相等"或"$\rho_f$ 在面上均匀"时才有 $\tilde R_h=R_h/\langle w\rangle$。故

$$
\tilde\lambda\;\approx\;\langle w\rangle\,\lambda,
\qquad \langle w\rangle=\frac{W}{n_f}=O(h^{3}).
\tag{3.3b}
$$

**这里的 $\approx$ 不能省略**：偏差正是推论 3.5 的档间失真（$\langle w\rangle_\rho$ 在
$[\min_fw_f,\max_fw_f]$ 内随解变化，非均匀网格上跨度可达 $8$ 倍），**调 $\lambda$ 无法消除**。
数值核验：E2a 是在 $w_f$ 与 $w_f/W$ 之间做的（严格族），实测目标函数差 $0$、解差 $5.1\times10^{-14}$（证据 C7）。

**含义**：把正则化从"不加权"改成"体积加权"时，**同一个物理强度对应的 $\lambda$ 大约要乘 $\langle w\rangle$**。
$\langle w\rangle$ 是平均面对偶体积，$\propto h^{3}$，故它随网格加密迅速变小：
生产网格（$1.0\sim1.4$ mm，$\langle w\rangle\approx5\ \mathrm{mm^{3}}$）约 $5$ 倍，
本实验 quick/full 配置分别为 $108$、$14$ 倍（证据 C7）。
**反之，若换了写法却沿用旧 $\lambda$，等效物理强度就会偏离 $\langle w\rangle$ 倍。**

### 推论 3.4（量纲与单位陷阱）

$\lambda$ 必须承担 $1/\text{体积}$ 的量纲（因为 $E_{\text{data}}$ 是 $\text{T}^{2}$ 而 $R_h$ 是
$[\rho]\times\text{体积}$）。因此

- 若把 $w_f$ 的体积单位由 $\mathrm{mm^{3}}$ 换成 $\mathrm{m^{3}}$（即 $w_f\mathbin{/}10^{9}$），
  则同一 $\lambda$ 对应的物理强度变化 $10^{9}$ 倍；
- 反过来，若代码里 $E_{\text{data}}$ 的几何用 m、而正则化权重用 $\mathrm{mm^{3}}$（当前生产代码的情形），
  就相当于在 $\lambda$ 里隐含了 $10^{9}$ 的因子。

*数值核验*（证据 C8）：$w_f$ 用 $\mathrm{m^{3}}$ 且沿用生产 $\lambda=10^{-10}$ 时，
等效强度落在 $\lambda_{bal}\times10^{-12.8}$（quick）$/\;10^{-11.5}$（full），
即深陷欠正则化区；实测还原误差 $60\%\sim200\%$、档间幅值偏置 $>1.4$。

### 推论 3.5（空间失真：非均匀网格上的**结构性**偏置）

由 (3.1)，不加权写法给每个求积块的权重是 $1$，等价于给**单位体积**的权重
$3/h(\mathbf x)^{3}$；而体积加权给的是 $3$（常数）。于是线性尺度 $2:1$ 的两档网格，
其**单位体积上的有效正则化强度相差 $2^{3}=8$ 倍**——这个因子

- 与 $\lambda$ 无关（调 $\lambda$ 只能整体缩放，不能消除档间差异）；
- 与 $|\nabla\mathbf M|$ 无关（是纯粹的网格几何事实）；
- 方向确定：**细网格单元受到的正则化更强** $\Rightarrow$ 细档被压得更平/更小。

这解释了为什么"换回不加权"不是解：它本身就是网格相关的泛函。

*数值核验*（证据 C9）：同一物理强度下，不加权写法的档间幅值偏置为 $0.157$，
体积加权为 $0.039$（$4$ 倍差距）；在同一 $\lambda$ 下把解搬到另外两张网格，
体积加权的误差只变 $1.07$ 倍，不加权变 $7.71$ 倍（证据 C11）。

---

## §4 为什么取值不当会"欠正则化"，并表现为"密小疏大"

### 4.1 变分问题的两层结构

记 $E(\mathbf m)=\|A\mathbf m-\mathbf b\|_2^{2}$，$R_h$ 如 (0.7)。$J_\lambda=E+\lambda R_h$。
数据项把解限制在**可行集**

$$
\mathcal D_\delta:=\bigl\{\mathbf m:\ \|A\mathbf m-\mathbf b\|_2\le\delta\bigr\},
\qquad \delta\simeq\sqrt{N_{\text{meas}}}\,\sigma\ \ (\sigma\ \text{为测量噪声标准差}),
\tag{4.1}
$$

之内（$\delta$ 是数据能达到的最小残差量级）。于是：

- **$\lambda$ 足够大**：解由 $R_h$ 在 $\mathcal D_\delta$ 上的极小解决定（正则化主导）；
- **$\lambda\to0$**：$R_h$ 不再起作用，解由 $\mathcal D_\delta$ 的**内部选择**决定，
  而该选择实际上由**算法的隐式范数**给出——$A$ 秩亏/病态时，从 $\mathbf m=\mathbf 0$ 出发的
  梯度型迭代（L-BFGS-B、LSQR、共轭梯度）收敛到**最小 $\ell^{2}$ 范数解**
  $\mathbf m^{\ell2}=A^{+}\mathbf b=\lim_{\mu\to0}(A^{\mathsf T}A+\mu I)^{-1}A^{\mathsf T}\mathbf b$
  （迭代向量都是 $A^{\mathsf T}A$ 的多项式作用在 $A^{\mathsf T}\mathbf b$ 上，其公共极限即该式）。

（当前代码恰是 $M_x=\mathbf 0$ 起步的 L-BFGS-B，所以 $\lambda$ 太小时落到 $\mathbf m^{\ell2}$。）

### 4.2 核心定理：最小 $\ell^{2}$ 范数解偏向大单元

**定理 4.1（$\ell^{2}$ 偏置定理）**
设 $A=\bar AD$（引理 1.1）。在数据不约束的（近）零空间方向上，

$$
\mathbf m^{\ell2}=\arg\min\bigl\{\|\mathbf m\|_2:\ \bar AD\mathbf m=\mathbf b\bigr\},
\tag{4.2}
$$

它等价于在**矩变量** $\mathbf q:=D\mathbf m$（$q_c=V_c\mathbf m_c$）下解

$$
\arg\min\Bigl\{\sum_c\bigl(q_c/V_c\bigr)^{2}:\ \bar A\mathbf q=\mathbf b\Bigr\}
\;=\;\arg\min\bigl\{\|D^{-1}\mathbf q\|_2:\ \bar A\mathbf q=\mathbf b\bigr\},
\tag{4.3}
$$

即每单元的"矩"被 $1/V_c^{2}$ 加权。可见 $\ell^{2}$ 准则**倾向于把矩放进 $V_c$ 大的单元**。
在只约束总量的退化情形（推论 1.2）可精确求解：

$$
\text{(i)}\ \ell^{2}:\ q_c\propto V_c^{2}\ \Rightarrow\ \mathbf m^{\ell2}_c=\frac{q_c}{V_c}\propto V_c ,
\qquad
\text{(ii)}\ V\text{-度量}:\ q_c\propto V_c\ \Rightarrow\ \mathbf m^{V}_c=\text{const}.
\tag{4.4}
$$

于是同一真值在两档网格上的幅值比分别为

$$
\frac{|\mathbf m|_{\text{粗}}}{|\mathbf m|_{\text{细}}}\Big|_{\ell^{2}}=\Bigl(\frac{h_{\text{粗}}}{h_{\text{细}}}\Bigr)^{3}=8,
\qquad
\frac{|\mathbf m|_{\text{粗}}}{|\mathbf m|_{\text{细}}}\Big|_{V}=1 .
\tag{4.5}
$$

*数值核验*（证据 C6、附录 D-②）：把测点放到 $200$ mm 远处（几何核 $\mathsf K$ 在磁体上近似常数，
(4.3) 的前提成立）后做**纯数据拟合**，实测幅值比 $=7.997$（理论 $8.0$）；
同一算例把测点放回正常距离（$24$ mm）时几何部分约束了该方向，比值降到 $1.716$，
**方向不变、量级被削弱**。体积度量的解析解为 $1.000$。

**注 4.2（定理的适用范围）** (4.5) 是**低分辨极限/上界**陈述，不是等式断言：
真实几何使 $A$ 的列不完全平行，偏置被削弱但不会反转（实测 $1.716\in(1,8)$）。
它给出的可检验预言是"**偏置方向恒为粗档偏高、细档偏低，且随网格对比度按体积比放大**"。

### 4.3 另一侧：$\lambda$ 太大会幅值收缩

**命题 4.3（摩擦上界）**
由 (0.4)，$\rho'(x)=x/\sqrt{1+(x/\epsilon)^{2}}\to\epsilon\,\mathrm{sgn}(x)$（$|x|\gg\epsilon$），
故深 L1 区每个面的**正则化"摩擦"幅值被 $ \lambda\epsilon w_f$ 封顶**，与 $|\mathbf m|$ 大小无关。
最优性条件

$$
\mathbf 0=2A^{\mathsf T}(A\mathbf m-\mathbf b)+\lambda G^{\mathsf T}W\rho'(G\mathbf m)
\tag{4.6}
$$

说明：**数据力必须超过该恒定的摩擦上限，幅值才能继续增长**；若 $\lambda\epsilon$ 太大，
解在到达真值幅值之前就被"刹住"，表现为幅值整体收缩（偏置单调低于 1）。
这与 (4.5) 的粗档偏高叠加，给出一个**双侧**约束：$\lambda$ 有下界（压住 $\ell^{2}$ 偏置）
与上界（不要触发收缩）。

### 4.4 $\lambda_{bal}$：把"取值不当"变成可计算的门槛

定义**力平衡标定**

$$
\boxed{\;\lambda_{bal}=\frac{\bigl\|\nabla E(\mathbf M_{\text{ref}})\bigr\|_2}
{\bigl\|\nabla R_h(\mathbf M_{\text{ref}})\bigr\|_2\big|_{\lambda=1}}\;}
\tag{4.7}
$$

（$\mathbf M_{\text{ref}}$ 为参考模型，量级正确即可；代码实现见 `b2m_core.calibrate_lambda`）。
它的意义是"数据拉力 $=$ 正则化摩擦"的临界强度。**关键性质是协变性**：$E$ 与 $R_h$ 各自都是
一致离散（§2），因此该比值在网格加密与单位变换下不变——这就是为什么 $\lambda_{bal}$
可以跨网格搬运，而"绝对数值"不能。

由 §4.2 与 §4.3，可用区间是

$$
\lambda\in\bigl[\lambda_{bal},\,10\,\lambda_{bal}\bigr]
\qquad\Longleftrightarrow\qquad k=\log_{10}(\lambda/\lambda_{bal})\in[0,+1],
\tag{4.8}
$$

下界由 $\ell^{2}$ 竞争给出，上界由收缩给出。**"宁大勿小"** 是 (4.5) 的直接后果：
低于下界时 $\ell^{2}$ 偏置按体积比放大，而高于上界只是整体收缩（更容易被识别和补救）。

### 4.5 综合：为什么"取值不当"就表现为"密小疏大"

把 §3 与 §4.2 拼起来：

| 情形 | 正则化项给单位体积的权重 | 零空间方向由谁决定 | 结果 |
|---|---|---|---|
| 体积加权 + $\lambda\gtrsim\lambda_{bal}$ | $3$（均匀） | 正则化项（体积公平） | 各档无偏 $\lvert\text{bias}-1\rvert\le0.05$ |
| 体积加权 + $\lambda\ll\lambda_{bal}$ | $3$（均匀） | $\ell^{2}$ 隐式范数 | 粗档偏高、细档偏低（定理 4.1） |
| 不加权 + 任意 $\lambda$ | $3/h(\mathbf x)^{3}$（**差 8 倍**） | 两项都偏向粗档（**叠加**） | 偏置比体积加权大数倍 |

**这解释了全部数值现象**：
(a) $\lambda$ 太小时两种写法都出现"密小疏大"（定理 4.1，与权重无关）；
(b) 同一物理强度下体积加权的偏置只有不加权的 $1/4$（推论 3.5：不加权的 8 倍档间权重与
$\ell^{2}$ 偏置方向相同，两者叠加）；
(c) $\lambda\gtrsim\lambda_{bal}$ 后偏置消失、且跨网格几乎不变（§2 的一致性 + 协变性）。

---

## §5 工程化推论

1. **归一化写法**（推荐）：$\displaystyle \hat R_h(\mathbf m)=\sum_f\frac{w_f}{W}\rho\bigl(|(G\mathbf m)_f|\bigr)$，
   $W=\sum_fw_f\approx3V_\Omega$。由 (2.9)，$W$ 与 $w_f$ 同量纲，故 $\hat R_h$ 对**长度单位不变**；
   由 (3.3)，$\hat\lambda=W\lambda$ 与**网格加密无关**。代码：
   ```python
   w = V_face; W = w.sum()
   loss_reg = lambda_reg * eps**2 * np.sum((w/W) * (np.sqrt(1+(gn/eps)**2) - 1))
   grad_reg = lambda_reg * (G_sub.T @ ((w/W)[:,None] * g * wt[:,None])).ravel()
   ```
2. **$\lambda$ 的定法**：先用 (4.7) 标定 $\lambda_{bal}$（参考模型只需量级正确），
   再取满足"残差 $\approx\sigma$"的**最大** $\lambda$（discrepancy principle，
   即 $\|\mathbf A\mathbf m-\mathbf b\|\lesssim\delta$，(4.1)）；两者一致性由 E5 实测确认
   （自动 $\lambda$ 与上帝视角最优 $\lambda$ 同为 $k=+1$，误差同为 $0.29\%$）。
3. **定理 ↔ 数值判据对照**

| 结论 | 定理/推论 | 数值证据 |
|---|---|---|
| 面对偶体积 $w_f=S_fh_f$、$G$ 系数 $\pm1/h$ | 引理 2.1 | C1 |
| 体积加权 $\to\int\lvert\nabla\mathbf M\rvert^{2}$，$O(h^{2})$ | 引理 2.2 | C2（E1） |
| 不加权 $\propto h^{-3}$、不是任何连续泛函的离散 | (2.5)、推论 3.2 | C2（$7.86\times/7.95\times$） |
| 权重按体积公平、$\sum w_f\approx3V_\Omega$ | 推论 2.5 | C3、C4 |
| 极限泛函是方向可分离的 $R_{\mathcal G}$ | 定理 2.3、注 2.4 | C5 |
| $\lambda$ 换算 $\tilde\lambda=\langle w\rangle\lambda$ | 推论 3.3 | C7（E2a，目标函数差 $0$） |
| $10^{9}$ 单位陷阱 | 推论 3.4 | C8（E2a） |
| $\ell^{2}$ 偏置 $\propto V_c$（$8:1$） | 定理 4.1 | C6（$7.997$） |
| 档间偏置随权重写法放大 | 推论 3.5 | C9（E2b，$0.157$ vs $0.039$） |
| 网格无关性 | §2 + 协变性 | C11（E4，$7.71\times$ vs $1.07\times$） |
| $\lambda$ 窗口 $[0,+1]$ | (4.8) | C10（E2a/E2b/E5） |

---

## §6 假设、边界与未证之处

1. **连续反问题本身不适定**，与离散化无关；本文只讨论"离散泛函是否忠实、$\lambda$ 是否可比"，
   不声称能恢复不可分辨的分量。
2. §2 的收敛结论假设：$\rho$ 凸、Lipschitz、$\rho(0)=0$；$\mathbf M$ 光滑且（或其法向导数）
   在 $\partial\Omega$ 上可忽略——代码只累加**内部面**，边界面的缺失带来 $O(h)$ 项；
   $\mathbf m$ 取单元中心值（中点求积），误差 $O(h^{2})$（严格依据见附录 A）。
3. **悬挂面**：`discretize` 对"粗单元对面 4 个细单元"的界面取 $S_f$ 为粗面面积、$h_f$ 为中心距；
   该约定使 $\sum w_f/V\to3$、二次泛函以 $O(h^{2})$ 收敛（C2/C3 已核验），
   但本文未逐面证明它对任意 $\rho$ 都是最优装配。
4. **定理 4.1 是低分辨理想化极限**（列方向近似平行、单约束）。真实几何是"削弱版"：
   实测 $1.716$ 落在 $(1,8)$ 内。它给出的是**方向性预言与上界**，不是等式。
5. **点偶极近似**（(1.1) 的中值求积）要求单元尺度 $\ll$ 测量距离；该误差与本文讨论的
   正则化标度问题独立，但会污染"还原度"的绝对值。
6. **窗口 $[0,+1]$ 只给出量级**，不是精确常数；上下界依赖于噪声水平、网格对比度与
   $\rho$ 的 $\epsilon$ 取值。$\lambda_{bal}$ 对参考模型幅值有弱依赖（深 L1 区近似线性、
   二次区近似无关），由 $10$ 倍窗口吸收。
7. **未覆盖**：各向同性 TV 的正确装配（注 2.4 给出了改法但未实现/未验证）、
   $L^{0}$/稀疏先验、$\epsilon$ 的最优选取、悬挂面外的其它非拟一致网格。
8. 单元中心取样 vs 面心取样的差异（(2.4) 的一阶项）在**均匀**网格上严格抵消，
   在非均匀网格上给出 $O(h)$ 项；实测收敛阶 $1.57\sim1.72$ 与其一致。

---

## 附录 A：用到 FEM 标准结果的地方

1. **$W=\mathrm{diag}(S_fh_f)$ 是质量集中（mass lumping）后的面质量矩阵。**
   对 $H(\mathrm{div})$ 相容的有限体积/混合有限元离散，弱形式 Dirichlet 能量
   $\int|\nabla M|^{2}\mathrm dV$ 的离散对应是 $(GM)^{\mathsf T}W(GM)$，
   其中 $W$ 为面对偶体积（$=$ 面面积 $\times$ 法向中心距）；对正交网格该矩阵对角，
   对角元即 $S_fh_f$。本文 (2.2) 只是把这个恒等式写成显式求和。
2. **插值/求积误差 $O(h^{2})$。** Ciarlet（*The Finite Element Method for Elliptic Problems*）的插值误差估计给出
   $\|\nabla M-\nabla M_h\|_{L^{2}}\le Ch|M|_{H^{2}}$；Strang & Fix（*An Analysis of the Finite Element
   Method*）的 Strang 第一引理把"变分犯罪"（variational crime）的误差
   $|a_h(M_h,v_h)-a(M,v_h)|$ 界住，从而能量误差为 $O(h^{2})$——
   这正是 (2.3) 与实测 $5.01\%\to1.68\%\to0.46\%$（比值约 $4$）的来源。
3. **质量集中（mass lumping）的收敛性**：对角化的面/单元质量矩阵保持 $O(h^{2})$ 能量精度
   （对正交网格的常系数问题）；否则会退化到 $O(h)$。本文 (2.9) 的 $\sum w_f\to3V_\Omega$
   是同一事实的几何推论。
4. **最小范数解与 Krylov 迭代**：从 $\mathbf 0$ 出发的梯度型方法在一致线性系统上的极限为
   $A^{+}\mathbf b=\lim_{\mu\to0}(A^{\mathsf T}A+\mu I)^{-1}A^{\mathsf T}\mathbf b$
   （Tikhonov 解族的极限），即最小 $\ell^{2}$ 范数解——(4.2) 的前提。

## 附录 B：符号 ↔ 代码对照

| 本文 | 代码（`B2M_newprob.py` 第 5 节 / `meshtest/b2m_core.py`） | 含义 |
|---|---|---|
| $\mathbf m$ | `M_sol_current.reshape(n_voxels, 3)` | 单元常值磁化强度 |
| $G$ | `mesh.cell_gradient` → `G_sub` | 弱形式梯度（元素 $\pm1/h_f$） |
| $W=\mathrm{diag}(w_f)$ | `mesh.get_face_inner_product().diagonal()` → `V_all` | 面质量矩阵 / 面对偶体积 |
| $w_f=S_fh_f$ | `V_face = V_all[internal_faces]` | 面权重（$\mathrm{mm^{3}}$） |
| $R_h$ | `loss_reg` | 体积加权正则化 |
| $\tilde R_h$ | 把 `V_face` 换成 `1` 的版本 | 不加权写法 |
| $\hat R_h$ | `V_face / W` | 归一化写法（推荐） |
| $\rho,\ \epsilon$ | `huber_epsilon**2*(sqrt(1+(gn/eps)**2)-1)` | Charbonnier/Huber |
| $\rho'/x$ 权重 | `weights = 1/sqrt(1+(gn/eps)**2)` | 梯度里的同一权重 |
| $\lambda$ | `lambda_reg`（当前 $10^{-10}$，$w_f$ 用 $\mathrm{mm^{3}}$） | 正则化强度 |
| $\lambda_{bal}$ | `b2m_core.calibrate_lambda` | 力平衡标定值 |
| $\delta,\sigma$ | `results.md` 中的 `σ`、`残差` | 噪声水平与残差 |

## 附录 C：数值证据汇总

（E1–E5 均来自 `meshtest/results.md`；C1、C5、C6 为本文写作时的受控核验，命令见附录 D。）

| 编号 | 内容 | 数值 |
|---|---|---|
| C1 | $G$ 系数与 $w_f$ | 系数 $=\pm1/h_f$；$w_f/h_f$ 反推面积 $=$ 面面积 |
| C2 | 体积加权 / 不加权 梯度能量收敛 | 加权 $5.01\%\to1.68\%\to0.46\%$；不加权 $338.8\%\to3467\%\to28628\%$；两者比值 $4.62\to36.28\to288.61$（每加密 $7.86\times$、$7.95\times$） |
| C3 | $\sum w_f/V_{\text{磁体}}$ | $2.31\to2.69\to2.84$（生产型网格，$6/3/1.5$ mm） |
| C4 | 权重按体积分配 | 细档占磁体体积 $18\%$、占 $\sum w_f$ $18\%$ |
| C5 | 极限泛函是方向可分离的 $R_{\mathcal G}$ | 深 L1：与 $R_{\mathcal G}$ 差 $5.30\%\to2.93\%\to2.05\%$；与 $R_{\text{iso}}$ 差 $44.7\%\to48.3\%\to49.6\%$（参照比 $1.53$）。二次区：两参照重合（比 $1.0000$），离散 $5.01\%\to1.68\%\to0.46\%$ |
| C6 | $\ell^{2}$ 偏置 $\propto V_c$ | 远场（$200$ mm）纯数据拟合 粗/细 $=7.997$（理论 $8.0$）；正常距离（$24$ mm）$=1.716$；$V$-度量解析 $=1.000$ |
| C7 | $\lambda$ 换算 | $\langle w\rangle=W/n_f$：quick $108.42$、full $13.93$、生产 $\approx5$；等效 $\lambda$ 下两约定目标函数差 $0$（$10^{-16}$）、解差 $5.1\times10^{-14}$ |
| C8 | $\mathrm{m^{3}}$ 陷阱 | 落到 $k=-12.8$（quick）$/-11.5$（full）；误差 $60\%\sim200\%$、档间偏置 $>1.4$ |
| C9 | 档间偏置（指纹起始 $\lambda$） | 不加权 $0.157$ vs 体积加权 $0.039$（比值 $0.25$）；更深处 $1.97$ vs $1.65$ |
| C10 | $\lambda$ 窗口 | quick 安全 $k=[-2,+1]$、full $k=[0,+3]$；判据取交集 $[0,+1]$；窗口内误差 $\le2\%$（full）、偏置 $\le0.05$ |
| C11 | 网格无关性 | 同一 $\lambda$ 搬到三张网格：体积加权 $1.07\times$、不加权 $7.71\times$ |
| C12 | $\lambda$ 选择规则 | discrepancy principle 选出的 $k=+1$ 与最优 $k=+1$ 相同、误差同为 $0.29\%$ |

## 附录 D：复现命令

解释器：`E:\Python\Miniforge\envs\discretize\python.exe`；工作目录 `meshtest/`。
（E1–E5：`python run_experiments.py --quick` 或 `--full`，结果见 `results.md`。）

**D-① $h^{-3}$ 律与 $\langle w\rangle$**（对应 C2/C3/C7）

```bash
python run_experiments.py --only E1     # 体积加权/不加权能量与收敛阶
python - <<'PY'
import numpy as np, sys; sys.path.insert(0,'.')
from run_experiments import setup, default_geom
domain, mbox = default_geom()
st = setup(domain, mbox, "split", (16,16,4), "volume_raw", "dense", "uniform", verbose=False)
W, nf = st["reg"]["W"], st["reg"]["n_faces"]
print(f"Σw_f={W:.4f} mm³  <w>=W/n_f={W/nf:.4f} mm³  Σw_f/V_magnet={W/st['vol'].sum():.4f}")
E_u = {16: 87.354, 32: 716.865, 64: 5776.035}; E_w = {16: 18.913, 32: 19.758, 64: 20.013}
r = {n: E_u[n]/E_w[n] for n in E_u}
print("不加权/加权比值:", {n: round(r[n],2) for n in r},
      "每加密倍率:", round(r[32]/r[16],2), round(r[64]/r[32],2))   # 期望 ≈ 8
PY
```

**D-③ 面梯度系数、面对偶体积与极限泛函**（对应 C1/C5）

```bash
python - <<'PY'
import numpy as np, sys; sys.path.insert(0,'.')
from run_experiments import setup, default_geom
from discretize import TreeMesh
st = setup(*default_geom(), "split", (16,16,4), "volume_raw", "dense", "uniform", verbose=False)
mesh = st["mesh"]; cc = mesh.cell_centers
G = mesh.cell_gradient.tocsr(); w_all = mesh.get_face_inner_product().diagonal()
for f in range(G.shape[0]):                       # C1: G 系数 = ±1/h, w_f/h = 面面积
    c = G.indices[G.indptr[f]:G.indptr[f+1]]; v = G.data[G.indptr[f]:G.indptr[f+1]]
    if len(c) != 2: continue
    d = cc[c[1]] - cc[c[0]]; h = abs(d[int(np.argmax(np.abs(d)))])
    print(f"G系数={v}  1/h={1/h:.5f}  w_f={w_all[f]:.3f}  反推S={w_all[f]/h:.3f}"); break
L = 4.0
def Mfun(X,Y,Z):
    r2 = X**2+Y**2+Z**2; w = np.where(r2 < 9.0, (1.0-r2/9.0)**2, 0.0)
    return np.stack([w, w*0.5, w*0.25], axis=-1)
def grads(X,Y,Z):
    r2 = X**2+Y**2+Z**2
    f = lambda c: np.where(r2 < 9.0, -4*c*(1.0-r2/9.0)/9.0, 0.0)
    g = np.stack([f(X), f(Y), f(Z)], axis=-1)
    return np.stack([g*1.0, g*0.5, g*0.25], axis=-2)     # [n, 3方向, 3分量]
rho = lambda x, eps: eps**2*(np.sqrt(1.0+(x/eps)**2)-1.0)
for eps in (1e-6, 1e3):                                  # C5: 深 L1 与二次区
    for nbase in (16, 32, 64):
        h = 2*L/nbase
        m = TreeMesh([[(h,nbase)],[(h,nbase)],[(h,nbase)]], origin=[-L,-L,-L], diagonal_balance=False)
        ml = m.max_level
        m.refine_bounding_box(np.array([[-L,-L,-L],[L,L,L]]), level=ml-1, finalize=False)
        m.refine(lambda c: ml if c.center[0] < 0 else ml-1, finalize=False)
        m.finalize(); m.number()
        c2 = m.cell_centers; V = m.cell_volumes
        Gs = m.cell_gradient.tocsr(); wf = m.get_face_inner_product().diagonal()
        g = Gs @ Mfun(c2[:,0], c2[:,1], c2[:,2]); it = np.diff(Gs.indptr) >= 2
        A = float(np.sum(wf[it]*rho(np.sqrt((g[it]**2).sum(1)), eps)))
        gr = grads(c2[:,0], c2[:,1], c2[:,2])
        B = float(np.sum(V*np.sum(rho(np.sqrt((gr**2).sum(-1)), eps), axis=-1)))   # R_G 方向可分离
        C = float(np.sum(V*rho(np.sqrt((gr**2).sum(-1).sum(-1)), eps)))           # R_iso
        print(f"eps={eps:g} nbase={nbase:3d} 离散={A:.6f} R_G={B:.6f} 误差={abs(A-B)/B*100:5.2f}%"
              f"  R_iso={C:.6f} 误差={abs(A-C)/C*100:5.2f}%  参照比={B/C:.4f}")
PY
```

**D-② 最小 $\ell^{2}$ 范数偏置 $\propto V_c$**（对应 C6）

```bash
python - <<'PY'
import numpy as np, sys; sys.path.insert(0,'.')
from run_experiments import setup, invert, default_geom
domain, mbox = default_geom()
# 测点放到 200mm 远处 ⇒ 几何核近似常数 ⇒ 定理 4.1 的前提成立
for zoff in (200.0, 24.0):
    st = setup(domain, mbox, "split", (16,16,4), "volume_raw", "oneside", "uniform",
               step=8.0, noise_frac=0.0, z_off=zoff, verbose=False)
    lv, vol = st["level"], st["vol"]
    M, _ = invert(st, st["lam_bal"]*1e-40, maxiter=1200)      # 几乎无正则化 = 纯数据拟合
    mm = np.linalg.norm(M, axis=1)
    r = {int(L): float(mm[lv==L].mean()) for L in np.unique(lv)}
    print(f"z_off={zoff:6.1f} 粗/细={r[1]/r[0]:6.3f}   (min-ℓ² 理论 => V比 = 8.0)")
# 解析对照：min ΣV_c m_c² 的解
v = vol; print("V-度量解析 粗/细 =", 1.000, " (=1，体积度量下无偏)")
PY
```
