# 数据集与结果

## 1. 测量数据集

三套数据，CSV 格式统一为 `x, y, z, Bx, By, Bz`（坐标 mm、磁场 µT）：

| 目录 | 文件 | 测点数 | 说明 | 谁在用 |
|---|---|---|---|---|
| `data/raw/specialU/` | `N_leg.csv` / `S_leg.csv` / `whole.csv` | 3696 / 2816 / 720 | 早期"特殊 U"磁铁，两腿分开扫 | `notebooks/B2M.ipynb` |
| `data/raw/UshapeNormal/` | `Large.csv` / `legs.csv` | 4200 / 7056 | 圆弧区 4 mm 间距、磁腿区 2 mm 间距 | `notebooks/B2M_noreg.py`、`notebooks/plot_2d_sections.py` |
| `data/raw/UshapeNormal_New_prob/` | `magnetic_field_{N,S,front,up,inside,curve}.csv` | 3968 / 4960 / 15360 / 13888 / 8432 / 2176，**合计 48784** | **小探头**重测：测点跨越了磁铁对称面，$z\in[0,40]$ mm | `notebooks/B2M_newprob.py`、`src/symmetry.py` |
| `data/raw/legacy/` | `U_normal_magnetic_field.csv`、`U_special_B.csv` | 975 / 975 | 更早的两组测量，仅存档 | — |

---

## 2. 磁铁几何

三套数据对应三块不同的磁铁（或同一块磁铁的不同测量），几何参数**不通用**：

| 参数（mm） | `specialU` | `UshapeNormal` | `UshapeNormal_New_prob` |
|---|---|---|---|
| $x$ 范围 | 104 – 186 | 111.5 – 193.5 | 96 – 186 ⚠️ |
| $y$ 范围 | 134 – 240 | 129.5 – 231.5 | 120 – 250 ⚠️ |
| $z$ 范围 | −28 – 0 | −35 – −3 | −8 – 24 ⚠️ |
| 外半径 $R_{\text{out}}$ | 41 | 41 | 45 ⚠️ |
| 内半径 $R_{\text{in}}$ | 22 | 22 | 18 ⚠️ |
| 厚度 $T$ | 28 | 32 | 32 ⚠️ |
| 基础体素 | 1.5 | 1.5 | 1.5 |
| 加密分界 $y_{\text{split}}$ | 144 | 144 | 144 |

几何模型：$U$ 形 = 上半圆环（$R_{\text{in}}\le r\le R_{\text{out}}$，圆心 $(X_c,Y_c)$）
+ 下半两条直腿（$y\le Y_c$ 的两条竖条），$Y_c=Y_{\max}-R_{\text{out}}$。

!!! danger "新数据集的几何是占位值（⚠️）"
    `UshapeNormal_New_prob` 的几何**没有实测**，是从扫描范围反推的粗估：

    - `magnetic_field_inside.csv` 覆盖 $x\in[125,157]$ ⇒ 槽宽 ≈33 mm、$X_c\approx141$；
    - `N`/`S` 扫描分别在 $x\le94$ / $x\ge189$ ⇒ 外侧面在 96 / 186 附近；
    - `front` 扫描 $y\le117$ ⇒ 磁腿下端面 $y\approx119$；
    - **$z$ 厚度 $T$ 无法从外部场反演**（$z$ 向均匀分布时外部场只依赖 $M\!\cdot\!T$），
      必须用卡尺实测。

    代码里 `GEOMETRY_IS_PLACEHOLDER = True` 会在每次运行时大声警告。
    填充实测尺寸后请把它置为 `False`。
    **本数据集的反演绝对幅值暂时不要引用**，但"对称面 $h=8$ mm""重合 37.7%"
    以及全部方法学结论（`meshtest/`）都不依赖这组占位尺寸。

---

## 3. 反演结果

### 3.1 测量磁场（`UshapeNormal`，扫描平面 $z=0$）

`Large.csv`（圆弧区，4 mm 间距，4200 点）与 `legs.csv`（磁腿区，2 mm 间距，7056 点）
合起来覆盖 U 形磁铁上方的一个平面。测量磁场的模值分布里 U 形轮廓清晰可辨：
两条磁腿外侧场强最高、槽内迅速衰减、圆弧顶部居中。

**这正是反演能成立的原因** —— 外部场的空间结构直接编码了内部磁化分布；
如果测点离得太远，几何核退化成常数，数据就只剩下"总磁矩"这一个信息量。

### 3.2 反演出的磁化强度（磁铁中面）

正则化反演给出的 $\lvert\mathbf M\rvert$ 在中面（$z\approx-19$ mm）上的分布，
叠加磁铁几何轮廓后可以看到磁化强度集中在 U 形的两条腿与圆弧上，槽内接近 0；
幅值在 $10^5$ A/m 量级，与钕铁硼永磁体的典型剩磁同量级。

### 3.3 正则化到底在做什么

同一批数据、同一个正演算子，去掉正则化项的纯最小二乘解对数据的拟合同样很好
（甚至更好），但幅值在空间上剧烈振荡、相邻单元符号乱跳 ——
这是"解不唯一，必须靠挑解规则"的直接证据，也是 `B2M_noreg.py` 存在的意义：
它专门用来展示"没有正则化会发生什么"。

### 3.4 配图由人工维护

本仓库**不做图片同步**。要往站点加图，把文件放进 `web/docs/`，
然后在页面里用 `![说明](文件名)` 引用即可（同目录，不需要 `../`）。
下面几张是现成的候选，源文件在 `figures/`（`.gitignore`，由脚本生成）：

| 建议文件名 | 内容 | 生成方式 |
|---|---|---|
| `B_field_section.png` | 测量磁场在 $z=0$ 平面的分布 | `python notebooks/plot_2d_sections.py` |
| `M_section.png` | 反演磁化强度在磁铁中面的分布 | 同上 |
| `M_Uspecial.png` | `specialU` 数据集的反演结果 | `notebooks/plot_2d_sections.ipynb` |
| `M_regularized_vs_noreg.png` | 正则化 vs 无正则化对照 | `python notebooks/B2M_noreg.py` |
| `symmetry_plane_detection.png` | $J(h)$ 凹坑与不确定带宽 | `python src/symmetry.py data/raw/UshapeNormal_New_prob` |
| `mirror_overlap.png` | 镜像重合点的分布与不一致量 | 同上 |
| `E1_discretization.png`、`E2_lambda_sweep.png`、`E5_lambda_rule.png` | 实验曲线 | `cd meshtest && python run_experiments.py --full` |

### 3.5 新数据集（小探头）目前的状态

| 项目 | 结果 |
|---|---|
| 对称面自动检测 | $h=8.000$ mm（原始最优点 7.698 mm，带宽 [7.55, 7.90] mm，已吸附到实测 $z$ 网格） |
| 数据装配 | 实测 48784 → 镜像补点 30379 → **最终 79163 点**（重合的 18405 个镜像副本被丢弃） |
| 几何一致性校验 | 任何测点（含镜像补点）都不允许落进磁铁实体内部，违反了直接报错 |
| 反演本身 | 流程跑通，但因为几何是占位值，**结果暂不发布** |

---

## 4. 结果文件在哪

`data/derived/` 里放的是**可由原始数据 + 脚本重新生成**的中间产物，
**不纳入 git**（见根目录 `.gitignore`），换台机器 clone 下来不会有它们：

```
data/derived/
├── specialU/                 Large.vtk / legs.vtk（csv2vtk 转出的结构化点云）
│                             inverted_M_results_{newmesh,newoperator}.vtu
├── UshapeNormal/             同上 + inverted_M_results_noreg.vtu（无正则化对照）
├── UshapeNormal_New_prob/    inverted_M_results_newmesh.vtu（⚠️ 占位几何）
│                             mag_structured_*.vtk（6 个扫描文件的场分布）
└── legacy/                   早期结果，仅存档
```

`.vtu` 是 `pyvista`/VTK 的 `UnstructuredGrid`，单元数据里带 `Mx, My, Mz`，
用 ParaView 或 `pyvista.read()` 直接打开。

---

## 5. 局限

1. **点偶极近似**：每个体素当点偶极子。测点离磁体表面太近（< 4~6 个单元尺度）时，
   近似误差会污染反演（`meshtest` 里踩过"真值处残差 = 90σ"的坑）。
2. **$z$ 向不可分辨**：若磁化沿 $z$ 均匀，外部场只依赖 $M\cdot T$，$M$ 与 $T$
   无法同时反演 —— 厚度必须实测。
3. **无真值**：真实磁铁没有"标准答案"，所以方法学结论全部建立在
   `meshtest/` 的合成真值对照实验上（独立 1 mm 细网格正演，避免 inverse crime）。
4. **未做**：各向异性 TV、L0/L1 稀疏先验、以及真实 U 形几何 + 真实点云上的
   $\lambda$ 标定值。
