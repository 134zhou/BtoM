# %% [markdown]
# # B2M（新版小探头数据 UshapeNormal_New_prob）
#
# 本 notebook 由 `B2M.ipynb` 改写而来，**只改"数据装配"这一块**：
#
# 1. **对称面高度 h 自动检测**（`src/symmetry.py`）：不再手工给 `Z_center = (Z_MIN+Z_MAX)/2`，
#    而是用测量数据自身的镜像一致性自动定出对称面高度，并给出不确定度与接受/拒绝判据。
# 2. **镜像重合处理**：小探头已经测到跨越对称面的数据，直接镜像会有 37.7% 的点与实测点
#    **坐标重合**。默认 `keep_measured`：实测值永远优先，重合处的镜像副本直接丢弃，
#    镜像只补进实测未覆盖的区域。
# 3. **z 边界用检测到的 h 重新锚定**：`Z_MIN/MAX = h ∓ T/2`（T 为卡尺实测厚度）。
# 4. **几何一致性校验**：任何测点（含镜像补点）都不允许落进磁铁实体内部。
#
# 敏感度矩阵（LinearOperator）、Huber 正则化、L-BFGS-B 求解、VTK 导出与原版完全一致。
#
# ## 关于镜像的符号约定
#
# 磁化强度位于 x-y 面内、且沿 z 均匀时，绕厚度中面 $z=h$ 有
#
# $$B_x(x,y,2h-z)=B_x(x,y,z),\quad B_y(x,y,2h-z)=B_y(x,y,z),\quad B_z(x,y,2h-z)=-B_z(x,y,z)$$
#
# 即 `Bz` 反号、`Bx/By` 不变 —— 与本项目原代码的约定一致（推导：把源按 $\pm z$ 成对配对，
# $M_z=0\Rightarrow \mathbf{M}\cdot\mathbf{dr}$ 与 $z$ 无关，故 $B_x/B_y$ 偶、$B_z$ 奇）。
# 注意这与"电流分布镜像不变"的情形符号相反，不能凭直觉写。

# %% [markdown]
# # 1. 导入要用到的库

# %%
import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr, LinearOperator
from scipy.optimize import minimize
import discretize
from discretize import TreeMesh
from numba import njit, prange
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# 本仓库根目录的公用模块（对称面检测 / 镜像去重）
# --- 仓库路径引导 -----------------------------------------------------------
# 向上找到含 .git 的仓库根，之后所有数据/输出路径都相对它解析，
# 这样在仓库里任何位置（notebook 目录、仓库根）运行都能找到同一批文件。
import sys as _sys, pathlib as _pl
_p = _pl.Path(globals()["__file__"]).resolve().parent if "__file__" in globals() else _pl.Path.cwd()
while not (_p / ".git").exists() and _p.parent != _p:
    _p = _p.parent
REPO = _p
RAW, DERIVED, FIGS = REPO / "data" / "raw", REPO / "data" / "derived", REPO / "figures"
_sys.path.insert(0, str(REPO / "src"))
# ---------------------------------------------------------------------------
from symmetry import (load_measurements, estimate_symmetry_plane, apply_mirror,
                      plot_symmetry_diagnostics, plot_mirror_overlap, default_tol)

# 直接以脚本方式运行（而非 notebook）时，Windows 控制台默认 GBK 会让中文/emoji 打印炸掉
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

FIGS.mkdir(parents=True, exist_ok=True)
(DERIVED / "UshapeNormal_New_prob").mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # 2. 参数配置
#
# ## 2.1 数据集与输出
# ## 2.2 对称面 / 镜像开关
# ## 2.3 磁铁几何
#
# > ⚠️ **几何参数是占位值**：下面这组数是从扫描范围反推的粗估，**必须替换成你实测的
# > 磁铁尺寸**。填充后把 `GEOMETRY_IS_PLACEHOLDER` 置为 `False`。
# > 反推依据：`magnetic_field_inside.csv` 覆盖 x∈[125,157]（⇒ 槽宽 ≈33 mm、X_c≈141）、
# > N/S 扫描分别在 x≤94 / x≥189（⇒ 外侧面在 96 / 186 附近）、front 扫描 y≤117
# > （⇒ 磁腿下端面 y≈119）。z 厚度 T 无法从外部场反演（对 z 均匀分布时外部场只依赖 M·T），
# > 必须用卡尺量。

# %%
# ==========================================
# 2.1 数据集与输出
# ==========================================
DATA_DIR = str(RAW / "UshapeNormal_New_prob")
csv_pattern = os.path.join(DATA_DIR, "*.csv")
output_vtk_base = str(DERIVED / "UshapeNormal_New_prob" / "inverted_M_results_newmesh")
B_unit_conversion = 1e-6            # 输入 B 的单位 -> T（这里输入是 uT）

# ==========================================
# 2.2 对称面 / 镜像
# ==========================================
SYMMETRY_MODE = "auto"              # "auto" 自动检测 | "manual" 用 MANUAL_H | "off" 不镜像
MANUAL_H = 8.0                      # 仅 SYMMETRY_MODE="manual" 时生效 [mm]
OVERLAP_POLICY = "keep_measured"    # "keep_measured" 实测优先 | "average" 重合对对称化平均
MIRROR_FILL = True                  # 是否把镜像点补进实测未覆盖的区域
RECENTER_Z_BY_H = True              # 用检测到的 h 重新锚定 Z_MIN/MAX = h ∓ T/2

# ==========================================
# 2.3 磁铁几何（单位：mm）—— 占位值，请替换
# ==========================================
GEOMETRY_IS_PLACEHOLDER = True
MAGNET_X_MIN, MAGNET_X_MAX = 96.0, 186.0    # 占位：X_c=141, R_out=45
MAGNET_Y_MIN, MAGNET_Y_MAX = 120.0, 250.0   # 占位：Y_c=Y_MAX-R_out=205
R_in = 18.0                                 # 占位：内半径。inside 扫描覆盖 x∈[125,157]，
                                            #       故槽半宽必须 >16（探头要留间隙），取 18
Z_THICKNESS = 32.0                          # 占位：卡尺实测厚度
MAGNET_Z_MIN, MAGNET_Z_MAX = -8.0, 24.0     # 占位：中面 8、厚度 32
GEOMETRY_MARGIN = 0.5                       # 几何校验余量 [mm]
GEOMETRY_CHECK = "error"                    # "error" 有测点落在磁铁内就报错 | "warn" 只告警

# 由上面派生
X_c = (MAGNET_X_MIN + MAGNET_X_MAX) / 2.0
R_out = (MAGNET_X_MAX - MAGNET_X_MIN) / 2.0
Y_c = MAGNET_Y_MAX - R_out

# ---- 自适应网格参数 ----
voxel_size_base = 1.5
Y_split = 144.0

# ---- 反演参数 ----
huber_epsilon = 4e3
lambda_reg = 1e-10
max_M = 1.5e6
max_iter = 50

if GEOMETRY_IS_PLACEHOLDER:
    print("=" * 78)
    print("⚠️  注意：磁铁几何仍是【占位值】，请替换成你实测的尺寸后再使用反演结果！")
    print(f"   X[{MAGNET_X_MIN},{MAGNET_X_MAX}] Y[{MAGNET_Y_MIN},{MAGNET_Y_MAX}] "
          f"R_in={R_in} T={Z_THICKNESS} Z[{MAGNET_Z_MIN},{MAGNET_Z_MAX}]")
    print("=" * 78)

# %% [markdown]
# # 3. 读取磁场数据
#
# 多文件读取（同原版）。这里额外统计"实测集内部的重复坐标"，但默认 `dedup="warn"`
# **不修改数据**，以便与旧结果保持一致。

# %%
P_meas, B_meas, meta = load_measurements(csv_pattern, dedup="warn")

# %% [markdown]
# # 4. 对称面自动检测 + 镜像增强
#
# `estimate_symmetry_plane` 的原理与判据：
#
# - 按 (x,y) 组成 z 列，用三次样条把每条列重采样到统一 z 网格；
# - 目标函数 $J(h)=\frac{1}{3}\sum_k \frac{\sum \Delta B_k^2}{\sum B_k^2}$，
#   其中 $\Delta \mathbf{B}=\mathbf{B}(2h-z)-(B_x,B_y,-B_z)(z)$ 只在"镜像点也落在测量 z 跨度内"的窗口统计；
# - $J\ll 1$ ⇒ 数据确实关于该平面镜像对称；$J\approx1$ ⇒ 完全没有这种对称性（单侧数据即如此）；
# - 接受判据：配对数足够、$J_{min}\le0.05$ 且 $J_{min}\le0.2\,J_{median}$、凹陷带宽 ≤2 mm；
# - 若最优值离"2h 落在测量 z 网格上"的候选在带宽内，则吸附过去，避免镜像插值误差；
# - 不接受时回落几何中面 `(Z_MIN+Z_MAX)/2`，行为与旧版完全一致。

# %%
fallback_h = 0.5 * (MAGNET_Z_MIN + MAGNET_Z_MAX)
sym = None

if SYMMETRY_MODE == "off":
    h_sym = fallback_h
    print(f"SYMMETRY_MODE='off'：不做镜像。")
elif SYMMETRY_MODE == "manual":
    h_sym = float(MANUAL_H)
    print(f"SYMMETRY_MODE='manual'：对称面 h = {h_sym:.3f} mm")
else:
    sym = estimate_symmetry_plane(P_meas, B_meas)
    if sym.accepted:
        h_sym = float(sym.h)
    else:
        h_sym = fallback_h
        print(f"⚠️  自动检测未通过，回落几何中面 h = {h_sym:.3f} mm")
        print(f"    原因: {sym.reason}")

print(f"\n采用对称面高度 h = {h_sym:.3f} mm")

# %%
# ---- 镜像 + 重合处理 ----
if SYMMETRY_MODE == "off":
    measured_points = P_meas / 1000.0
    B_measured = B_meas * B_unit_conversion
    mirror_rep = None
    print(f"不镜像：{len(measured_points)} 个测点")
else:
    P_all, B_all, mirror_rep = apply_mirror(
        P_meas, B_meas, h_sym,
        policy=OVERLAP_POLICY, fill=MIRROR_FILL, meta=meta,
        n_duplicate_measured=meta.get("n_dup", 0))
    measured_points = P_all / 1000.0           # m
    B_measured = B_all * B_unit_conversion     # T
    # 镜像补点的权重：本版本不加权，全部为 1（重合处已按策略去重，不存在重复计权）

B_data = B_measured.ravel()
# .ravel() 默认行优先：索引 = 3*i + k，k∈{0,1,2} 对应 x,y,z
print(f"数据装配完成：{len(B_measured)} 个测点，B_data 模长最大值 = {np.max(np.abs(B_data)):.6e} T")

# %%
# ---- 诊断图 ----
if sym is not None:
    plot_symmetry_diagnostics(sym, P_meas, B_meas,
                              str(FIGS / "symmetry_plane_detection.png"),
                              fallback_h=fallback_h, title_suffix=f"  [{DATA_DIR}]")
if mirror_rep is not None and mirror_rep.n_overlap:
    plot_mirror_overlap(mirror_rep, P_meas, B_meas,
                        str(FIGS / "mirror_overlap.png"))

# %% [markdown]
# # 5. 几何与 z 边界锚定 + 一致性校验
#
# `RECENTER_Z_BY_H=True` 时用检测到的 h 重新锚定 z 边界：`Z_MIN/MAX = h ∓ T/2`。
# 这样磁铁的厚度中面正好落在 z 基网格的格边界上（`nz_base` 是 2 的幂，必为偶数）。
#
# 校验规则：**任何测点（含镜像补点）都不允许落进磁铁实体内部**——探头不可能在磁铁里，
# 而且偶极子正演在源区内不成立。

# %%
if RECENTER_Z_BY_H:
    MAGNET_Z_MIN = h_sym - Z_THICKNESS / 2.0
    MAGNET_Z_MAX = h_sym + Z_THICKNESS / 2.0
    print(f"已按 h 重新锚定 z 边界：Z ∈ [{MAGNET_Z_MIN:.3f}, {MAGNET_Z_MAX:.3f}] mm"
          f"（厚度 {Z_THICKNESS}，中面 {(MAGNET_Z_MIN + MAGNET_Z_MAX) / 2:.3f}）")
else:
    mid = 0.5 * (MAGNET_Z_MIN + MAGNET_Z_MAX)
    print(f"沿用给定 z 边界：Z ∈ [{MAGNET_Z_MIN}, {MAGNET_Z_MAX}] mm，"
          f"中面 {mid:.3f}，与检测到的 h={h_sym:.3f} 相差 {abs(mid - h_sym):.3f} mm")


def u_solid_mask(points_mm, margin=0.0):
    """点是否落在 U 形磁铁实体内部（含 z 范围），margin 为向外扩张量 [mm]"""
    x = points_mm[:, 0]
    y = points_mm[:, 1]
    z = points_mm[:, 2]
    in_z = (z >= MAGNET_Z_MIN - margin) & (z <= MAGNET_Z_MAX + margin)
    dist = np.hypot(x - X_c, y - Y_c)
    in_arc = (y > Y_c) & (dist >= R_in - margin) & (dist <= R_out + margin)
    in_legs = ((y >= MAGNET_Y_MIN - margin) & (y <= Y_c) &
               ((x >= X_c - R_out - margin) & (x <= X_c - R_in + margin) |
                (x >= X_c + R_in - margin) & (x <= X_c + R_out + margin)))
    return in_z & (in_arc | in_legs)


bad = u_solid_mask(measured_points * 1000.0, GEOMETRY_MARGIN)
n_bad = int(bad.sum())
if n_bad:
    src = np.array(["实测"] * len(P_meas) + ["镜像补点"] * (len(measured_points) - len(P_meas))) \
        if len(measured_points) != len(P_meas) else np.array(["实测"] * len(P_meas))
    msg = [f"⚠️  几何校验失败：{n_bad} / {len(measured_points)} 个测点落在磁铁实体内部"
           f"（余量 {GEOMETRY_MARGIN} mm），说明几何参数与本次坐标系不一致。",
           f"  其中实测点 {int((bad & (src == '实测')).sum())} 个，镜像补点 "
           f"{int((bad & (src == '镜像补点')).sum())} 个。"]
    for name, mask in [("实测", src == "实测"), ("镜像补点", src == "镜像补点")]:
        if (bad & mask).any():
            q = measured_points[bad & mask] * 1000.0
            msg.append(f"  {name}: x[{q[:, 0].min():.1f},{q[:, 0].max():.1f}] "
                       f"y[{q[:, 1].min():.1f},{q[:, 1].max():.1f}] "
                       f"z[{q[:, 2].min():.1f},{q[:, 2].max():.1f}] mm")
    msg.append("  请检查 MAGNET_X_MIN/MAX、MAGNET_Y_MIN/MAX、R_in、Z_THICKNESS。")
    print("\n".join(msg))
    if GEOMETRY_CHECK == "error":
        raise RuntimeError("几何参数与测量数据不一致，已停止（把 GEOMETRY_CHECK 设为 'warn' "
                           "可只告警继续）。")
else:
    print(f"✅ 几何校验通过：{len(measured_points)} 个测点全部在磁铁实体之外"
          f"（余量 {GEOMETRY_MARGIN} mm）。")

# %% [markdown]
# # 6. 可变网格
#
# 这里需要注意，返回 -1 是最细的网格（basemesh），返回 0 是最粗的网格。
# 与原版一致：先按 bbox 细化到次细层级，再按 U 形几何与 y 坐标决定细化级别。

# %%
print("正在构建 TreeMesh 自适应网格...")

nx_desired = int(np.ceil((MAGNET_X_MAX - MAGNET_X_MIN) / voxel_size_base))
ny_desired = int(np.ceil((MAGNET_Y_MAX - MAGNET_Y_MIN) / voxel_size_base))
nz_desired = int(np.ceil((MAGNET_Z_MAX - MAGNET_Z_MIN) / voxel_size_base))


def next_pow2(n):
    return 2 ** int(np.ceil(np.log2(n)))


nx_base = next_pow2(nx_desired)
ny_base = next_pow2(ny_desired)
nz_base = next_pow2(nz_desired)

dx = (MAGNET_X_MAX - MAGNET_X_MIN) / nx_base
dy = (MAGNET_Y_MAX - MAGNET_Y_MIN) / ny_base
dz = (MAGNET_Z_MAX - MAGNET_Z_MIN) / nz_base

mesh = TreeMesh(
    [[(dx, nx_base)], [(dy, ny_base)], [(dz, nz_base)]],
    origin=[MAGNET_X_MIN, MAGNET_Y_MIN, MAGNET_Z_MIN],
    diagonal_balance=False
)

max_level = mesh.max_level
absolute_fine_level = max_level
absolute_coarse_level = max_level - 1

BBox = np.array([[MAGNET_X_MIN, MAGNET_Y_MIN, MAGNET_Z_MIN],
                 [MAGNET_X_MAX, MAGNET_Y_MAX, MAGNET_Z_MAX]])
mesh.refine_bounding_box(BBox, level=absolute_coarse_level, finalize=False)

print(f"基础网格实际单元尺寸: dx={dx:.2f}, dy={dy:.2f}, dz={dz:.2f} mm")
print(f"基础网格单元数: {nx_base}*{ny_base}*{nz_base}")
print(f"z 厚度中面放在格边界上: {(0.5 * (MAGNET_Z_MIN + MAGNET_Z_MAX) - MAGNET_Z_MIN) / dz:.3f} 个格")


def u_inside_xy(x, y):
    """U 形磁铁在 x-y 平面内的实体判定（与原版公式一致）"""
    if y > Y_c:
        dist = np.sqrt((x - X_c) ** 2 + (y - Y_c) ** 2)
        return (R_in <= dist) and (dist <= R_out)
    return (MAGNET_Y_MIN <= y <= Y_c) and \
           ((X_c - R_out <= x <= X_c - R_in) or (X_c + R_in <= x <= X_c + R_out))


def refine_func(cell):
    x, y, z = cell.center
    if not u_inside_xy(x, y):
        return 0
    return absolute_fine_level if y < Y_split else absolute_coarse_level


mesh.refine(refine_func, finalize=False)
mesh.finalize()
mesh.number()

all_cc = mesh.cell_centers          # mm
all_vol = mesh.cell_volumes         # mm^3

is_magnet = np.array([u_inside_xy(x, y) for x, y in zip(all_cc[:, 0], all_cc[:, 1])],
                     dtype=bool)

magnet_idx = np.asarray(is_magnet).nonzero()[0]
voxel_points_mm = all_cc[magnet_idx]
voxel_points = voxel_points_mm / 1000.0
voxel_vol = all_vol[magnet_idx] / 1000 ** 3
n_voxels = len(voxel_points)

print(f"自适应网格：总单元 {mesh.nC}，磁铁单元 {n_voxels}")
print(f"最小单元尺寸约 {min(all_vol) ** (1 / 3):.1f} mm，最大约 {max(all_vol) ** (1 / 3):.1f} mm")

# %% [markdown]
# # 7. 构建敏感度矩阵 $A$
#
# $B = A \cdot M$，$A$ 是 $3M \times 3N$ 矩阵。用 `LinearOperator` 隐式表示（matvec/rmatvec 用
# numba 并行实现），不显式构建稠密矩阵以节省内存。**与原版完全相同。**

# %%
print("正在构建敏感度矩阵 LinearOperator...")


@njit(parallel=True)
def matvec(M_vec, meas, voxels, vols):
    M = meas.shape[0]
    N = voxels.shape[0]
    out = np.zeros(M * 3)
    mu0_4pi = 1e-7
    for i in prange(M):
        xi = meas[i, 0]; yi = meas[i, 1]; zi = meas[i, 2]
        bx = 0.0; by = 0.0; bz = 0.0
        for j in range(N):
            dx = xi - voxels[j, 0]; dy = yi - voxels[j, 1]; dz = zi - voxels[j, 2]
            r2 = dx * dx + dy * dy + dz * dz
            if r2 < 1e-8:
                r2 = 1e-8
            r = np.sqrt(r2); r3 = r2 * r; r5 = r3 * r2
            coef = mu0_4pi * vols[j]
            mx = M_vec[3 * j]; my = M_vec[3 * j + 1]; mz = M_vec[3 * j + 2]
            md = mx * dx + my * dy + mz * dz
            k = coef * (3.0 * md / r5)
            bx += k * dx - coef * mx / r3
            by += k * dy - coef * my / r3
            bz += k * dz - coef * mz / r3
        out[3 * i] = bx; out[3 * i + 1] = by; out[3 * i + 2] = bz
    return out


@njit(parallel=True)
def rmatvec(B_vec, meas, voxels, vols):
    M = meas.shape[0]
    N = voxels.shape[0]
    out = np.zeros(N * 3)
    mu0_4pi = 1e-7
    for j in prange(N):
        xj = voxels[j, 0]; yj = voxels[j, 1]; zj = voxels[j, 2]
        coef = mu0_4pi * vols[j]
        mx = 0.0; my = 0.0; mz = 0.0
        for i in range(M):
            dx = meas[i, 0] - xj; dy = meas[i, 1] - yj; dz = meas[i, 2] - zj
            r2 = dx * dx + dy * dy + dz * dz
            if r2 < 1e-8:
                r2 = 1e-8
            r = np.sqrt(r2); r3 = r2 * r; r5 = r3 * r2
            bx = B_vec[3 * i]; by = B_vec[3 * i + 1]; bz = B_vec[3 * i + 2]
            bd = bx * dx + by * dy + bz * dz
            k = coef * (3.0 * bd / r5)
            mx += k * dx - coef * bx / r3
            my += k * dy - coef * by / r3
            mz += k * dz - coef * bz / r3
        out[3 * j] = mx; out[3 * j + 1] = my; out[3 * j + 2] = mz
    return out


M_all = len(measured_points)
N_all = n_voxels
A_op = LinearOperator((M_all * 3, N_all * 3),
                      matvec=lambda v: matvec(np.ascontiguousarray(v, dtype=np.float64),
                                              measured_points, voxel_points, voxel_vol),
                      rmatvec=lambda w: rmatvec(np.ascontiguousarray(w, dtype=np.float64),
                                                measured_points, voxel_points, voxel_vol),
                      dtype=np.float64)

# %% [markdown]
# # 8. 构建正则化矩阵
#
# 磁化强度不允许突变，因此把梯度作为惩罚项；为了突出接缝，梯度大于阈值时降低权重
# （类 Huber / Charbonnier）。$G$ 来自 `mesh.cell_gradient`，$W$ 来自 `face_inner_product`。
# **与原版完全相同。**

# %%
print("正在构建正则化矩阵（cell_gradient + face_inner_product）...")

G_full = mesh.cell_gradient
W = mesh.get_face_inner_product()
if not isinstance(G_full, csr_matrix):
    G_full = G_full.tocsr()
V_all = W.diagonal()

row_start = G_full.indptr
col_indices = G_full.indices
internal_faces = []
for face_idx in range(G_full.shape[0]):
    cols = col_indices[row_start[face_idx]:row_start[face_idx + 1]]
    if len(cols) >= 2:
        if is_magnet[cols].all():
            internal_faces.append(face_idx)

internal_faces = np.array(internal_faces)
n_valid_faces = len(internal_faces)
V_face = V_all[internal_faces]
print(f"内部约束面数量（含悬挂面）: {n_valid_faces}")
print(f"面对偶体积 V_face: min={V_face.min():.3f}, max={V_face.max():.3f} mm^3"
      f"（{len(np.unique(np.round(V_face, 3)))} 个取值）")

G_sub = G_full[internal_faces, :][:, magnet_idx]

# %% [markdown]
# # 9. 求解（L-BFGS-B + Huber 边缘保持正则化）
#
# $$\phi_m = \lambda \sum_f \epsilon^2\left(\sqrt{1+\left(\frac{\|\nabla M\|_f}{\epsilon}\right)^2}-1\right)$$

# %%
print("开始 L-BFGS-B 高效反演...")

M_init = np.zeros(3 * n_voxels)
bounds = [(-max_M, max_M)] * (3 * n_voxels)


def objective_and_gradient(M_vec):
    M_sol_current = M_vec.reshape(n_voxels, 3)

    residual = A_op.matvec(M_vec) - B_data
    loss_data = np.sum(residual ** 2)
    grad_data = 2.0 * (A_op.rmatvec(residual))

    grad = G_sub @ M_sol_current
    grad_norm = np.sqrt(np.sum(grad ** 2, axis=1))
    loss_reg = lambda_reg * huber_epsilon ** 2 * np.sum(
        V_face * (np.sqrt(1.0 + (grad_norm / huber_epsilon) ** 2) - 1.0))

    weights = 1.0 / np.sqrt(1.0 + (grad_norm / huber_epsilon) ** 2)
    grad_reg = lambda_reg * (G_sub.T @ (V_face[:, None] * grad * weights[:, None])).ravel()

    return loss_data + loss_reg, grad_data + grad_reg


res_opt = minimize(
    fun=objective_and_gradient,
    x0=M_init,
    jac=True,
    method="L-BFGS-B",
    bounds=bounds,
    options={"maxiter": max_iter, "ftol": 1e-30, "gtol": 1e-30}
)

M_sol = res_opt.x.reshape(n_voxels, 3)
print(f"L-BFGS-B 反演完成！最大反演磁化强度: {np.max(np.abs(M_sol)):.3e} A/m")
print(f"最终目标函数值: {res_opt.fun:.6e}")

# %% [markdown]
# # 10. 导出 VTK
#
# 空气区域磁化强度保持 0，与原版一致。

# %%
M_full = np.zeros((mesh.nC, 3))
M_full[magnet_idx, :] = M_sol
model_dict = {
    "M_vector": M_full,
    "Mx": M_full[:, 0],
    "My": M_full[:, 1],
    "Mz": M_full[:, 2],
}
print(f"正在直接写入 VTK (UnstructuredGrid) 文件: {output_vtk_base}.vtu ...")
mesh.write_vtk(output_vtk_base, models=model_dict)
print(f"✅ VTK 导出成功！文件已写入 {output_vtk_base}.vtu")

# %% [markdown]
# # 11. 结果速览
#
# 取磁铁中面附近一层的 $|\mathbf{M}|$ 分布画出来，便于快速检查。

# %%
mag = np.linalg.norm(M_sol, axis=1)
z_mag = voxel_points_mm[:, 2]
z_layers = np.unique(np.round(z_mag, 1))
z_target = 0.5 * (MAGNET_Z_MIN + MAGNET_Z_MAX)
if len(z_layers):
    z_layer = z_layers[np.argmin(np.abs(z_layers - z_target))]
    sel = np.isclose(z_mag, z_layer, atol=0.15) & (mag > 0)
    fig, ax = plt.subplots(figsize=(6.4, 5.4))
    sc = ax.scatter(voxel_points_mm[sel, 0], voxel_points_mm[sel, 1], c=mag[sel] / 1000.0,
                    s=4, cmap="inferno")
    fig.colorbar(sc, ax=ax, label="磁化强度 |M| / (kA/m)")
    th = np.linspace(0, np.pi, 200)
    ax.plot(X_c + R_out * np.cos(th), Y_c + R_out * np.sin(th), "k-", lw=1.0)
    ax.plot(X_c + R_in * np.cos(th), Y_c + R_in * np.sin(th), "k-", lw=1.0)
    for x0, x1 in ((X_c - R_out, X_c - R_in), (X_c + R_in, X_c + R_out)):
        ax.plot([x0, x0], [Y_c, MAGNET_Y_MIN], "k-", lw=1.0)
        ax.plot([x1, x1], [Y_c, MAGNET_Y_MIN], "k-", lw=1.0)
        ax.plot([x0, x1], [MAGNET_Y_MIN, MAGNET_Y_MIN], "k-", lw=1.0)
    ax.set_aspect("equal")
    ax.set_xlabel("x / mm")
    ax.set_ylabel("y / mm")
    ax.set_title(f"|M| 分布（z = {z_layer:g} mm）")
    fig.tight_layout()
    fig.savefig(FIGS / "M_newprob_section.png", dpi=200)
    plt.close(fig)
    print(f"已保存 {FIGS}/M_newprob_section.png（z = {z_layer:g} mm，{sel.sum()}/{len(z_mag[z_mag == z_layer])} 个单元 |M|>0）")

# %% [markdown]
# # 12.（可选）三种数据装配模式的对比
#
# 把 `RUN_MODE_COMPARISON` 打开会各跑一次同参数反演，用来回答"镜像到底值不值得"：
#
# - `off`：完全不镜像，只用实测数据；
# - `keep_measured`：镜像补空区，但与实测点重合的镜像副本一律丢弃（默认策略）；
# - `average`：重合对做对称化平均（等价于两条数据各 0.5 权），会强行施加对称性。
#
# 对比相对残差、`max|M|` 和中面 |M| 分布。**注意很慢**：每个模式都要跑 `COMPARE_MAXITER` 次迭代。

# %%
RUN_MODE_COMPARISON = False
COMPARE_MAXITER = 20

if RUN_MODE_COMPARISON:
    def build_operator(points_m):
        return LinearOperator((len(points_m) * 3, n_voxels * 3),
                              matvec=lambda v: matvec(np.ascontiguousarray(v, dtype=np.float64),
                                                      points_m, voxel_points, voxel_vol),
                              rmatvec=lambda w: rmatvec(np.ascontiguousarray(w, dtype=np.float64),
                                                        points_m, voxel_points, voxel_vol),
                              dtype=np.float64)

    def solve_with(points_m, B_t, maxiter=COMPARE_MAXITER):
        A_cmp = build_operator(points_m)
        b_cmp = B_t.ravel()

        def fg(M_vec):
            r = A_cmp.matvec(M_vec) - b_cmp
            g = G_sub @ M_vec.reshape(n_voxels, 3)
            gn = np.sqrt(np.sum(g ** 2, axis=1))
            wgt = 1.0 / np.sqrt(1.0 + (gn / huber_epsilon) ** 2)
            loss = np.sum(r ** 2) + lambda_reg * huber_epsilon ** 2 * np.sum(
                V_face * (np.sqrt(1.0 + (gn / huber_epsilon) ** 2) - 1.0))
            grad = 2.0 * A_cmp.rmatvec(r) + \
                lambda_reg * (G_sub.T @ (V_face[:, None] * g * wgt[:, None])).ravel()
            return loss, grad

        res_cmp = minimize(fg, np.zeros(3 * n_voxels), jac=True, method="L-BFGS-B",
                           bounds=[(-max_M, max_M)] * (3 * n_voxels),
                           options={"maxiter": maxiter, "ftol": 1e-30, "gtol": 1e-30})
        M_cmp = res_cmp.x.reshape(n_voxels, 3)
        misfit = np.linalg.norm(A_cmp.matvec(res_cmp.x) - b_cmp) / np.linalg.norm(b_cmp)
        return M_cmp, misfit

    configs = [("off", None), ("keep_measured", "keep_measured"), ("average", "average")]
    # 本节自带截面选择，不依赖第 11 节（顺序执行时结果一致）
    z_mag_cmp = voxel_points_mm[:, 2]
    z_layers_cmp = np.unique(np.round(z_mag_cmp, 1))
    z_mid_cmp = 0.5 * (MAGNET_Z_MIN + MAGNET_Z_MAX)
    z_layer_cmp = z_layers_cmp[np.argmin(np.abs(z_layers_cmp - z_mid_cmp))]
    fig, axes = plt.subplots(1, len(configs), figsize=(4.6 * len(configs), 4.6))
    vmax = 0.0
    out = {}
    for ax, (name, pol) in zip(np.atleast_1d(axes), configs):
        if pol is None:
            Pm, Bt = P_meas / 1000.0, B_meas * B_unit_conversion
        else:
            P2, B2, _ = apply_mirror(P_meas, B_meas, h_sym, policy=pol,
                                     fill=MIRROR_FILL, verbose=False)
            Pm, Bt = P2 / 1000.0, B2 * B_unit_conversion
        M_cmp, misfit = solve_with(Pm, Bt)
        m_cmp = np.linalg.norm(M_cmp, axis=1)
        out[name] = (m_cmp, misfit, len(Pm))
        vmax = max(vmax, np.percentile(m_cmp, 99.5))
        print(f"{name:14s} 测点 {len(Pm):6d}  相对残差 {misfit:.4e}  "
              f"max|M| {np.abs(M_cmp).max():.3e} A/m")
    for ax, (name, pol) in zip(np.atleast_1d(axes), configs):
        m_cmp, misfit, npts = out[name]
        sel = np.isclose(z_mag_cmp, z_layer_cmp, atol=0.15) & (m_cmp > 0)
        sc = ax.scatter(voxel_points_mm[sel, 0], voxel_points_mm[sel, 1],
                        c=m_cmp[sel] / 1000.0, s=4, cmap="inferno", vmin=0, vmax=vmax / 1000.0)
        th = np.linspace(0, np.pi, 200)
        ax.plot(X_c + R_out * np.cos(th), Y_c + R_out * np.sin(th), "k-", lw=1.0)
        ax.plot(X_c + R_in * np.cos(th), Y_c + R_in * np.sin(th), "k-", lw=1.0)
        for x0, x1 in ((X_c - R_out, X_c - R_in), (X_c + R_in, X_c + R_out)):
            ax.plot([x0, x0], [Y_c, MAGNET_Y_MIN], "k-", lw=1.0)
            ax.plot([x1, x1], [Y_c, MAGNET_Y_MIN], "k-", lw=1.0)
            ax.plot([x0, x1], [MAGNET_Y_MIN, MAGNET_Y_MIN], "k-", lw=1.0)
        ax.set_aspect("equal")
        ax.set_xlabel("x / mm")
        ax.set_title(f"{name}\n测点 {npts}，相对残差 {misfit:.2e}")
    fig.colorbar(sc, ax=list(np.atleast_1d(axes)), label="|M| / (kA/m)", pad=0.02)
    fig.suptitle(f"三种数据装配模式对比（h = {h_sym:.2f} mm, z = {z_layer_cmp:g} mm）")
    fig.savefig(FIGS / "M_mode_comparison.png", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"已保存 {FIGS}/M_mode_comparison.png")
