# %% [markdown]
# # 无正则化最小二乘反演（B2M_noreg）
#
# 本 notebook 是 `B2M.ipynb` 的**无正则化矩阵版本**：保留数据读取、Z 镜像增强、
# TreeMesh 自适应网格与偶极子正演，但**去掉 Huber 正则化项与 `[-max_M, max_M]` 边界**，
# 直接求解线性最小二乘问题
#
# $$\min_{\mathbf{M}} \|\mathbf{A}\mathbf{M} - \mathbf{B}\|^2$$
#
# 目的是说明：没有人为施加的限制（平滑正则化），反演结果会发散、振荡，无法给出
# 物理上合理的磁化强度，从而证明 `B2M.ipynb` 中正则化帮助收敛的必要性。

# %%
# ============================================
# 1. 导入与参数配置（与 B2M.ipynb 一致）
# ============================================
import glob
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import lsqr, LinearOperator
import discretize
from discretize import TreeMesh
from numba import njit, prange
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pyvista as pv
import os

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False
os.makedirs("figures", exist_ok=True)

csv_pattern = "UshapeNormal/*.csv"
output_vtk_base = "UshapeNormal/inverted_M_results_noreg"
B_unit_conversion = 1e-6

MAGNET_X_MIN, MAGNET_X_MAX = 111.5, 193.5
MAGNET_Y_MIN, MAGNET_Y_MAX = 129.5, 231.5
MAGNET_Z_MIN, MAGNET_Z_MAX = -35, -3
X_c = (MAGNET_X_MIN + MAGNET_X_MAX) / 2.0
R_out = (MAGNET_X_MAX - MAGNET_X_MIN) / 2.0
Y_c = MAGNET_Y_MAX - R_out
R_in = 22
voxel_size_base = 1.5
Y_split = 144.0
max_M = 1.5e6  # 仅用于正则化版本对照

# %% [markdown]
# # 2. 读取磁场数据 + Z 镜像增强
#
# 与 `B2M.ipynb` 完全相同：多文件读取、Z 向镜像，`B_data` 展平为
# $[B_{x0}, B_{y0}, B_{z0}, B_{x1}, \dots]$。

# %%
csv_files = glob.glob(csv_pattern)
df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)
measured_points_origin = df[['x', 'y', 'z']].values / 1000.0
B_measured_origin = df[['Bx', 'By', 'Bz']].values * B_unit_conversion
Z_center_m = ((MAGNET_Z_MIN + MAGNET_Z_MAX) / 2.0) / 1000.0
measured_points_mirrored = measured_points_origin.copy()
measured_points_mirrored[:, 2] = 2 * Z_center_m - measured_points_origin[:, 2]
B_measured_mirrored = B_measured_origin.copy()
B_measured_mirrored[:, 2] *= -1.0
measured_points = np.vstack([measured_points_origin, measured_points_mirrored])
B_measured = np.vstack([B_measured_origin, B_measured_mirrored])
B_data = B_measured.ravel()
print(f"数据合并完成，共 {len(B_measured)} 个测点（含镜像）。")

# %% [markdown]
# # 3. TreeMesh 自适应网格 + 磁铁掩膜
#
# 与 `B2M.ipynb` 完全相同。

# %%
def next_pow2(n):
    return 2 ** int(np.ceil(np.log2(n)))

nx_base = next_pow2(int(np.ceil((MAGNET_X_MAX - MAGNET_X_MIN) / voxel_size_base)))
ny_base = next_pow2(int(np.ceil((MAGNET_Y_MAX - MAGNET_Y_MIN) / voxel_size_base)))
nz_base = next_pow2(int(np.ceil((MAGNET_Z_MAX - MAGNET_Z_MIN) / voxel_size_base)))
dx = (MAGNET_X_MAX - MAGNET_X_MIN) / nx_base
dy = (MAGNET_Y_MAX - MAGNET_Y_MIN) / ny_base
dz = (MAGNET_Z_MAX - MAGNET_Z_MIN) / nz_base
mesh = TreeMesh([[(dx, nx_base)], [(dy, ny_base)], [(dz, nz_base)]],
                origin=[MAGNET_X_MIN, MAGNET_Y_MIN, MAGNET_Z_MIN], diagonal_balance=False)
max_level = mesh.max_level
mesh.refine_bounding_box(np.array([[MAGNET_X_MIN, MAGNET_Y_MIN, MAGNET_Z_MIN],
                                   [MAGNET_X_MAX, MAGNET_Y_MAX, MAGNET_Z_MAX]]),
                         level=max_level - 1)

def refine_func(cell):
    x, y, z = cell.center
    if y > Y_c:
        dist = np.sqrt((x - X_c) ** 2 + (y - Y_c) ** 2)
        inside = (R_in <= dist <= R_out)
    else:
        inside = (MAGNET_Y_MIN <= y <= Y_c) and \
                 ((X_c - R_out <= x <= X_c - R_in) or (X_c + R_in <= x <= X_c + R_out))
    if not inside:
        return 0
    return max_level if y < Y_split else max_level - 1

mesh.refine(refine_func)
mesh.number()
all_cc = mesh.cell_centers
all_vol = mesh.cell_volumes
is_magnet = np.zeros(mesh.nC, dtype=bool)
for i, (x, y, z) in enumerate(all_cc):
    if y > Y_c:
        dist = np.sqrt((x - X_c) ** 2 + (y - Y_c) ** 2)
        is_magnet[i] = (R_in <= dist <= R_out)
    else:
        is_magnet[i] = (MAGNET_Y_MIN <= y <= Y_c) and \
                       ((X_c - R_out <= x <= X_c - R_in) or (X_c + R_in <= x <= X_c + R_out))
magnet_idx = np.asarray(is_magnet).nonzero()[0]
voxel_points_mm = all_cc[magnet_idx]
voxel_points = voxel_points_mm / 1000.0
voxel_vol = all_vol[magnet_idx] / 1000 ** 3
n_voxels = len(voxel_points)
print(f"自适应网格：总单元 {mesh.nC}，磁铁单元 {n_voxels}")

# %% [markdown]
# # 4. 敏感度矩阵 A（偶极子正演，LinearOperator）
#
# 不显式构建 $3M \times 3N$ 的稠密矩阵（约 26 GB），改用 `LinearOperator`
# 隐式表示，`matvec`/`rmatvec` 用 numba 加速的偶极子核。

# %%
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
# # 5. 无正则化最小二乘求解
#
# 用 `scipy.sparse.linalg.lsqr` 直接求解 $\min\|\mathbf{A}\mathbf{M}-\mathbf{B}\|^2$，
# **不加正则化、不加边界**。

# %%
print("开始无正则化最小二乘求解（lsqr）...")
sol = lsqr(A_op, B_data, atol=1e-9, btol=1e-9, conlim=1e12, iter_lim=500)
M_vec_noreg = sol[0]
istop, itn, r1norm, r2norm, anorm, acond, arnorm, xnorm = sol[1:9]
M_noreg = M_vec_noreg.reshape(n_voxels, 3)
print(f"lsqr 结束：istop={istop}，迭代 {itn} 次")
print(f"残差 ||A M - B|| = {r1norm:.6e}（初始 {r2norm:.6e}）")
print(f"条件数估计 acond = {acond:.3e}")
print(f"无正则化 max|M| = {np.max(np.abs(M_noreg)):.3e} A/m")
print(f"（正则化版本 B2M.ipynb 的 max|M| ≈ 2.77e5 A/m）")

# %% [markdown]
# # 6. 导出无正则化结果至 VTK
#
# 与 `B2M.ipynb` 相同的映射与导出方式（空气单元磁化强度置 0）。

# %%
M_full = np.zeros((mesh.nC, 3))
M_full[magnet_idx, :] = M_noreg
model_dict = {
    "M_vector": M_full,
    "Mx": M_full[:, 0],
    "My": M_full[:, 1],
    "Mz": M_full[:, 2],
}
mesh.write_vtk(output_vtk_base, models=model_dict)
print(f"已导出无正则化结果：{output_vtk_base}.vtu")

# %% [markdown]
# # 7. 对比：正则化 vs 无正则化（同一 z 截面）
#
# 在磁铁中面附近同一层（z = -20 mm）比较两种反演的 $|\mathbf{M}|$ 分布，
# 直观展示无正则化结果的发散与振荡。

# %%
def slice_at_z(cc, Mx, My, Mz, z_target):
    """取 |M|>0 的磁铁单元中离 z_target 最近的一层，返回 (x, y, |M|, Mx, My)。"""
    mag = np.sqrt(Mx**2 + My**2 + Mz**2)
    sel_mag = mag > 1.0
    z_mag = cc[sel_mag, 2]
    z_layers = np.unique(np.round(z_mag, 1))
    z_layer = z_layers[np.argmin(np.abs(z_layers - z_target))]
    sel = np.isclose(cc[:, 2], z_layer, atol=0.15) & sel_mag
    return cc[sel, 0], cc[sel, 1], mag[sel], Mx[sel], My[sel], z_layer

# 正则化结果（读 B2M.ipynb 导出的 VTU）
mesh_reg = pv.read("UshapeNormal/inverted_M_results_test.vtu")
cc_reg = mesh_reg.cell_centers().points
xr, yr, mr, mxr, myr, zl = slice_at_z(cc_reg,
                                      mesh_reg.cell_data["Mx"],
                                      mesh_reg.cell_data["My"],
                                      mesh_reg.cell_data["Mz"], -19.0)

# 无正则化结果（本 notebook 求得，按 magnet_idx 对应到 cell）
cc_nr = voxel_points_mm
xn, yn, mn, mxn, myn, _ = slice_at_z(cc_nr, M_noreg[:, 0], M_noreg[:, 1], M_noreg[:, 2], zl)

vmax = max(mr.max(), mn.max())
gi = np.linspace(min(xr.min(), xn.min()), max(xr.max(), xn.max()), 300)
gj = np.linspace(min(yr.min(), yn.min()), max(yr.max(), yn.max()), 360)
GX, GY = np.meshgrid(gi, gj)
grid_reg = griddata((xr, yr), mr, (GX, GY), method="linear")
grid_nr = griddata((xn, yn), mn, (GX, GY), method="linear")

def draw_u(ax):
    th = np.linspace(0, np.pi, 200)
    ax.plot(X_c + R_out * np.cos(th), Y_c + R_out * np.sin(th), color="k", lw=1.0)
    ax.plot(X_c + R_in * np.cos(th), Y_c + R_in * np.sin(th), color="k", lw=1.0)
    for x0, x1 in ((X_c - R_out, X_c - R_in), (X_c + R_in, X_c + R_out)):
        ax.plot([x0, x0], [Y_c, MAGNET_Y_MIN], color="k", lw=1.0)
        ax.plot([x1, x1], [Y_c, MAGNET_Y_MIN], color="k", lw=1.0)
        ax.plot([x0, x1], [MAGNET_Y_MIN, MAGNET_Y_MIN], color="k", lw=1.0)

fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.0))
im_r = axes[0].pcolormesh(GX, GY, grid_reg / 1000.0, cmap="inferno", vmin=0, vmax=vmax / 1000.0)
axes[0].set_title("正则化反演（B2M.ipynb）")
axes[0].set_aspect("equal")
draw_u(axes[0])
im_n = axes[1].pcolormesh(GX, GY, grid_nr / 1000.0, cmap="inferno", vmin=0, vmax=vmax / 1000.0)
axes[1].set_title("无正则化反演（本 notebook）")
axes[1].set_aspect("equal")
draw_u(axes[1])
for ax in axes:
    ax.set_xlabel("x / mm")
    ax.set_ylabel("y / mm")
fig.colorbar(im_n, ax=axes, label="磁化强度 |M| / (kA/m)", pad=0.02)
fig.suptitle(f"z = {zl:g} mm 截面：正则化 vs 无正则化")
fig.subplots_adjust(left=0.07, right=0.96, top=0.90, bottom=0.12, wspace=0.22)
fig.savefig("figures/M_regularized_vs_noreg.png", dpi=300)
fig.savefig("figures/M_regularized_vs_noreg.pdf")
plt.close(fig)
print(f"已保存 figures/M_regularized_vs_noreg.png / .pdf")
print(f"正则化 max|M| = {mr.max():.3e} A/m；无正则化 max|M| = {mn.max():.3e} A/m")
