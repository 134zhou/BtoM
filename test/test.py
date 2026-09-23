import numpy as np
import pandas as pd
from scipy.optimize import lsq_linear
import discretize
from discretize import TensorMesh

# ==========================================
# 1. 参数配置
# ==========================================
csv_filename = "./whole.csv"
output_m_csv = "inverted_M_results.csv"
B_unit_conversion = 1e-6

# 磁铁边界 (mm)
MAGNET_X_MIN, MAGNET_X_MAX = 106, 188
MAGNET_Y_MIN, MAGNET_Y_MAX = 127, 231
MAGNET_Z_MIN, MAGNET_Z_MAX = -28, 0

# 圆弧 U 形几何中心
X_c = (MAGNET_X_MIN + MAGNET_X_MAX) / 2.0
R_out = (MAGNET_X_MAX - MAGNET_X_MIN) / 2.0
Y_c = MAGNET_Y_MAX - R_out
R_in = 22

# 网格分辨率 (mm)
voxel_size = 5.0
voxel_volume = (voxel_size / 1000.0)**3   # m³

# ==========================================
# 2. 读取测量数据并 Z 轴镜像增强
# ==========================================
print("正在读取测量数据...")
try:
    df = pd.read_csv(csv_filename)
except FileNotFoundError:
    print(f"未找到 {csv_filename}，生成虚拟数据...")
    exit(1)

meas_pts_orig = df[['x', 'y', 'z']].values / 1000.0   # 转为米
b_meas_orig = df[['Bx', 'By', 'Bz']].values * B_unit_conversion

# Z 向镜像
Z_center_m = ((MAGNET_Z_MIN + MAGNET_Z_MAX) / 2.0) / 1000.0
meas_pts_mirrored = meas_pts_orig.copy()
meas_pts_mirrored[:, 2] = 2 * Z_center_m - meas_pts_orig[:, 2]

b_meas_mirrored = b_meas_orig.copy()
b_meas_mirrored[:, 0] = b_meas_orig[:, 0]
b_meas_mirrored[:, 1] = b_meas_orig[:, 1]
b_meas_mirrored[:, 2] = -b_meas_orig[:, 2]

meas_pts = np.vstack([meas_pts_orig, meas_pts_mirrored])
b_meas = np.vstack([b_meas_orig, b_meas_mirrored])
b_meas_vector = b_meas.ravel()

# ==========================================
# 3. 用 discretize 创建均匀立方体网格
# ==========================================
print("正在用 discretize.TensorMesh 构建网格...")

# 计算每个维度的单元数量（确保整除，必要时调整 voxel_size 或边界）
nx = int(round((MAGNET_X_MAX - MAGNET_X_MIN) / voxel_size))
ny = int(round((MAGNET_Y_MAX - MAGNET_Y_MIN) / voxel_size))
nz = int(round((MAGNET_Z_MAX - MAGNET_Z_MIN) / voxel_size))

# 每个方向的单元尺寸列表（全为 voxel_size）
hx = [voxel_size] * nx
hy = [voxel_size] * ny
hz = [voxel_size] * nz

# 创建网格，原点设为磁铁包围盒的最小角
mesh = TensorMesh([hx, hy, hz], origin=[MAGNET_X_MIN, MAGNET_Y_MIN, MAGNET_Z_MIN])

# 所有体素中心坐标 (mm)，形状 (nC, 3)
all_voxels_mm = mesh.cell_centers   # 注意：discretize 返回的坐标单位与 h 一致，此处为 mm

# ==========================================
# 4. 几何掩膜：区分磁铁实体与空气
# ==========================================
print("正在标识磁铁实体网格...")

is_magnet_mask = np.zeros(mesh.nC, dtype=bool)
# 按体素逐个判断（速度足够快）
for i, (vx, vy, vz) in enumerate(all_voxels_mm):
    if vy > Y_c:
        dist_to_center = np.sqrt((vx - X_c)**2 + (vy - Y_c)**2)
        is_magnet = (R_in <= dist_to_center <= R_out)
    else:
        in_y_range = (MAGNET_Y_MIN <= vy <= Y_c)
        in_left_leg = (X_c - R_out <= vx <= X_c - R_in)
        in_right_leg = (X_c + R_in <= vx <= X_c + R_out)
        is_magnet = in_y_range and (in_left_leg or in_right_leg)
    is_magnet_mask[i] = is_magnet

magnet_cells_mm = all_voxels_mm[is_magnet_mask]   # 仅磁铁实体的坐标（mm）
voxel_pts = magnet_cells_mm / 1000.0              # 转换为米，用于反演
num_voxels = len(voxel_pts)

print(f"总网格数: {mesh.nC}，磁铁实体数: {num_voxels}")

# ==========================================
# 5. 构建敏感度矩阵（与原代码一致）
# ==========================================
def build_matrix_optimized(meas_pts, voxels, vol):
    M, N = len(meas_pts), len(voxels)
    A = np.zeros((M * 3, N * 3))
    mu0_4pi = 1e-7
    coef = mu0_4pi * vol

    dr = meas_pts[:, None, :] - voxels[None, :, :]
    r_len = np.linalg.norm(dr, axis=2)
    r_len[r_len < 1e-4] = 1e-4

    r3 = r_len**3
    r5 = r_len**5
    dx, dy, dz = dr[:, :, 0], dr[:, :, 1], dr[:, :, 2]

    for col_idx, d_coord in enumerate([dx, dy, dz]):
        term1 = 3 * d_coord[:, :, None] * dr / r5[:, :, None]
        identity_term = np.zeros((3,))
        identity_term[col_idx] = 1.0
        term2 = identity_term[None, None, :] / r3[:, :, None]
        B_contrib = coef * (term1 - term2)

        for k in range(3):
            A[k::3, col_idx::3] = B_contrib[:, :, k]
    return A

print("正在计算敏感度矩阵...")
A = build_matrix_optimized(meas_pts, voxel_pts, voxel_volume)

# ==========================================
# 6. 最小二乘反演 + 方向平滑
# ==========================================
print("正在执行反演...")
lam = 1e-16
A_reg = np.vstack([A, np.sqrt(lam) * np.eye(3 * num_voxels)])
b_reg = np.concatenate([b_meas_vector, np.zeros(3 * num_voxels)])

max_M = 1.5e6
M_vector = np.linalg.lstsq(A_reg, b_reg, rcond=None)[0]
M_vector = np.clip(M_vector, -max_M, max_M)
M_estimated = M_vector.reshape((num_voxels, 3))

# 邻域方向平滑
print("正在执行方向平滑...")
M_smoothed = M_estimated.copy()
smooth_radius_mm = 12.0
for i, r_v1 in enumerate(voxel_pts):
    dists = np.linalg.norm(voxel_pts - r_v1, axis=1) * 1000.0
    neighbors = np.where(dists <= smooth_radius_mm)[0]
    if len(neighbors) > 0:
        M_smoothed[i] = np.mean(M_estimated[neighbors], axis=0)
M_estimated = M_smoothed
print("反演完成！")

# ==========================================
# 7. 重组完整网格并导出 CSV
# ==========================================
print(f"正在保存结果到 {output_m_csv}...")

# 全网格磁化向量（初始为0）
full_M = np.zeros((mesh.nC, 3))
# 将反演结果填入磁铁位置
full_M[is_magnet_mask] = M_estimated

# 导出（坐标用 mm 便于原 VTK 脚本）
output_df = pd.DataFrame({
    'x': all_voxels_mm[:, 0],
    'y': all_voxels_mm[:, 1],
    'z': all_voxels_mm[:, 2],
    'Mx': full_M[:, 0],
    'My': full_M[:, 1],
    'Mz': full_M[:, 2]
})
output_df.to_csv(output_m_csv, index=False)
print("✅ CSV 文件保存成功！")