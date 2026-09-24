# %% [markdown]
# # 磁场数据与磁化强度的二维截面
#
# 本文生成文章配图所需的两张二维截面图：
# 1. 测量磁场数据在扫描平面（z = 0）上的分布；
# 2. 反演磁化强度在磁铁中面（z ≈ -19 mm）上的分布。
#
# 数据全部来自真实测量/反演结果，不做示意数据。

# %%
# ============================================
# 导入与绘图配置
# ============================================
import os
import numpy as np
import pandas as pd
import glob
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
import pyvista as pv

# 中文标注
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

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
FIGS.mkdir(parents=True, exist_ok=True)

# %% [markdown]
# # 磁铁几何（与 B2M.ipynb 一致）
#
# U 形磁铁（UshapeNormal）：X 111.5–193.5、Y 129.5–231.5、Z -35–-3（单位 mm）。
# 圆弧中心 (X_c, Y_c)，外半径 R_out、内半径 R_in。

# %%
MAGNET_X_MIN, MAGNET_X_MAX = 111.5, 193.5
MAGNET_Y_MIN, MAGNET_Y_MAX = 129.5, 231.5
MAGNET_Z_MIN, MAGNET_Z_MAX = -35.0, -3.0
X_c = (MAGNET_X_MIN + MAGNET_X_MAX) / 2.0
Y_c = MAGNET_Y_MAX - (MAGNET_X_MAX - MAGNET_X_MIN) / 2.0
R_out = (MAGNET_X_MAX - MAGNET_X_MIN) / 2.0
R_in = 22.0


def draw_u_outline(ax, color="k", lw=1.2):
    """在 x-y 平面绘制 U 形磁铁轮廓（外弧 + 内弧 + 两条磁腿）。"""
    theta = np.linspace(0.0, np.pi, 200)
    ax.plot(X_c + R_out * np.cos(theta), Y_c + R_out * np.sin(theta), color=color, lw=lw)
    ax.plot(X_c + R_in * np.cos(theta), Y_c + R_in * np.sin(theta), color=color, lw=lw)
    for x0, x1 in ((X_c - R_out, X_c - R_in), (X_c + R_in, X_c + R_out)):
        ax.plot([x0, x0], [Y_c, MAGNET_Y_MIN], color=color, lw=lw)
        ax.plot([x1, x1], [Y_c, MAGNET_Y_MIN], color=color, lw=lw)
        ax.plot([x0, x1], [MAGNET_Y_MIN, MAGNET_Y_MIN], color=color, lw=lw)

# %% [markdown]
# # 图 A：磁场数据二维截面（z = 0 扫描平面）
#
# `data/raw/UshapeNormal/Large.csv`（圆弧区，4 mm 间距）与 `data/raw/UshapeNormal/legs.csv`（磁腿区，2 mm 间距）
# 在 y=150/152 处衔接，共同覆盖完整 U 形扫描区。取 z = 0（传感器最低扫描面），
# 用线性插值绘制到统一规则网格。

# %%
csv_files = glob.glob(str(RAW / "UshapeNormal" / "*.csv"))
df = pd.concat([pd.read_csv(f) for f in csv_files], ignore_index=True)
d0 = df[df["z"] == 0]
x, y = d0["x"].values, d0["y"].values
Bx, By, Bz = d0["Bx"].values, d0["By"].values, d0["Bz"].values
Bmag = np.sqrt(Bx**2 + By**2 + Bz**2)  # μT

gi = np.linspace(x.min(), x.max(), 400)
gj = np.linspace(y.min(), y.max(), 500)
GX, GY = np.meshgrid(gi, gj)
Bgrid = griddata((x, y), Bmag, (GX, GY), method="linear")

fig, ax = plt.subplots(figsize=(6.6, 5.4))
im = ax.pcolormesh(GX, GY, Bgrid, cmap="viridis", shading="auto")
cb = fig.colorbar(im, ax=ax, pad=0.02)
cb.set_label("磁场强度 |B| / μT")

# 面内方向（下采样，单位化后显示方向）
xs, ys = x[::4], y[::4]
ux, uy = Bx[::4], By[::4]
um = np.sqrt(ux**2 + uy**2)
um[um == 0] = 1.0
ax.quiver(xs, ys, ux / um, uy / um, color="white", scale=30.0, width=0.004,
          headwidth=3.0, headlength=4.0, alpha=0.85)

draw_u_outline(ax)
ax.set_xlabel("x / mm")
ax.set_ylabel("y / mm")
ax.set_title("测量磁场 |B| 分布（z = 0 平面）")
ax.set_aspect("equal")
fig.tight_layout()
fig.savefig(FIGS / "B_field_section.png", dpi=300)
fig.savefig(FIGS / "B_field_section.pdf")
plt.close(fig)
print(f"已保存 {FIGS}/B_field_section.png / .pdf")

# %% [markdown]
# # 图 B：磁化强度二维截面（z ≈ -19 mm 磁铁中面）
#
# 直接读取 B2M.ipynb 导出的正则化反演结果 `data/derived/UshapeNormal/inverted_M_results_test.vtu`（TreeMesh），
# 取 z 接近 -19 mm（Z 镜像对称中面）的磁铁单元，插值绘制。

# %%
mesh_vtu = pv.read(str(DERIVED / "UshapeNormal" / "inverted_M_results_test.vtu"))
cc = mesh_vtu.cell_centers().points  # (nC, 3)，单位 mm
Mx = mesh_vtu.cell_data["Mx"]
My = mesh_vtu.cell_data["My"]
Mz = mesh_vtu.cell_data["Mz"]
Mmag = np.sqrt(Mx**2 + My**2 + Mz**2)  # A/m

# 磁铁单元 z 中心为 2 mm 间距的偶数层（…-20、-18…），对称中面 z=-19 位于两层之间。
# 自动选取最接近 -19 的那一层（镜像对称，上下两层 |M| 相同）。
z_target = -19.0
z_mag = cc[Mmag > 1.0, 2]
z_layers = np.unique(np.round(z_mag, 1))
z_layer = z_layers[np.argmin(np.abs(z_layers - z_target))]
sel = np.isclose(cc[:, 2], z_layer, atol=0.15) & (Mmag > 1.0)  # 仅磁铁内部
xc, yc = cc[sel, 0], cc[sel, 1]
mc = Mmag[sel]
mxc, myc = Mx[sel], My[sel]
print(f"z = {z_layer} mm 截面（对称中面附近）：{sel.sum()} 个磁铁单元，|M| 范围 [{mc.min():.1f}, {mc.max():.1f}] A/m")

gi2 = np.linspace(xc.min(), xc.max(), 400)
gj2 = np.linspace(yc.min(), yc.max(), 500)
GX2, GY2 = np.meshgrid(gi2, gj2)
Mgrid = griddata((xc, yc), mc, (GX2, GY2), method="linear")

# 仅磁铁实体内部填色：U 形内凹区域与外部置 NaN（透明），避免插值把 U 形内空气区也涂色
dist2 = np.sqrt((GX2 - X_c) ** 2 + (GY2 - Y_c) ** 2)
in_arc = (GY2 > Y_c) & (R_in <= dist2) & (dist2 <= R_out)
in_legs = (GY2 <= Y_c) & (GY2 >= MAGNET_Y_MIN) & (
    ((GX2 >= X_c - R_out) & (GX2 <= X_c - R_in)) |
    ((GX2 >= X_c + R_in) & (GX2 <= X_c + R_out))
)
magnet_mask = in_arc | in_legs
Mgrid = np.where(magnet_mask, Mgrid, np.nan)

fig2, ax2 = plt.subplots(figsize=(6.6, 5.4))
im2 = ax2.pcolormesh(GX2, GY2, Mgrid / 1000.0, cmap="inferno", shading="auto")
cb2 = fig2.colorbar(im2, ax=ax2, pad=0.02)
cb2.set_label("磁化强度 |M| / (kA/m)")

xc_s, yc_s = xc[::4], yc[::4]
mux, muy = mxc[::4], myc[::4]
mum = np.sqrt(mux**2 + muy**2)
mum[mum == 0] = 1.0
ax2.quiver(xc_s, yc_s, mux / mum, muy / mum, color="white", scale=22.0, width=0.004,
           headwidth=3.0, headlength=4.0, alpha=0.85)

draw_u_outline(ax2)
ax2.set_xlabel("x / mm")
ax2.set_ylabel("y / mm")
ax2.set_title(f"反演磁化强度 |M| 分布（z = {z_layer:g} mm 截面）")
ax2.set_aspect("equal")
fig2.tight_layout()
fig2.savefig(FIGS / "M_section.png", dpi=300)
fig2.savefig(FIGS / "M_section.pdf")
plt.close(fig2)
print(f"已保存 {FIGS}/M_section.png / .pdf")
