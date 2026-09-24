# -*- coding: utf-8 -*-
"""meshtest 共用核心：网格 / 正演核 / 正则化变体 / 目标函数 / 还原度度量。

设计原则
--------
1. **与生产代码逐行同构**。正演 numba 核、Huber(Charbonnier) 正则化、L-BFGS-B 目标函数
   都从 B2M_newprob.py 原样搬过来；E0 实验会用"逐字复制生产公式"的函数与本模块对照，
   要求机器精度一致。任何"harness 与生产不一致"都会让后续结论失效。
2. **正则化变体集中在一处**，靠 ``variant`` 切换，避免手改公式造成标度漂移：

   ==============  ================================================================
   variant         公式（w_f = 面对偶体积 S·h；W = Σ w_f ≈ 3·V_magnet）
   ==============  ================================================================
   none            R = λ Σ_f ρ(|GM|_f)                      旧写法，每面权重 1
   volume_raw      R = λ Σ_f w_f ρ(|GM|_f)                  w 用网格单位(mm³)
   volume_norm     R = λ Σ_f (w_f/W) ρ(|GM|_f)              **推荐**：加密/单位都不变
   volume_m3       R = λ Σ_f (w_f/1e9) ρ(|GM|_f)            m³ 写法(专门复现 1e9 陷阱)
   ==============  ================================================================

   其中 ρ(x) = ε²(√(1+(x/ε)²) − 1)，ε = huber_epsilon。

3. 单位约定与生产一致：几何用 mm，正演时体素体积换算成 m³，正则化权重用网格原始单位(mm³)。
"""

from __future__ import annotations

import numpy as np
from numba import njit, prange
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import LinearOperator
from scipy.optimize import minimize

__all__ = [
    "LEVEL_NAMES", "REG_VARIANTS",
    "matvec", "rmatvec",
    "make_forward", "build_mesh", "magnet_cells", "magnet_levels",
    "VARIABLES", "variable_scale", "operator_coef", "to_m", "check_gradient",
    "build_regularizer", "chi", "dchi",
    "objective_and_gradient", "production_objective", "solve",
    "metrics", "volume_weighted_norm", "describe_mesh",
    "data_term", "reg_term", "calibrate_lambda", "make_reference_model",
]

REG_VARIANTS = ("none", "volume_raw", "volume_norm", "volume_m3")

LEVEL_NAMES = {0: "最细", 1: "次细", 2: "粗两级", 3: "粗三级"}


# --------------------------------------------------------------------------------------
# 正演核：逐字复制自 B2M_newprob.py（勿改，改了 E0 会报错）
# --------------------------------------------------------------------------------------
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


def make_forward(meas_m, voxel_points_m, coef):
    """返回 (matvec_f, rmatvec_f)：测点/体素都用米。

    第三个参数是**算子系数**，由调用方按所选的求解变量决定（见 ``operator_coef``）：
    优化 m 时传体积 $v_{\\rm SI}$；优化 $x=\\sqrt{v_{\\rm SI}}\\,m$ 时传 $\\sqrt{v_{\\rm SI}}$。
    两种情况下都必须与回代换算（``variable_scale``）成对使用，否则幅值会差一个常数因子。
    """
    meas_m = np.ascontiguousarray(meas_m, dtype=np.float64)
    voxel_points_m = np.ascontiguousarray(voxel_points_m, dtype=np.float64)
    coef = np.ascontiguousarray(coef, dtype=np.float64)

    def mv(v):
        return matvec(np.ascontiguousarray(v, dtype=np.float64),
                      meas_m, voxel_points_m, coef)

    def rv(w):
        return rmatvec(np.ascontiguousarray(w, dtype=np.float64),
                       meas_m, voxel_points_m, coef)

    return mv, rv


# --------------------------------------------------------------------------------------
# 求解变量：x = scale * m
#   "m"         : x = m            （历史写法）
#   "sqrt_v_m"  : x = √(v_SI) · m  （E7 采用的写法：求解器隐式最小化的 Σx² 正好是 ∫m²dV）
# 两条不变量（都踩过坑）：
#   I1 正演系数与回代换算必须成对：coef=√v_SI 配 m=x/√v_SI；且 √v 必须用 SI（m³）
#   I2 梯度里只有正则化项要乘 ∂m/∂x = 1/scale，数据项不能再除一次
# --------------------------------------------------------------------------------------
VARIABLES = ("m", "sqrt_v_m")


def variable_scale(voxel_vol_m3, variable="m"):
    """x = scale · m 中的 scale。$\\sqrt{v}$ 用 SI（输入体积须为 m³）。"""
    v = np.asarray(voxel_vol_m3, float)
    if variable == "m":
        return np.ones_like(v)
    if variable == "sqrt_v_m":
        return np.sqrt(v)
    raise ValueError(f"未知 variable {variable!r}（可选 {VARIABLES}）")


def operator_coef(voxel_vol_m3, variable="m"):
    """传给正演核的系数（与 ``variable_scale`` 成对）。"""
    v = np.asarray(voxel_vol_m3, float)
    if variable == "m":
        return v
    if variable == "sqrt_v_m":
        return np.sqrt(v)
    raise ValueError(f"未知 variable {variable!r}（可选 {VARIABLES}）")


def to_m(x, scale):
    """变量 x → 磁化强度 m = x / scale。"""
    s = np.asarray(scale, float)
    return np.asarray(x, float).reshape(-1, 3) / s[:, None]


def check_gradient(x, b, matvec_f, rmatvec_f, G_sub, w, eps, lam, scale=None,
                   n_probe=12, h=None, seed=0):
    """中心差分校验解析梯度，返回**抽样分量上的相对范数误差** ‖fd−g‖/‖g‖（不变量 I2 的回归测试）。

    用相对范数而不是逐分量相对误差：后者的分母可能接近 0（梯度各分量量级差很大时），
    会把浮点噪声放大成假失败。

    ⚠ 探测点要用**非零**的 x：在 x=0 处 Charbonnier 的导数恒为 0，λ 项不参与，
    测不出"数据项被多除了一次 scale"这个 bug（正是踩过的那个）。
    """
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(x), size=min(n_probe, len(x)), replace=False)
    if h is None:
        # 差分步长按变量自身的量级取（两种变量的尺度差 10⁴ 量级，固定步长会一边太粗一边太细）
        # 步长按变量量级取：本问题的 ‖梯度‖ 很小（~1e-9），步长太小时中心差分会
        # 被 E 的浮点舍入淹没（相对误差 ~ eps·E/(h·‖g‖)）。取 1e-3·|x| 兼顾两者。
        h = 1e-3 * max(float(np.mean(np.abs(x))), 1e-30)

    def f(v):
        return objective_and_gradient(v, b, matvec_f, rmatvec_f, G_sub, w, eps, lam,
                                      scale=scale)[0]

    _, g = objective_and_gradient(x, b, matvec_f, rmatvec_f, G_sub, w, eps, lam, scale=scale)
    fd = np.empty(len(idx))
    for k, i in enumerate(idx):
        xp = x.copy(); xp[i] += h
        xm = x.copy(); xm[i] -= h
        fd[k] = (f(xp) - f(xm)) / (2 * h)
    denom = max(float(np.linalg.norm(g[idx])), 1e-30)
    return float(np.linalg.norm(fd - g[idx]) / denom)


# --------------------------------------------------------------------------------------
# 网格与几何
# --------------------------------------------------------------------------------------
def _next_pow2(n):
    return 2 ** int(np.ceil(np.log2(n)))


def build_mesh(domain, nbase, style="split", magnet_box=None, split_axis=0,
               split_frac=0.5, mag_levels=1, diagonal_balance=False, fine_side="low",
               verbose=False):
    """构造实验网格。

    Parameters
    ----------
    domain : (Lx, Ly, Lz)，原点取 (0,0,0)，单位 mm
    nbase  : (nx, ny, nz) 基础网格单元数（会被向上取到 2 的幂）
    style  : "uniform"  全网格同一档
             "split"    split_axis 方向上 split_frac 之前细一档、之后粗一档
             "checker"  按 x/y 的 2×2 块交替细/粗（检验各向与交错界面）
             "production" 磁体 bbox 统一细化到次细档，磁体内部再按 split 分档
             （后三种都制造 1 档 = 线尺寸 2× = 体积 8× 的密度对比）
    fine_side : "low"（默认，低坐标一侧细）或 "high"（高坐标一侧细）；
             用来造出"互为粗细互换"的两张网格（E7 的互换不变性测试）
    mag_levels : 最高档与次高档之间的档数差（1 表示差一档）
    """
    from discretize import TreeMesh

    if fine_side not in ("low", "high"):
        raise ValueError(f"fine_side 只能取 'low'/'high'，收到 {fine_side!r}")

    Lx, Ly, Lz = domain
    nx, ny, nz = (_next_pow2(int(np.ceil(n))) for n in nbase)
    hx, hy, hz = Lx / nx, Ly / ny, Lz / nz
    mesh = TreeMesh([[(hx, nx)], [(hy, ny)], [(hz, nz)]],
                    origin=[0.0, 0.0, 0.0], diagonal_balance=diagonal_balance)
    ml = mesh.max_level
    fine, coarse = ml, ml - mag_levels

    box = np.array(magnet_box, dtype=float) if magnet_box is not None else \
        np.array([[0.0, 0.0, 0.0], [Lx, Ly, Lz]])

    def pick(c, ax, thr):
        """低坐标一侧细还是高坐标一侧细，由 fine_side 决定。"""
        return fine if ((c[ax] < thr) == (fine_side == "low")) else coarse

    if style == "uniform":
        mesh.refine_bounding_box(box, level=fine, finalize=False)
    elif style in ("split", "checker", "production"):
        mesh.refine_bounding_box(box, level=coarse, finalize=False)
        ax = split_axis
        thr = split_frac * domain[ax]

        def fn(cell):
            c = cell.center
            inside = all(box[0][k] - 1e-9 <= c[k] <= box[1][k] + 1e-9 for k in range(3))
            if not inside:
                return 0
            if style == "split":
                return pick(c, ax, thr)
            if style == "checker":
                i = int(np.floor(c[0] / (2 * hx)))
                j = int(np.floor(c[1] / (2 * hy)))
                return fine if (i + j) % 2 == 0 else coarse
            # production：split_axis 前半细、后半粗（与 B2M 的 y<Y_split 同构）
            return pick(c, ax, thr)

        mesh.refine(fn, finalize=False)
    else:
        raise ValueError(f"未知 style: {style}")

    mesh.finalize()
    mesh.number()
    if verbose:
        print(describe_mesh(mesh))
    return mesh


def magnet_cells(mesh, magnet_box):
    box = np.asarray(magnet_box, float)
    cc = mesh.cell_centers
    inside = np.all((cc >= box[0] - 1e-9) & (cc <= box[1] + 1e-9), axis=1)
    return np.where(inside)[0]


def magnet_levels(cell_vol):
    """按体积把单元分档，返回**序号**：0 = 最细（体积最小），1,2,... 依次更粗。

    用序号而不是 log2 比值，避免出现 0/3 这种看着像 bug 的标签。
    """
    v = np.asarray(cell_vol, float)
    u = np.unique(np.round(v, 9))
    idx = {val: i for i, val in enumerate(u)}
    return np.array([idx[val] for val in np.round(v, 9)], dtype=int)


def describe_mesh(mesh):
    v = mesh.cell_volumes
    u = np.unique(np.round(v, 6))
    return (f"TreeMesh nC={mesh.nC} nF={mesh.nF} max_level={mesh.max_level} "
            f"单元体积档位={len(u)} 档 {np.round(u, 3)} mm³")


# --------------------------------------------------------------------------------------
# 正则化
# --------------------------------------------------------------------------------------
def chi(x, eps):
    """Charbonnier / 平滑 Huber：eps²(√(1+(x/eps)²) − 1)"""
    return eps ** 2 * (np.sqrt(1.0 + (x / eps) ** 2) - 1.0)


def dchi(x, eps):
    """dχ/dx 中与 x 相乘的权重：1/√(1+(x/eps)²)（使 Gᵀ(w·g·wt) 与 dχ 一致）"""
    return 1.0 / np.sqrt(1.0 + (x / eps) ** 2)


def build_regularizer(mesh, in_magnet, variant="volume_norm"):
    """构造梯度算子与权重。

    Returns
    -------
    dict(G_sub, w, W, V_magnet, n_faces, cell_vol)
      G_sub     : (n_internal_faces, n_magnet_cells) 弱形式梯度
      w         : 每个内部面的**最终**权重（已按 variant 归一化）
      W         : w 的原始总量（Σ 未归一的面对偶体积），用于报告与换算
      V_magnet  : 磁体单元体积之和（mm³）
    """
    if variant not in REG_VARIANTS:
        raise ValueError(f"未知 variant: {variant}")
    G = mesh.cell_gradient
    if not isinstance(G, csr_matrix):
        G = G.tocsr()
    w_all = mesh.get_face_inner_product().diagonal()      # 面对偶体积 S·h（mm³）

    is_mag = np.asarray(in_magnet, bool)
    n_adj = np.diff(G.indptr)
    rows = np.where(n_adj >= 2)[0]
    # 只保留"相邻单元全部在磁体内"的内部面
    keep = np.array([bool(is_mag[G.indices[G.indptr[f]:G.indptr[f + 1]]].all())
                     for f in rows], dtype=bool)
    faces = rows[keep]

    magnet_idx = np.where(is_mag)[0]
    G_sub = G[faces, :][:, magnet_idx]
    V_f = w_all[faces]                                    # mm³
    W = float(V_f.sum())
    V_magnet = float(mesh.cell_volumes[magnet_idx].sum())

    if variant == "none":
        w = np.ones_like(V_f)
    elif variant == "volume_raw":
        w = V_f.copy()
    elif variant == "volume_norm":
        w = V_f / W if W > 0 else V_f
    elif variant == "volume_m3":
        w = V_f / 1e9
    return dict(G_sub=G_sub, w=w, w_face=V_f, W=W, V_magnet=V_magnet,
                n_faces=len(faces), magnet_idx=magnet_idx)


# --------------------------------------------------------------------------------------
# 目标函数与求解
# --------------------------------------------------------------------------------------
def objective_and_gradient(M_vec, b, matvec_f, rmatvec_f, G_sub, w, eps, lam, scale=None):
    """与生产 notebook 同构：数据项 Σ r² + λ Σ_f w_f χ(|GM|_f)。

    ``scale`` 不为 None 时，传进来的向量是**求解变量** $x$（$m=x/\\text{scale}$）：

    * 数据项：正演算子已经按所选变量构造好（系数 = ``operator_coef``），
      所以梯度就是 $2A_x^{\\mathsf T}r$，**不能再除 scale**；
    * 正则化项：写成 $m$ 的梯度后再乘 $\\partial m/\\partial x=1/\\text{scale}$。

    这两条就是"梯度链式法则"的落点；实测把数据项也除一次会让 L-BFGS-B 第 0 次迭代
    就 `ABNORMAL`、解恒为 0（见 ``check_gradient`` 的回归测试）。
    """
    r = matvec_f(M_vec) - b
    if scale is None:
        M = M_vec.reshape(-1, 3)
        div = None
    else:
        M = to_m(M_vec, scale)
        div = np.repeat(np.asarray(scale, float), 3)
    g = G_sub @ M
    gn = np.sqrt(np.sum(g ** 2, axis=1))
    loss = float(np.sum(r ** 2) + lam * np.sum(w * chi(gn, eps)))
    wt = dchi(gn, eps)
    grad = 2.0 * rmatvec_f(r)
    if lam != 0.0:
        grad_reg = lam * (G_sub.T @ (w[:, None] * g * wt[:, None])).ravel()
        grad = grad + (grad_reg / div if div is not None else grad_reg)
    return loss, grad


def production_objective(M_vec, b, matvec_f, rmatvec_f, G_sub, V_face, eps, lam):
    """逐字复制 B2M_newprob.py / B2M.ipynb 第 5 节的写法（volume_raw + 原始 λ）。

    只用于 E0 一致性锚定：本模块的 volume_raw 必须与它机器精度一致。
    """
    r = matvec_f(M_vec) - b
    M = M_vec.reshape(-1, 3)
    grad = G_sub @ M
    grad_norm = np.sqrt(np.sum(grad ** 2, axis=1))
    loss_reg = lam * eps ** 2 * np.sum(V_face * (np.sqrt(1.0 + (grad_norm / eps) ** 2) - 1.0))
    weights = 1.0 / np.sqrt(1.0 + (grad_norm / eps) ** 2)
    grad_reg = lam * (G_sub.T @ (V_face[:, None] * grad * weights[:, None])).ravel()
    return float(np.sum(r ** 2) + loss_reg), 2.0 * rmatvec_f(r) + grad_reg


def solve(b, matvec_f, rmatvec_f, G_sub, w, eps, lam, n_voxels, max_M=1.5e6,
          maxiter=50, x0=None, ftol=1e-30, gtol=1e-30, scale=None):
    """L-BFGS-B 求解（与生产相同的求解器与边界）。

    ``scale`` 不为 None 时求解变量是 $x=\\text{scale}\\cdot m$，逐格边界取
    $|x_i|\\le$ ``max_M``$\\cdot\\text{scale}_i$（等价于物理上限 $|m|\\le$ ``max_M``）；
    返回的仍是 $x$ 的 (n,3) 形状，用 ``to_m`` 换算回磁化强度。
    """
    fg = lambda x: objective_and_gradient(x, b, matvec_f, rmatvec_f, G_sub, w, eps, lam,
                                          scale=scale)
    if scale is None:
        bounds = [(-max_M, max_M)] * (3 * n_voxels)
    else:
        s = np.asarray(scale, float)
        bounds = [(-max_M * si, max_M * si) for si in s for _ in range(3)]
    res = minimize(fg, np.zeros(3 * n_voxels) if x0 is None else x0, jac=True,
                   method="L-BFGS-B", bounds=bounds,
                   options={"maxiter": maxiter, "ftol": ftol, "gtol": gtol})
    return res.x.reshape(n_voxels, 3), res


# --------------------------------------------------------------------------------------
# 度量
# --------------------------------------------------------------------------------------
def volume_weighted_norm(field, cell_vol):
    """sqrt(Σ V_c |f_c|²)。field 可为 (n,3) 或 (n,)。"""
    f = np.asarray(field, float)
    v = np.asarray(cell_vol, float)
    if f.ndim == 1:
        return float(np.sqrt(np.sum(v * f ** 2)))
    return float(np.sqrt(np.sum(v[:, None] * f ** 2)))


def data_term(M_vec, b, matvec_f):
    """Σ r²（数据项，不含 λ）。"""
    r = matvec_f(M_vec) - b
    return float(np.sum(r ** 2))


def reg_term(M_vec, G_sub, w, eps):
    """Σ_f w_f χ(|GM|_f)（正则化项在 λ=1 时的值）。"""
    g = G_sub @ M_vec.reshape(-1, 3)
    gn = np.sqrt(np.sum(g ** 2, axis=1))
    return float(np.sum(w * chi(gn, eps)))


def calibrate_lambda(M_ref, b, matvec_f, rmatvec_f, G_sub, w, eps):
    """标定 λ：在参考模型处让「数据项梯度」与「正则化项梯度」等大。

    为什么不用 λ = 数据项/正则化项：对**分段常数真值**（均匀块）梯度恒为 0 ⇒ 正则化项 = 0，
    标定退化。改用**力的平衡**：λ_bal = ‖∇E_data(M_ref)‖ / ‖∇E_reg(M_ref)‖|_{λ=1}，
    这正是"数据想把 M 拉向真值"与"TV/Huber 摩擦"相抵的强度。深 L1 区里 λε 相当于与梯度
    大小无关的常数摩擦力：λ 过大会让 M 在到达真值前就被刹住（幅值整体收缩），λ 过小则
    网格依赖与噪声拟合暴露出来 —— 这两侧正是要扫出指纹的地方。

    参考模型由 ``make_reference_model`` 给出（M_true + 角尺度约 1/4 磁体尺寸的正弦扰动），
    对均匀块/接缝块两种真值都有良好且可复现的梯度尺度。
    """
    M_ref = np.asarray(M_ref, float)
    r = matvec_f(M_ref.ravel()) - b
    gd = 2.0 * rmatvec_f(r)
    g = G_sub @ M_ref.reshape(-1, 3)
    gn = np.sqrt(np.sum(g ** 2, axis=1))
    wt = dchi(gn, eps)
    gr = (G_sub.T @ (w[:, None] * g * wt[:, None])).ravel()
    ng_d, ng_r = float(np.linalg.norm(gd)), float(np.linalg.norm(gr))
    if ng_r <= 0.0:
        return float("inf")
    return ng_d / ng_r


def make_reference_model(points_mm, M_true=None, amp=None, amp_frac=0.3, waves=4.0):
    """参考模型：M_true + 沿 x 的正弦扰动（幅值 = amp_frac × 真值 rms）。

    生产里没有真值，可以只给 ``amp``（材料饱和磁化强度的量级猜测）构造纯正弦参考模型
    —— λ_bal 只需要一个"合理量级"的参考模型，不要求精确。这样同一套标定流程在
    合成实验与生产数据上完全一致。
    """
    p = np.asarray(points_mm, float)
    span = float(np.ptp(p[:, 0])) or 1.0
    if M_true is None:
        if amp is None:
            raise ValueError("M_true 与 amp 至少给一个")
        M_ref = np.zeros((len(p), 3))
    else:
        M_true = np.asarray(M_true, float)
        M_ref = M_true.copy()
        if amp is None:
            amp = amp_frac * float(np.sqrt(np.mean(np.sum(M_true ** 2, axis=1))))
    M_ref[:, 0] += amp * np.sin(2.0 * np.pi * waves * (p[:, 0] - p[:, 0].min()) / span)
    return M_ref


def metrics(M_rec, M_true, cell_vol, b, matvec_f, level=None, M_true_level=None):
    """还原度指标（全部体积加权）。

    rel_err     : ||M_rec − M_true||_V / ||M_true||_V   —— "还原度"主判据
    bias_level  : 每档 Σ V|M_rec| / Σ V|M_true|（=1 无偏）—— "密小疏大"的直接度量
    maxM_ratio  : max|M_rec| / max|M_true|
    data_rms    : 模型场与数据之差的 rms
    """
    M_rec = np.asarray(M_rec, float)
    M_true = np.asarray(M_true, float)
    cell_vol = np.asarray(cell_vol, float)
    out = {}
    out["rel_err"] = volume_weighted_norm(M_rec - M_true, cell_vol) / \
        max(volume_weighted_norm(M_true, cell_vol), 1e-30)
    m_rec = np.linalg.norm(M_rec, axis=1)
    m_true = np.linalg.norm(M_true, axis=1)
    out["maxM_ratio"] = float(m_rec.max() / max(m_true.max(), 1e-30))
    r = matvec_f(M_rec.ravel()) - b
    out["data_rms"] = float(np.sqrt(np.mean(r ** 2)))
    if level is not None:
        level = np.asarray(level)
        bias = {}
        for lv in np.unique(level):
            m = level == lv
            num = float(np.sum(cell_vol[m] * m_rec[m]))
            den = float(np.sum(cell_vol[m] * m_true[m]))
            bias[int(lv)] = num / max(den, 1e-30)
        out["bias_level"] = bias
    return out
