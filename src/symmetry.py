# -*- coding: utf-8 -*-
"""磁铁 Z 对称面自动检测 + 镜像重合数据处理（B2M 项目公用模块）。

背景
----
B2M 系列 notebook 用测量到的外部磁场反推磁铁内部磁化强度。历史做法是把测点沿 Z 向
镜像一份（z -> 2*Z_center - z, Bz -> -Bz）来补充数据、改善收敛，其中 Z_center 取
(MAGNET_Z_MIN + MAGNET_Z_MAX)/2，是手工给定的常数。

换成小探头之后，测点已经跨越了磁铁的对称面，镜像点会与实测点**坐标重合**（本仓库新数据
UshapeNormal_New_prob 在 h=8 mm 时重合 18405/48784 = 37.7%），既重复计权又互相矛盾。
本模块提供两件事：

1. ``estimate_symmetry_plane``：用测量数据自身自动定出对称面高度 h（不需要知道几何参数）。
   原理：对一块磁化强度位于 x-y 面内、沿 z 均匀的磁铁，绕厚度中面 z=h 有

       Bx(x,y,2h-z) = Bx(x,y,z)
       By(x,y,2h-z) = By(x,y,z)
       Bz(x,y,2h-z) = -Bz(x,y,z)

   （推导：把源按 ±z 成对配对，Mz=0 ⇒ M·dr 与 z 无关，故 Bx/By 偶、Bz 奇。注意这与
   "电流分布镜像不变"的情形符号相反，不能凭直觉写。）于是可以定义镜像一致性目标

       J(h) = mean_k( sum ΔB_k^2 / sum B_k^2 ),  ΔB = B(2h-z) - (Bx, By, -Bz)(z)

   在数据 z 跨度内扫描 h：J 取极小值的 h 就是对称面。J 是无量纲相对失配，J<<1 表示
   数据确实关于该平面镜像对称，J≈1 表示完全没有这种对称性（单侧数据即如此）。

2. ``apply_mirror``：按检测出的 h 做镜像，并处理重合点。默认 policy="keep_measured"，
   即"实测值永远优先"：镜像点若落在已有实测点附近就直接丢弃，只在实测未覆盖的区域
   保留镜像补点。可选 policy="average" 则把重合对做对称化平均
   B <- (B(p) + R·B(sigma(p)))/2（对精确重合同等于两条数据各 0.5 权，不需要给目标函数
   加权重向量），并把剩余镜像点补进空区。

用法
----
    python src/symmetry.py data/raw/UshapeNormal_New_prob   # 检测 + 报告 + 出图
    python src/symmetry.py data/raw/UshapeNormal_New_prob --policy average
    python src/symmetry.py --selftest                          # 合成数据自测，不需要外部文件

作为库使用：

    from symmetry import load_measurements, estimate_symmetry_plane, apply_mirror
    P, B, meta = load_measurements("data/raw/UshapeNormal_New_prob/*.csv")
    res = estimate_symmetry_plane(P, B)
    P2, B2, rep = apply_mirror(P, B, res.h)
"""

from __future__ import annotations

import argparse
import glob
import os
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.interpolate import CubicSpline
from scipy.spatial import cKDTree

__all__ = [
    "SymmetryResult",
    "MirrorReport",
    "load_measurements",
    "estimate_symmetry_plane",
    "estimate_noise",
    "apply_mirror",
    "default_tol",
    "legacy_mirror",
    "regression_check",
    "plot_symmetry_diagnostics",
    "plot_mirror_overlap",
    "selftest",
]

# 物理约定：镜像把 Bz 取反，Bx/By 保持不变
MIRROR_SIGN = np.array([1.0, 1.0, -1.0])

# 检测器默认参数（理由见 estimate_symmetry_plane 的 docstring）
DEFAULTS = dict(
    dz_coarse=0.5,      # 粗扫重采样网格步长 [mm]
    dz_fine=0.25,       # 精扫重采样网格步长 [mm]
    step_coarse=0.25,   # 粗扫 h 步长 [mm]
    step_fine=0.05,     # 精扫 h 步长 [mm]
    coarse_cols=800,    # 粗扫使用的列数（随机抽样，控制耗时）
    min_col_pts=5,      # 一条 z 列至少多少个测点
    min_col_span=6.0,   # 一条 z 列至少多长 [mm]
    min_pairs=200,      # 有效配对数下限
    min_cols=50,        # 有效列数下限
    j_max=0.05,         # 接受判据：J(h*) <= j_max
    j_ratio_max=0.2,    # 接受判据：J(h*) <= j_ratio_max * median(J)
    band_max=2.0,       # 接受判据：J<1.1*Jmin 的带宽 <= band_max [mm]
)


# --------------------------------------------------------------------------------------
# 数据结构
# --------------------------------------------------------------------------------------
@dataclass
class SymmetryResult:
    """对称面检测结果。"""

    h: float                      # 最终采用的对称面高度 [mm]
    h_raw: float                  # 精扫得到的最优值 [mm]
    h_snap: float | None          # 吸附到"镜像落回测量 z 网格"的候选 [mm]
    j_min: float                  # J(h_raw)
    band: tuple[float, float]     # J < 1.1*J_min 的 h 区间
    n_pairs: int                  # h_raw 处的有效配对数
    n_cols: int                   # 有效列数
    accepted: bool                # 是否接受自动检测结果
    reason: str                   # 判定说明
    curve: np.ndarray             # (n, 3): h, J, n_pairs（粗+精扫曲线）
    exact_table: list = field(default_factory=list)  # 纯网格配对交叉校验 [(h, rms_per_comp, n)]
    j_median: float = float("nan")
    n_used_cols: int = 0          # 参与检测的总列数
    comp_ratio: np.ndarray = field(default_factory=lambda: np.full(3, np.nan))  # h* 处各分量相对失配
    comp_snr: np.ndarray = field(default_factory=lambda: np.full(3, np.nan))    # 各分量信噪比粗估

    @property
    def band_width(self) -> float:
        return float(self.band[1] - self.band[0])

    def summary(self) -> str:
        lines = [
            "================ 对称面自动检测 ================",
            f"  检测结果 h_raw = {self.h_raw:.3f} mm   (J = {self.j_min:.4e})",
            f"  不确定带宽  = [{self.band[0]:.2f}, {self.band[1]:.2f}] mm "
            f"(宽 {self.band_width:.2f} mm, J<1.1*Jmin)",
        ]
        if self.h_snap is not None:
            lines.append(f"  网格吸附候选 h_snap = {self.h_snap:.3f} mm "
                         f"(|Δ| = {abs(self.h_raw - self.h_snap):.3f} mm)")
        lines += [
            f"  配对数 = {self.n_pairs}  有效列 = {self.n_cols}  总列 = {self.n_used_cols}",
            f"  J 中位数(全扫描范围) = {self.j_median:.4e}",
        ]
        if np.isfinite(self.comp_ratio).any():
            lines.append("  h* 处各分量相对失配 J_k (越小越对称, ~1 表示该分量无信息):")
            lines.append("    " + "  ".join(
                f"{n}={v:.3e}" for n, v in zip("Bx By Bz".split(), self.comp_ratio)))
        if np.isfinite(self.comp_snr).any():
            lines.append("  各分量信噪比粗估 (rms(B)/二阶差分噪声):")
            lines.append("    " + "  ".join(
                f"{n}={v:.0f}" for n, v in zip("Bx By Bz".split(), self.comp_snr)))
        lines += [
            f"  采用 h = {self.h:.3f} mm   accepted = {self.accepted}",
            f"  判定: {self.reason}",
        ]
        if self.exact_table:
            lines.append("  纯网格配对交叉校验（无插值误差）:")
            for h, rms, n in self.exact_table:
                lines.append(f"    h = {h:6.2f} mm  n = {n:6d}  "
                             f"rms|ΔB| = [{rms[0]:8.1f} {rms[1]:8.1f} {rms[2]:8.1f}] uT")
        lines.append("===============================================")
        return "\n".join(lines)


@dataclass
class MirrorReport:
    """镜像 + 去重结果统计。"""

    n_measured: int
    n_mirrored_total: int
    n_overlap: int              # 镜像点落在实测点 tol 内（重合）的数量
    n_fill: int                 # 保留下来的镜像补点数
    n_final: int
    tol: float
    policy: str
    fill: bool
    h: float
    asym_rms: np.ndarray        # 重合处 |B_measured - R*B_mirror| 的 rms（3 分量）[uT]
    asym_mean: np.ndarray       # 同上 mean
    asym_rel: np.ndarray        # 相对该分量的全局 rms
    b_rms: np.ndarray           # 重合点实测 B 的全局 rms（3 分量）[uT]
    worst: list = field(default_factory=list)   # [(x, y, z, dBx, dBy, dBz)]
    per_source: dict = field(default_factory=dict)
    n_duplicate_measured: int = 0   # 实测集内部坐标重复数（load_measurements 报告）
    noise_rms: np.ndarray = field(default_factory=lambda: np.full(3, np.nan))  # 二阶差分噪声估计 [uT]

    def summary(self) -> str:
        r = self.asym_rms
        lines = [
            "================ 镜像 / 重合处理 ================",
            f"  对称面 h = {self.h:.3f} mm   策略 = {self.policy}   补空区 = {self.fill}",
            f"  去重容差 tol = {self.tol:.3f} mm",
            f"  实测点 {self.n_measured}  镜像候选 {self.n_mirrored_total}"
            f"  其中重合 {self.n_overlap}  补点 {self.n_fill}",
            f"  最终点数 {self.n_final}",
        ]
        if self.n_duplicate_measured:
            lines.append(f"  注意: 实测集内部有 {self.n_duplicate_measured} 个重复坐标"
                         f"（默认只报警、不修改，以保持旧结果可复现）")
        if self.n_overlap:
            lines.append("  重合处不一致量 (实测 vs 镜像):")
            for k, name in enumerate("Bx By Bz".split()):
                extra = ""
                if np.isfinite(self.noise_rms[k]) and self.noise_rms[k] > 0:
                    extra = f"  = {r[k] / self.noise_rms[k]:5.1f}x 噪声({self.noise_rms[k]:.0f} uT)"
                lines.append(f"    {name}: rms = {r[k]:8.1f} uT  mean = {self.asym_mean[k]:8.1f} uT"
                             f"  (相对该分量 rms {self.asym_rel[k] * 100:6.2f}%){extra}")
            if self.per_source:
                lines.append("  按来源文件分组 (rms |ΔB| / uT):")
                for name, v in self.per_source.items():
                    lines.append(f"    {name:32s} n={v[0]:6d}  "
                                 f"[{v[1][0]:8.1f} {v[1][1]:8.1f} {v[1][2]:8.1f}]")
            if self.worst:
                lines.append("  最差 10 个重合点:")
                for x, y, z, dbx, dby, dbz in self.worst:
                    lines.append(f"    (x={x:7.2f}, y={y:7.2f}, z={z:6.2f})  "
                                 f"ΔB = [{dbx:9.1f} {dby:9.1f} {dbz:9.1f}] uT")
        lines.append("================================================")
        return "\n".join(lines)


# --------------------------------------------------------------------------------------
# 数据读取
# --------------------------------------------------------------------------------------
def load_measurements(pattern: str, dedup: str = "warn", verbose: bool = True):
    """读取一个或多个 csv（列: x,y,z,Bx,By,Bz），返回 (P[M,3] mm, B[M,3] uT, meta)。

    Parameters
    ----------
    pattern : glob 模式，例如 "data/raw/UshapeNormal_New_prob/*.csv"
    dedup   : "warn"（默认，只统计并告警，不修改数据，保证旧结果可复现）
              "merge"（把坐标重复的点做算术平均）
              "off"（什么都不做）
    """
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"没有匹配到任何 csv: {pattern}")
    frames = [pd.read_csv(f) for f in files]
    file_sizes = [len(fr) for fr in frames]
    df = pd.concat(frames, ignore_index=True)
    need = {"x", "y", "z", "Bx", "By", "Bz"}
    missing = need - set(df.columns)
    if missing:
        raise ValueError(f"{pattern} 缺少列: {sorted(missing)}")

    if dedup == "off":
        out = df
        n_dup = 0
    else:
        key = list(zip(df["x"].to_numpy(), df["y"].to_numpy(), df["z"].to_numpy()))
        seen = {}
        dups = []
        for i, k in enumerate(key):
            if k in seen:
                dups.append(i)
            else:
                seen[k] = i
        n_dup = len(dups)
        if dedup == "merge" and dups:
            g = df.groupby(["x", "y", "z"], as_index=False)[["Bx", "By", "Bz"]].mean()
            out = g
        else:
            out = df
    P = out[["x", "y", "z"]].to_numpy(dtype=float)
    B = out[["Bx", "By", "Bz"]].to_numpy(dtype=float)
    meta = dict(files=files, file_sizes=file_sizes, n_dup=n_dup, n_raw=len(df), dedup=dedup)
    if verbose:
        print(f"读取 {len(files)} 个文件: {[os.path.basename(f) for f in files]}")
        print(f"  测点 {len(P)} 个    x[{P[:, 0].min():.1f},{P[:, 0].max():.1f}] "
              f"y[{P[:, 1].min():.1f},{P[:, 1].max():.1f}] z[{P[:, 2].min():.1f},{P[:, 2].max():.1f}] mm")
        if n_dup:
            print(f"  警告: 实测集内部有 {n_dup} 个重复 (x,y,z) 坐标"
                  f"（dedup='{dedup}'，默认不修改数据）")
    return P, B, meta


# --------------------------------------------------------------------------------------
# 对称面检测
# --------------------------------------------------------------------------------------
def _build_columns(P, B, min_pts, min_span):
    """按 (x,y) 分组为 z 列；列内先对重复 z 取平均，再按 z 排序。"""
    order = np.lexsort((P[:, 2], P[:, 1], P[:, 0]))
    P = P[order]
    B = B[order]
    cols = []
    start = 0
    n = len(P)
    while start < n:
        end = start + 1
        while end < n and P[end, 0] == P[start, 0] and P[end, 1] == P[start, 1]:
            end += 1
        zs = P[start:end, 2]
        Bs = B[start:end]
        # 同一列内的重复 z：取平均，避免 CubicSpline 报 "x must be strictly increasing"
        uz, inv = np.unique(zs, return_inverse=True)
        if len(uz) != len(zs):
            acc = np.zeros((len(uz), 3))
            cnt = np.zeros(len(uz))
            np.add.at(acc, inv, Bs)
            np.add.at(cnt, inv, 1.0)
            Bs = acc / cnt[:, None]
        if len(uz) >= min_pts and (uz[-1] - uz[0]) >= min_span:
            cols.append((P[start, 0], P[start, 1], uz, Bs))
        start = end
    return cols


def _component_snr(P, B):
    """各分量信噪比粗估 = rms(B) / 二阶差分噪声估计。"""
    noise = estimate_noise(P, B)
    rms = np.sqrt((B ** 2).mean(0))
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(noise > 0, rms / noise, np.inf)


def estimate_noise(P, B, min_pts=5):
    """用列内 z 方向的二阶差分估计测量噪声 (3 分量) [uT]。

    sigma^2 = sum( d2^2 ) / (6 * N)，d2 = B(z+h) - 2B(z) + B(z-h)，只统计等间距的 z 列。
    """
    cols = _build_columns(P, B, min_pts, 0.0)
    acc = np.zeros(3)
    n = 0
    for _x, _y, zs, Bs in cols:
        d = np.diff(zs)
        if d.size == 0 or d.mean() <= 0 or d.std() > 1e-3 * d.mean():
            continue
        d2 = Bs[2:] - 2.0 * Bs[1:-1] + Bs[:-2]
        acc += (d2 ** 2).sum(0)
        n += d2.shape[0]
    if n == 0:
        return np.full(3, np.nan)
    return np.sqrt(acc / (6.0 * n))


def _resample(cols, z_grid, dz):
    """把每条列用三次样条重采样到统一 z 网格上，返回 (V[n,K,3], spans[n,2])，网格外为 NaN。"""
    V = np.full((len(cols), len(z_grid), 3), np.nan)
    spans = np.zeros((len(cols), 2))
    for i, (_x, _y, zs, Bs) in enumerate(cols):
        spans[i] = (zs[0], zs[-1])
        m = (z_grid >= zs[0] - 1e-9) & (z_grid <= zs[-1] + 1e-9)
        for k in range(3):
            V[i, m, k] = CubicSpline(zs, Bs[:, k])(z_grid[m])
    return V, spans


def _evaluator(z_grid, V, spans):
    """返回 J(h) 的可调用对象（矢量化：一次算完所有列）。"""
    z0 = spans[:, 0][:, None]
    z1 = spans[:, 1][:, None]
    Z = z_grid[None, :]
    dz = z_grid[1] - z_grid[0]
    finite = ~np.isnan(V)

    def J(h, with_pairs=False):
        t = (2.0 * h - z_grid - z_grid[0]) / dz
        i0 = np.clip(np.floor(t).astype(int), 0, len(z_grid) - 2)
        w = (t - i0)[None, :, None]
        Vm = V[:, i0, :] * (1.0 - w) + V[:, i0 + 1, :] * w
        valid = (~np.isnan(Vm)).all(2) & finite.all(2)
        tw = t[None, :]
        # 有效窗口：镜像采样点必须落在全局网格内（避免端部外插），
        # 且 z 与 2h-z 都落在该列自身的 z 跨度内（保证每对只统计一次）
        window = ((tw >= 0.0) & (tw <= len(z_grid) - 1.0)
                  & (Z >= np.maximum(z0, 2.0 * h - z1)) & (Z <= np.minimum(z1, 2.0 * h - z0)))
        m = valid & window
        n = int(m.sum())
        if n < 200:
            out = (np.nan, 0, np.full(3, np.nan))
            return out if with_pairs else out[0]
        A = Vm[m]
        T = V[m] * MIRROR_SIGN          # 目标值 (Bx, By, -Bz)
        num = (A - T) ** 2
        den = T ** 2
        ratio = num.sum(0) / np.maximum(den.sum(0), 1e-12)   # 每个分量的相对失配
        val = float(np.mean(ratio))
        return (val, n, ratio) if with_pairs else val

    return J


def _exact_pair_table(P, B, h_candidates, tol=1e-6):
    """纯网格配对交叉校验：只统计镜像点精确落在实测坐标上的配对（无插值误差）。"""
    tree = cKDTree(P)
    rows = []
    for h in h_candidates:
        Pm = P.copy()
        Pm[:, 2] = 2.0 * h - P[:, 2]
        d, idx = tree.query(Pm, k=1, workers=-1)
        m = d <= tol
        if m.sum() < 50:
            continue
        target = B[m] * MIRROR_SIGN
        delta = B[idx[m]] - target
        rows.append((float(h), np.sqrt((delta ** 2).mean(0)), int(m.sum())))
    rows.sort(key=lambda r: r[1].mean())
    return rows


def estimate_symmetry_plane(P, B, snap=True, verbose=True, **kw):
    """从测量数据自动估计 Z 对称面高度。

    Parameters
    ----------
    P, B : 测点 [mm] 与磁场 [uT]
    snap : 是否把结果吸附到"2h 落在实测 z 网格上"的候选（可免掉镜像插值误差）
    其余关键字见模块顶部 ``DEFAULTS``。

    Returns
    -------
    SymmetryResult
    """
    cfg = dict(DEFAULTS)
    cfg.update(kw)
    rng = np.random.default_rng(0)

    cols = _build_columns(P, B, cfg["min_col_pts"], cfg["min_col_span"])
    if verbose:
        print(f"检测器: 有效 z 列 {len(cols)} 条（要求 >= {cfg['min_col_pts']} 点且跨度 >= "
              f"{cfg['min_col_span']} mm）")
    if not cols:
        return SymmetryResult(np.nan, np.nan, None, np.nan, (np.nan, np.nan), 0, 0,
                              False, "没有可用的 z 列，无法检测；请退回几何中面",
                              np.empty((0, 3)), [], float("nan"), 0)

    z_lo = min(c[2][0] for c in cols)
    z_hi = max(c[2][-1] for c in cols)

    # ---- 粗扫（抽样列，抑制耗时）----
    idx = np.arange(len(cols))
    if len(cols) > cfg["coarse_cols"]:
        idx = np.sort(rng.choice(len(cols), cfg["coarse_cols"], replace=False))
    zg_c = np.arange(z_lo, z_hi + 1e-9, cfg["dz_coarse"])
    Vc, Sc = _resample([cols[i] for i in idx], zg_c, cfg["dz_coarse"])
    Jc = _evaluator(zg_c, Vc, Sc)
    hs_c = np.arange(z_lo, z_hi + 1e-9, cfg["step_coarse"])
    vals_c = np.array([Jc(h) for h in hs_c])
    if np.all(np.isnan(vals_c)):
        return SymmetryResult(np.nan, np.nan, None, np.nan, (np.nan, np.nan), 0, 0,
                              False, "粗扫找不到任何有效镜像配对（数据不跨越对称面）；"
                                     "请退回几何中面 (MAGNET_Z_MIN+MAGNET_Z_MAX)/2",
                              np.empty((0, 3)), [], float("nan"), len(cols))
    i_c = int(np.nanargmin(vals_c))
    j_med = float(np.nanmedian(vals_c))
    if not np.isfinite(vals_c[i_c]) or vals_c[i_c] > 0.3:
        return SymmetryResult(np.nan, np.nan, None, np.nan, (np.nan, np.nan), 0, 0,
                              False, f"粗扫最小值 J={vals_c[i_c]:.3f} 过高（≈1 表示无镜像对称性），"
                                     f"疑似单侧数据；请退回几何中面",
                              np.column_stack([hs_c, vals_c, np.zeros_like(hs_c)]), [], j_med,
                              len(cols))

    # ---- 精扫（全部列，只在粗扫凹坑附近）----
    lo = max(z_lo, hs_c[i_c] - 2.0)
    hi = min(z_hi, hs_c[i_c] + 2.0)
    zg_f = np.arange(z_lo, z_hi + 1e-9, cfg["dz_fine"])
    Vf, Sf = _resample(cols, zg_f, cfg["dz_fine"])
    Jf = _evaluator(zg_f, Vf, Sf)
    hs_f = np.arange(lo, hi + 1e-9, cfg["step_fine"])
    res_f = [Jf(h, True) for h in hs_f]
    vals_f = np.array([r[0] for r in res_f])
    pairs_f = np.array([r[1] for r in res_f])
    ratios_f = np.array([r[2] for r in res_f])
    if np.all(np.isnan(vals_f)):
        return SymmetryResult(np.nan, np.nan, None, np.nan, (np.nan, np.nan), 0, 0,
                              False, "精扫无有效配对；请退回几何中面",
                              np.empty((0, 3)), [], j_med, len(cols))
    i_f = int(np.nanargmin(vals_f))
    h_raw = float(hs_f[i_f])
    j_min = float(vals_f[i_f])
    n_pairs = int(pairs_f[i_f])
    comp_ratio = ratios_f[i_f]
    good = vals_f < 1.1 * j_min
    band = (float(hs_f[good].min()), float(hs_f[good].max()))

    # 抛物线细化（仅在凹坑内部，避免被边界带偏）
    if 0 < i_f < len(hs_f) - 1:
        ya, yb, yc = vals_f[i_f - 1], vals_f[i_f], vals_f[i_f + 1]
        den = ya - 2 * yb + yc
        if np.isfinite(den) and den > 0:
            h_raw = float(hs_f[i_f] - 0.5 * (yc - ya) / den * cfg["step_fine"])

    # ---- 吸附到实测 z 网格 ----
    z_planes = np.unique(np.round(P[:, 2], 6))
    cand = np.unique(np.round(0.5 * (z_planes[:, None] + z_planes[None, :]), 6))
    cand = cand[(cand >= z_lo) & (cand <= z_hi)]
    h_snap = None
    if snap and len(cand):
        h_snap = float(cand[np.argmin(np.abs(cand - h_raw))])
    snap_tol = max(0.3, band[1] - band[0])
    h_use = h_raw
    snap_used = False
    if h_snap is not None and abs(h_raw - h_snap) <= snap_tol:
        h_use = h_snap
        snap_used = True

    # ---- 接受判据 ----
    ratio = j_min / j_med if j_med > 0 else np.inf
    reasons = []
    accepted = True
    if n_pairs < cfg["min_pairs"]:
        accepted = False
        reasons.append(f"配对数 {n_pairs} < {cfg['min_pairs']}")
    if len(cols) < cfg["min_cols"]:
        accepted = False
        reasons.append(f"有效 z 列 {len(cols)} < {cfg['min_cols']}（数据太少）")
    if np.isfinite(vals_c[i_c]) and vals_c[i_c] <= 0.3 and (band[1] - band[0]) > cfg["band_max"]:
        accepted = False
        reasons.append(f"凹陷带宽 {band[1] - band[0]:.2f} mm > {cfg['band_max']} mm（不尖锐）")
    if j_min > cfg["j_max"]:
        accepted = False
        reasons.append(f"J={j_min:.4f} > {cfg['j_max']}（镜像不一致）")
    if ratio > cfg["j_ratio_max"]:
        accepted = False
        reasons.append(f"J/J_median={ratio:.3f} > {cfg['j_ratio_max']}（凹陷不够深）")

    # ---- 纯网格配对交叉校验（只在 h_raw 附近 ±3 mm）----
    near = cand[np.abs(cand - h_raw) <= 3.0]
    exact_table = _exact_pair_table(P, B, near) if len(near) else []
    exact_note = ""
    if exact_table:
        exact_note = ("  纯网格配对最优 h = %.2f mm" % exact_table[0][0])

    if accepted:
        reason = (f"接受: J={j_min:.4e} << J_median={j_med:.3f}, 配对数 {n_pairs}, "
                  f"带宽 {band[1] - band[0]:.2f} mm"
                  + (f"；已吸附到测量 z 网格 h={h_use:.3f} mm" if snap_used else "")
                  + exact_note)
    else:
        reason = ("不接受(" + "; ".join(reasons) + ")，回落几何中面；"
                  + "该数据可能不跨越对称面或磁铁本身不对称")

    curve_c = np.column_stack([hs_c, vals_c, np.zeros_like(hs_c)])
    curve_f = np.column_stack([hs_f, vals_f, pairs_f.astype(float)])
    curve = np.vstack([curve_c, curve_f])

    res = SymmetryResult(
        h=float(h_use) if accepted else float("nan"),
        h_raw=h_raw, h_snap=h_snap, j_min=j_min, band=band, n_pairs=n_pairs,
        n_cols=len(cols), accepted=accepted, reason=reason, curve=curve,
        exact_table=exact_table, j_median=j_med, n_used_cols=len(cols),
        comp_ratio=np.asarray(comp_ratio, float),
        comp_snr=_component_snr(P, B),
    )
    if verbose:
        print(res.summary())
    return res


# --------------------------------------------------------------------------------------
# 镜像 + 重合处理
# --------------------------------------------------------------------------------------
def legacy_mirror(P, B, h):
    """旧版（B2M.ipynb）的无条件镜像公式，用于回归对照。

    z -> 2h - z, Bz -> -Bz，不去重、不检查重合。
    """
    Pm = np.asarray(P, float).copy()
    Pm[:, 2] = 2.0 * h - Pm[:, 2]
    Bm = np.asarray(B, float) * MIRROR_SIGN
    return np.vstack([P, Pm]), np.vstack([B, Bm])


def same_point_set(P1, B1, P2, B2, atol=1e-9):
    """两组 (P, B) 是否为同一点集（允许行顺序不同）。返回 (bool, 说明)。"""
    def key(P, B):
        a = np.column_stack([P, B])
        return a[np.lexsort(a.T[::-1])]

    A = key(np.asarray(P1, float), np.asarray(B1, float))
    C = key(np.asarray(P2, float), np.asarray(B2, float))
    if A.shape != C.shape:
        return False, f"点数/列数不同: {A.shape} vs {C.shape}"
    if not np.allclose(A, C, atol=atol, rtol=0):
        bad = int((~np.isclose(A, C, atol=atol, rtol=0)).any(1).sum())
        return False, f"有 {bad} 行不一致（最大差 {np.abs(A - C).max():.3e}）"
    return True, f"完全一致（{len(A)} 行，仅行顺序可能不同）"


def regression_check(P, B, h, policy="keep_measured", tol=None, verbose=True):
    """回归校验：保证"旧数据集跑出来与今天一致"。

    - 无重合时：新结果必须与旧的无条件镜像公式**逐元素完全一致**（仅行顺序可不同）；
    - 有重合时：新结果必须精确等于"实测点 + 未落入实测 tol 邻域的镜像点"。
    """
    Pn, Bn, rep = apply_mirror(P, B, h, policy=policy, tol=tol, fill=True, verbose=False)
    Pm = np.asarray(P, float).copy()
    Pm[:, 2] = 2.0 * h - Pm[:, 2]
    d, _ = cKDTree(P).query(Pm, k=1, workers=-1)
    keep = d > rep.tol

    if rep.n_overlap == 0:
        Po, Bo = legacy_mirror(P, B, h)
        ok, info = same_point_set(Pn, Bn, Po, Bo)
        info = "无重合 -> 应与旧公式完全一致; " + info
    elif policy == "keep_measured":
        Pe = np.vstack([P, Pm[keep]])
        Be = np.vstack([B, (np.asarray(B, float) * MIRROR_SIGN)[keep]])
        ok, info = same_point_set(Pn, Bn, Pe, Be)
        info = f"重合 {rep.n_overlap} 个镜像副本已剔除; " + info
    else:
        # average：数值由对称化决定，这里只校验点集（数值对称性由 selftest 校验）
        Pe = np.vstack([P, Pm[keep]])
        ok, info = same_point_set(Pn, np.zeros_like(Bn), Pe, np.zeros_like(Pe))
        info = f"重合 {rep.n_overlap} 个镜像副本已对称化; 点集: " + info

    if verbose:
        print(f"[回归校验] {'PASS' if ok else 'FAIL'}: {info}")
    return ok


def default_tol(P) -> float:
    """默认去重容差 = 0.5 x 最近邻距离中位数（自适应各种扫描网格间距）。"""
    tree = cKDTree(P)
    if len(P) > 20000:
        idx = np.random.default_rng(0).choice(len(P), 20000, replace=False)
        d, _ = tree.query(P[idx], k=2, workers=-1)
    else:
        d, _ = tree.query(P, k=2, workers=-1)
    return float(0.5 * np.median(d[:, 1]))


def apply_mirror(P, B, h, policy="keep_measured", tol=None, fill=True, meta=None,
                 n_duplicate_measured=0, verbose=True):
    """按对称面 h 做镜像增强，并处理与实测点重合的部分。

    Parameters
    ----------
    policy : "keep_measured"（默认）镜像点与实测点重合时丢弃镜像副本，实测值原样保留；
             "average" 重合对做对称化平均 B <- (B(p) + R*B(sigma(p)))/2。
    tol    : 去重容差 [mm]，None 时用 ``default_tol``
    fill   : 是否把落在实测未覆盖区域的镜像点补进数据集

    Returns
    -------
    P2, B2, MirrorReport
    """
    P = np.asarray(P, float)
    B = np.asarray(B, float)
    if policy not in ("keep_measured", "average"):
        raise ValueError(f"未知 policy: {policy}")
    if tol is None:
        tol = default_tol(P)

    Pm = P.copy()
    Pm[:, 2] = 2.0 * h - P[:, 2]
    Bm = B * MIRROR_SIGN

    tree = cKDTree(P)
    d, idx = tree.query(Pm, k=1, workers=-1)
    overlap = d <= tol

    # ---- 重合处不一致量诊断 ----
    asym_rms = asym_mean = asym_rel = b_rms = np.zeros(3)
    worst = []
    per_source = {}
    if overlap.any():
        target = B[idx[overlap]]
        delta = Bm[overlap] - target
        asym_rms = np.sqrt((delta ** 2).mean(0))
        asym_mean = delta.mean(0)
        b_rms = np.sqrt((target ** 2).mean(0))
        asym_rel = asym_rms / np.maximum(b_rms, 1e-12)
        order = np.argsort(-np.abs(delta).max(1))[:10]
        sel = np.where(overlap)[0][order]
        worst = [(float(P[i, 0]), float(P[i, 1]), float(P[i, 2]),
                  float(delta[order, 0][j]), float(delta[order, 1][j]),
                  float(delta[order, 2][j])) for j, i in enumerate(sel)]
        if meta and meta.get("files"):
            src = _source_of(P, meta)
            for name in np.unique(src[overlap]) if src is not None else []:
                m = overlap & (src == name)
                if m.sum() == 0:
                    continue
                dd = Bm[m] - B[idx[m]]
                per_source[name] = (int(m.sum()), np.sqrt((dd ** 2).mean(0)))

    # ---- 实测值（可选对称化）----
    Bk = B.copy()
    if policy == "average" and overlap.any():
        acc = B.copy()
        cnt = np.ones(len(P))
        np.add.at(acc, idx[overlap], Bm[overlap])
        np.add.at(cnt, idx[overlap], 1.0)
        Bk = acc / cnt[:, None]

    # ---- 组装 ----
    keep_fill = (~overlap) if fill else np.zeros(len(P), bool)
    if keep_fill.any():
        P2 = np.vstack([P, Pm[keep_fill]])
        B2 = np.vstack([Bk, Bm[keep_fill]])
    else:
        P2, B2 = P.copy(), Bk.copy()

    rep = MirrorReport(
        n_measured=len(P), n_mirrored_total=len(P), n_overlap=int(overlap.sum()),
        n_fill=int(keep_fill.sum()), n_final=len(P2), tol=float(tol), policy=policy,
        fill=bool(fill), h=float(h), asym_rms=asym_rms, asym_mean=asym_mean,
        asym_rel=asym_rel, b_rms=b_rms, worst=worst, per_source=per_source,
        n_duplicate_measured=int(n_duplicate_measured),
        noise_rms=estimate_noise(P, B),
    )
    if verbose:
        print(rep.summary())
    return P2, B2, rep


def _source_of(P, meta):
    """给每个测点标出它来自哪个 csv（按文件顺序拼接的顺序推断）。"""
    files = meta.get("files") or []
    if not files or meta.get("dedup") == "merge":
        return None
    sizes = meta.get("file_sizes")
    if not sizes:
        return None
    src = np.empty(len(P), dtype=object)
    s = 0
    for f, n in zip(files, sizes):
        src[s:s + n] = os.path.basename(f)
        s += n
    if s != len(P):
        return None
    return src


# --------------------------------------------------------------------------------------
# 诊断图
# --------------------------------------------------------------------------------------
def plot_symmetry_diagnostics(res: SymmetryResult, P, B, out_png, fallback_h=None,
                              title_suffix=""):
    """J(h) 曲线 + 典型 z 列镜像对照 + 纯网格配对校验。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
    plt.rcParams["axes.unicode_minus"] = False

    fig = plt.figure(figsize=(14.5, 9.5))
    gs = fig.add_gridspec(2, 3, hspace=0.32, wspace=0.28)

    # (a) J(h)
    ax = fig.add_subplot(gs[0, :2])
    cur = res.curve
    ax.semilogy(cur[:, 0], np.maximum(cur[:, 1], 1e-8), ".", ms=3, color="tab:blue")
    if np.isfinite(res.h_raw):
        ax.axvline(res.h_raw, color="tab:red", lw=1.5,
                   label=f"检测 h_raw={res.h_raw:.2f} mm (J={res.j_min:.3e})")
    if res.h_snap is not None:
        ax.axvline(res.h_snap, color="tab:green", ls="--", lw=1.2,
                   label=f"网格吸附 h_snap={res.h_snap:.2f} mm")
    if np.isfinite(res.h_raw):
        ax.axvspan(res.band[0], res.band[1], color="tab:red", alpha=0.15,
                   label=f"不确定带宽 {res.band_width:.2f} mm")
    ax.axhline(DEFAULTS["j_max"], color="k", ls=":", lw=1.0,
               label=f"接受阈值 J={DEFAULTS['j_max']}")
    if fallback_h is not None and np.isfinite(fallback_h):
        ax.axvline(fallback_h, color="gray", ls="-.", lw=1.2,
                   label=f"几何中面回落 {fallback_h:.2f} mm")
    ax.set_xlabel("候选对称面高度 h / mm")
    ax.set_ylabel("镜像一致性失配 J(h)")
    ax.set_title("对称面自动检测" + title_suffix)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="best")

    # (b) 纯网格配对校验
    ax = fig.add_subplot(gs[0, 2])
    if res.exact_table:
        hh = [r[0] for r in res.exact_table]
        rr = [r[1].mean() for r in res.exact_table]
        ax.plot(hh, rr, "o-", color="tab:purple", ms=5)
        best = min(res.exact_table, key=lambda r: r[1].mean())[0]
        ax.axvline(best, color="tab:purple", ls="--", lw=1.0, label=f"最优 {best:.2f} mm")
        ax.set_xlabel("h / mm")
        ax.set_ylabel("rms |ΔB| / μT")
        ax.set_title("纯网格配对（无插值误差）")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "无精确网格配对", ha="center", va="center")
        ax.set_axis_off()

    # (c) 三条典型 z 列的镜像对照
    cand_cols = []
    for (x, y), g in pd.DataFrame(P, columns=["x", "y", "z"]).assign(
            Bx=B[:, 0], By=B[:, 1], Bz=B[:, 2]).groupby(["x", "y"]):
        if len(g) >= 10 and (g["z"].max() - g["z"].min()) >= 10:
            cand_cols.append((x, y, len(g)))
    cand_cols.sort(key=lambda t: -t[2])
    picks = cand_cols[:3]

    dfP = pd.DataFrame(P, columns=["x", "y", "z"])
    for j, (x, y, _n) in enumerate(picks):
        if j >= 3:
            break
        ax = fig.add_subplot(gs[1, j])
        g = dfP[(dfP.x == x) & (dfP.y == y)].sort_values("z")
        i = g.index.to_numpy()
        zz = P[i, 2]
        for k, (name, c) in enumerate(zip("Bx By Bz".split(), ["tab:blue", "tab:orange", "tab:green"])):
            ax.plot(zz, B[i, k], "o-", ms=3, color=c, label=name)
            if np.isfinite(res.h):
                zm = 2 * res.h - zz
                mm = (zm >= zz.min()) & (zm <= zz.max())
                ax.plot(zm[mm], MIRROR_SIGN[k] * B[i, k][mm], "x", ms=4, color=c, alpha=0.7)
        if np.isfinite(res.h):
            ax.axvline(res.h, color="tab:red", lw=1.2)
        ax.set_xlabel("z / mm")
        ax.set_ylabel("B / μT")
        ax.set_title(f"列 (x={x:g}, y={y:g})\n实线=实测，×=镜像预测", fontsize=9)
        ax.grid(alpha=0.3)
        if j == 0:
            ax.legend(fontsize=8)

    fig.suptitle("磁铁 Z 对称面检测诊断", fontsize=13)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    import matplotlib.pyplot as plt2
    plt2.close(fig)
    print(f"已保存 {out_png}")
    return out_png


def plot_mirror_overlap(rep: MirrorReport, P, B, out_png, repeatability=None):
    """重合点统计 + 不一致量空间分布。"""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
    plt.rcParams["axes.unicode_minus"] = False

    Pm = P.copy()
    Pm[:, 2] = 2.0 * rep.h - P[:, 2]
    Bm = B * MIRROR_SIGN
    d, idx = cKDTree(P).query(Pm, k=1, workers=-1)

    fig = plt.figure(figsize=(14.5, 9.0))
    gs = fig.add_gridspec(2, 3, hspace=0.34, wspace=0.28)

    # (a) 每个 z 平面的重合比例
    ax = fig.add_subplot(gs[0, 0])
    zp = np.unique(np.round(P[:, 2], 6))
    frac = [float((d[np.isclose(P[:, 2], z, atol=1e-6)] <= rep.tol).mean()) * 100 for z in zp]
    ax.bar(zp, frac, width=1.6, color="tab:blue")
    ax.axvline(rep.h, color="tab:red", lw=1.2, label=f"h={rep.h:.2f} mm")
    ax.set_xlabel("z / mm")
    ax.set_ylabel("镜像重合比例 / %")
    ax.set_title("重合点按 z 分布")
    ax.grid(alpha=0.3, axis="y")
    ax.legend(fontsize=8)

    # (b) 重合数 vs tol
    ax = fig.add_subplot(gs[0, 1])
    tols = np.linspace(0.02, 2.5, 60)
    cnt = [(d <= t).sum() for t in tols]
    ax.plot(tols, cnt, "-", color="tab:purple")
    ax.axvline(rep.tol, color="tab:red", ls="--", lw=1.2, label=f"采用 tol={rep.tol:.2f} mm")
    ax.set_xlabel("去重容差 tol / mm")
    ax.set_ylabel("被丢弃的镜像点数")
    ax.set_title("重合数对 tol 的敏感性")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)

    # (c) 不一致量 vs |B|
    ax = fig.add_subplot(gs[0, 2])
    m = d <= rep.tol
    if m.any():
        delta = Bm[m] - B[idx[m]]
        mag = np.linalg.norm(B[idx[m]], axis=1) + 1e-9
        dm = np.linalg.norm(delta, axis=1)
        ax.loglog(mag, np.maximum(dm, 1e-3), ".", ms=2, alpha=0.25, color="tab:orange")
        bins = np.logspace(np.log10(mag.min()), np.log10(mag.max()), 24)
        bi = np.digitize(mag, bins)
        bx = [mag[bi == k].mean() for k in range(1, len(bins)) if (bi == k).sum() > 5]
        by = [np.sqrt((dm[bi == k] ** 2).mean()) for k in range(1, len(bins)) if (bi == k).sum() > 5]
        ax.loglog(bx, by, "o-", color="k", ms=4, label="分箱 rms")
        if repeatability:
            ax.axhline(repeatability, color="tab:green", ls="--", lw=1.2,
                       label=f"扫描重复性 ≈{repeatability:.0f} μT")
        ax.set_xlabel("|B| / μT")
        ax.set_ylabel("|ΔB| / μT")
        ax.set_title("重合处不一致量")
        ax.grid(alpha=0.3, which="both")
        ax.legend(fontsize=8)
    else:
        ax.text(0.5, 0.5, "无重合点", ha="center", va="center")
        ax.set_axis_off()

    # (d) 最差不一致点的空间位置（x-y，按 z 着色）
    ax = fig.add_subplot(gs[1, :2])
    if m.any():
        delta = Bm[m] - B[idx[m]]
        dm = np.linalg.norm(delta, axis=1)
        thr = np.quantile(dm, 0.99)
        sel = dm >= thr
        sc = ax.scatter(P[idx[m]][sel, 0], P[idx[m]][sel, 1], c=P[idx[m]][sel, 2],
                        s=8, cmap="turbo")
        fig.colorbar(sc, ax=ax, label="z / mm")
        ax.set_xlabel("x / mm")
        ax.set_ylabel("y / mm")
        ax.set_aspect("equal")
        ax.set_title(f"不一致量最大的 1% 重合点（|ΔB| ≥ {thr:.0f} μT）")
        ax.grid(alpha=0.3)

    # (e) 各分量直方图
    ax = fig.add_subplot(gs[1, 2])
    if m.any():
        delta = Bm[m] - B[idx[m]]
        for k, (name, c) in enumerate(zip("Bx By Bz".split(), ["tab:blue", "tab:orange", "tab:green"])):
            ax.hist(delta[:, k], bins=80, histtype="step", color=c, label=name)
        ax.axvline(0, color="k", lw=0.8)
        ax.set_yscale("log")
        ax.set_xlabel("ΔB = 镜像 − 实测 / μT")
        ax.set_ylabel("计数")
        ax.set_title("重合处残差分布")
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)

    fig.suptitle(f"镜像重合诊断（h={rep.h:.3f} mm, {rep.policy}, tol={rep.tol:.2f} mm）", fontsize=13)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    print(f"已保存 {out_png}")
    return out_png


# --------------------------------------------------------------------------------------
# 自测
# --------------------------------------------------------------------------------------
def _synthetic_dataset(h0=0.0, noise=0.0, seed=0, one_sided=False, shift=0.0):
    """在 z=h0 平面内放几个面内磁偶极子，生成严格关于 z=h0 镜像对称的场。

    偶极子刻意放在网格点之外，避免奇点；量级调到与真实数据同量级（~1e4 uT）。
    """
    rng = np.random.default_rng(seed)
    # (位置 x, y, 面内方向, 偶极矩大小)
    dip = [(-19.0, 1.0, (1.0, 0.2), 1.4e7),
           (19.0, 1.0, (1.0, 0.0), 1.4e7),
           (1.0, 19.0, (0.0, 1.0), -1.0e7)]
    xs = np.arange(-30.0, 30.1, 2.0)
    ys = np.arange(-30.0, 30.1, 2.0)
    zs = np.arange(-12.0, 12.1, 1.0)
    X, Y, Zg = np.meshgrid(xs, ys, zs, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Zg.ravel()])
    B = np.zeros_like(pts)
    for dx, dy, (mx, my), moment in dip:
        m = np.array([mx, my, 0.0])
        m = m / np.hypot(m[0], m[1]) * moment
        u = pts - np.array([dx, dy, 0.0])
        r = np.maximum(np.linalg.norm(u, axis=1), 1e-3)
        B += 1e-7 * (3 * (u @ m)[:, None] * u / r[:, None] ** 5
                     - m / r[:, None] ** 3) * 1e6
    if noise:
        B += rng.normal(0, noise, B.shape)
    P = pts.copy()
    P[:, 2] += h0 + shift
    if one_sided:
        m = P[:, 2] >= (h0 + shift)
        P, B = P[m], B[m]
    return P, B


def selftest(verbose=True):
    """合成数据自测：已知对称面、平移、单侧拒绝、噪声鲁棒性。"""
    ok = True

    def check(name, cond, info):
        nonlocal ok
        ok = ok and cond
        print(f"  [{'PASS' if cond else 'FAIL'}] {name}: {info}")

    print("=== symmetry.py 自测 ===")
    P, B = _synthetic_dataset(h0=0.0, noise=5.0, seed=1)
    r = estimate_symmetry_plane(P, B, verbose=False)
    check("人造对称面 h0=0", r.accepted and abs(r.h - 0.0) <= 0.2,
          f"h={r.h:.3f} (h_raw={r.h_raw:.3f}) J={r.j_min:.3e} accepted={r.accepted}")

    P, B = _synthetic_dataset(h0=0.0, shift=3.7, noise=5.0, seed=2)
    r = estimate_symmetry_plane(P, B, verbose=False)
    check("z 刚性平移 +3.7", r.accepted and abs(r.h - 3.7) <= 0.25,
          f"h={r.h:.3f} J={r.j_min:.3e}")

    P, B = _synthetic_dataset(h0=0.0, one_sided=True, noise=5.0, seed=3)
    r = estimate_symmetry_plane(P, B, verbose=False)
    check("单侧数据必须拒绝", not r.accepted, f"accepted={r.accepted} reason={r.reason[:60]}")

    P, B = _synthetic_dataset(h0=0.0, noise=200.0, seed=4)
    r = estimate_symmetry_plane(P, B, verbose=False)
    check("200 μT 噪声", r.accepted and abs(r.h) <= 0.25, f"h={r.h:.3f} J={r.j_min:.3e}")

    P, B = _synthetic_dataset(h0=-2.0, noise=20.0, seed=5)
    r = estimate_symmetry_plane(P, B, verbose=False)
    check("人造对称面 h0=-2", r.accepted and abs(r.h + 2.0) <= 0.2, f"h={r.h:.3f}")

    # 镜像/去重逻辑
    P, B = _synthetic_dataset(h0=7.0, noise=10.0, seed=6)
    P2, B2, rep = apply_mirror(P, B, 7.0, policy="keep_measured", verbose=False)
    d, _ = cKDTree(P2).query(P2, k=2, workers=-1)
    check("keep_measured 无重合", rep.n_overlap > 0 and (d[:, 1] > rep.tol).all(),
          f"重合={rep.n_overlap} 补点={rep.n_fill} 最终={rep.n_final} min间距={d[:,1].min():.3f}")
    check("重合处不一致量已统计", float(np.abs(rep.asym_rms).max()) > 0,
          f"asym_rms={np.round(rep.asym_rms, 2)} μT")

    P3, B3, rep3 = apply_mirror(P, B, 7.0, policy="average", verbose=False)
    # 对称化后：每个测点镜像位置上的值应与其镜像预测完全一致（差 ~1e-12）
    tree3 = cKDTree(P3)
    Pm = P3.copy()
    Pm[:, 2] = 2.0 * 7.0 - P3[:, 2]
    d3, i3 = tree3.query(Pm, k=1, workers=-1)
    ok3 = d3 <= rep3.tol
    resid = np.abs(B3[i3][ok3] - (B3 * MIRROR_SIGN)[ok3]).max()
    check("average 对称化精确成立", rep3.n_overlap > 0 and ok3.all() and resid < 1e-6,
          f"重合={rep3.n_overlap} 有配对={ok3.sum()}/{len(P3)} 最大镜像失配={resid:.3e} μT")

    print(f"=== 自测{'全部通过' if ok else '存在失败项'} ===")
    return ok


# --------------------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------------------
def main(argv=None):
    ap = argparse.ArgumentParser(description="磁铁 Z 对称面自动检测 + 镜像重合处理")
    ap.add_argument("folder", nargs="?", default=None,
                    help="数据文件夹（读取其中的 *.csv），例如 data/raw/UshapeNormal_New_prob")
    ap.add_argument("--pattern", default=None, help="直接给 glob 模式，覆盖 folder")
    _REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--out", default=os.path.join(_REPO, "figures"),
                    help="诊断图输出目录（默认仓库根的 figures/）")
    ap.add_argument("--policy", default="keep_measured", choices=["keep_measured", "average"],
                    help="重合点处理策略")
    ap.add_argument("--tol", type=float, default=None, help="去重容差 [mm]，默认自适应")
    ap.add_argument("--no-fill", action="store_true", help="不补镜像点到未测区域")
    ap.add_argument("--manual-h", type=float, default=None, help="跳过检测，直接指定 h [mm]")
    ap.add_argument("--fallback-h", type=float, default=None,
                    help="检测失败时的回落值（几何中面），仅用于在图上标注")
    ap.add_argument("--selftest", action="store_true", help="运行合成数据自测")
    ap.add_argument("--regression", action="store_true",
                    help="与旧版无条件镜像公式做回归对照")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args(argv)

    if args.selftest:
        return 0 if selftest() else 1

    pattern = args.pattern or (os.path.join(args.folder, "*.csv") if args.folder else None)
    if pattern is None:
        ap.error("请给出 folder 或 --pattern")

    P, B, meta = load_measurements(pattern, verbose=not args.quiet)
    if args.manual_h is not None:
        h = float(args.manual_h)
        res = SymmetryResult(h, h, None, float("nan"), (h, h), 0, 0, True,
                             "手动指定", np.empty((0, 3)))
        print(f"手动指定对称面 h = {h} mm")
    else:
        res = estimate_symmetry_plane(P, B, verbose=not args.quiet)
        if not res.accepted:
            fb = args.fallback_h
            if fb is None:
                raise SystemExit("检测未通过且没有给出 --fallback-h，请手工确认对称面高度。")
            h = float(fb)
            print(f"检测未通过，回落几何中面 h = {h} mm")
        else:
            h = res.h

    P2, B2, rep = apply_mirror(P, B, h, policy=args.policy, tol=args.tol, fill=not args.no_fill,
                               meta=meta, n_duplicate_measured=meta.get("n_dup", 0),
                               verbose=not args.quiet)

    os.makedirs(args.out, exist_ok=True)
    tag = os.path.basename(args.folder.rstrip("/\\")) if args.folder else "data"
    plot_symmetry_diagnostics(res, P, B, os.path.join(args.out, "symmetry_plane_detection.png"),
                              fallback_h=args.fallback_h, title_suffix=f"  [{tag}]")
    plot_mirror_overlap(rep, P, B, os.path.join(args.out, "mirror_overlap.png"))

    print("\n可直接粘贴进 notebook 的参数：")
    print(f"  SYMMETRY_H = {h:.4f}   # mm "
          f"({'手动指定' if args.manual_h is not None else 'symmetry.py 自动检测'})")
    print(f"  OVERLAP_POLICY = '{args.policy}'   # 重合 {rep.n_overlap} 个镜像副本"
          f"{'已丢弃' if args.policy == 'keep_measured' else '已对称化平均'}")
    print(f"  数据装配: 实测 {rep.n_measured} -> 最终 {rep.n_final} 点")

    if args.regression:
        ok = regression_check(P, B, h, policy=args.policy, tol=args.tol)
        return 0 if ok else 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
