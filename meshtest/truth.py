# -*- coding: utf-8 -*-
"""合成真值与合成磁场：先造一个已知的磁化强度分布，正演出磁场，再交给反演。

为什么要独立细网格
------------------
真值场用一个 **0.5 mm 的独立细网格**（而不是反演网格）正演出，避免 inverse crime ——
否则"用同一套离散化正演再反演"会得到虚高的还原度，掩盖正则化本身的偏差。

两套真值（按需求：均匀块 + 带接缝块）
------------------------------------
- ``uniform``：整块均匀磁化 M = (M0, 0, 0)。最硬的判据：任何按网格档位的偏置都会暴露。
- ``seam``   ：沿 **y** 分成两块，|M| 不同（1.25e5 / 0.7e5 A/m，同向），对比度 1.79:1。
              接缝故意与网格的细化分界（沿 x）**正交**：这样每一档网格内部都同时含有
              高/低磁化两半，"按档偏置"就不可能与真值的空间结构混淆 —— 否则测到的
              只是"接缝恰好落在细档"这种假信号（已经踩过这个坑）。
"""

from __future__ import annotations

import os
import numpy as np

from b2m_core import matvec

__all__ = ["make_truth", "make_measurements", "synthesize", "TRUTH_KINDS"]

TRUTH_KINDS = ("uniform", "seam")

M0_DEFAULT = 1.0e5          # A/m
SEAM_RATIO = (1.0e5, -1.0e5)   # 接缝两侧的 Mx（反向：真实拼接磁铁的接缝）

CACHE_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


def _fine_points(magnet_box, h=1.0):
    """磁体内部 h mm 均匀网格（真值离散化，与反演网格无关；默认 1mm，比反演网格细 4 倍以上）。"""
    box = np.asarray(magnet_box, float)
    axes = [np.arange(box[0][k] + h / 2, box[1][k], h) for k in range(3)]
    X, Y, Z = np.meshgrid(*axes, indexing="ij")
    pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    vol = np.full(len(pts), (h / 1000.0) ** 3)      # m³
    return pts, vol


def make_truth(kind, magnet_box, m0=M0_DEFAULT, seam_frac=0.5, h=1.0):
    """返回 (fine_pts_mm, fine_vol_m3, M_fine, meta)。

    M_fine : (n, 3) A/m，定义在 0.5 mm 细网格的单元中心（用于正演）
    """
    if kind not in TRUTH_KINDS:
        raise ValueError(f"未知真值类型 {kind}（可选 {TRUTH_KINDS}）")
    pts, vol = _fine_points(magnet_box, h)
    M = np.zeros_like(pts)
    box = np.asarray(magnet_box, float)
    if kind == "uniform":
        M[:, 0] = m0
        meta = dict(kind=kind, m_left=m0, m_right=m0, seam_x=None)
    else:
        y_seam = box[0][1] + seam_frac * (box[1][1] - box[0][1])
        low = pts[:, 1] < y_seam
        M[low, 0] = SEAM_RATIO[0]
        M[~low, 0] = SEAM_RATIO[1]
        meta = dict(kind=kind, m_left=SEAM_RATIO[0], m_right=SEAM_RATIO[1], seam_x=y_seam,
                    seam_axis=1)
    return pts, vol, M, meta


def make_measurements(domain, mode="dense", z_off=24.0, step=4.0, mag_box=None):
    """测点：磁体外部的扫描点（模拟真实探头只能在磁体外测）。

    mode : "dense"  顶面 + 四个侧面（良态）
           "sparse" 顶面一层、点距粗（中等欠定）
           "oneside" 只有远处单一平面（高度欠定：磁化的 z 分布几乎不可解，
                     正则化在很宽的 λ 区间内都在"塑造"解，指纹最容易显形）
    """
    Lx, Ly, Lz = domain
    if mode == "dense":
        gx = np.arange(0.0, Lx + 1e-9, step)
        gy = np.arange(0.0, Ly + 1e-9, step)
        gz = np.arange(0.0, Lz + 1e-9, step)
        X, Y = np.meshgrid(gx, gy, indexing="ij")
        top = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, Lz + z_off)])
        X, Z = np.meshgrid(gx, gz, indexing="ij")
        side_y0 = np.column_stack([X.ravel(), np.full(X.size, -z_off), Z.ravel()])
        side_y1 = np.column_stack([X.ravel(), np.full(X.size, Ly + z_off), Z.ravel()])
        Y, Z = np.meshgrid(gy, gz, indexing="ij")
        side_x0 = np.column_stack([np.full(Y.size, -z_off), Y.ravel(), Z.ravel()])
        side_x1 = np.column_stack([np.full(Y.size, Lx + z_off), Y.ravel(), Z.ravel()])
        pts = np.vstack([top, side_y0, side_y1, side_x0, side_x1])
    elif mode == "oneside":
        gx = np.arange(0.0, Lx + 1e-9, 2 * step)
        gy = np.arange(0.0, Ly + 1e-9, 2 * step)
        X, Y = np.meshgrid(gx, gy, indexing="ij")
        pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, Lz + z_off)])
    elif mode == "sparse":
        gx = np.arange(0.0, Lx + 1e-9, 2 * step)
        gy = np.arange(0.0, Ly + 1e-9, 2 * step)
        X, Y = np.meshgrid(gx, gy, indexing="ij")
        pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, Lz + z_off)])
    else:
        raise ValueError(f"未知 mode {mode}")
    return pts.astype(np.float64)


def synthesize(kind, domain, magnet_box, meas_mode="dense", step=4.0,
               noise_frac=0.01, seed=0, m0=M0_DEFAULT, h=1.0, z_off=24.0, verbose=True):
    """生成一套完整的合成算例。

    Returns
    -------
    dict(B_meas, B_true, meas, M_fine, fine_pts, fine_vol, sigma, truth_meta, truth_fn)
      truth_fn(points_mm) -> M_true (n,3)：在任意点上取真值（用于与反演单元中心比对）
    """
    fine_pts, fine_vol, M_fine, tmeta = make_truth(kind, magnet_box, m0=m0, h=h)
    meas = make_measurements(domain, mode=meas_mode, z_off=z_off, step=step)
    B_true = matvec(np.ascontiguousarray(M_fine.ravel()), meas / 1000.0,
                    fine_pts / 1000.0, fine_vol)
    rng = np.random.default_rng(seed)
    sigma = noise_frac * float(np.sqrt((B_true ** 2).mean()))
    B_meas = B_true + rng.normal(0.0, sigma, B_true.shape)

    box = np.asarray(magnet_box, float)

    def truth_fn(points_mm):
        p = np.asarray(points_mm, float)
        M = np.zeros_like(p)
        if kind == "uniform":
            M[:, 0] = m0
        else:
            ax = tmeta.get("seam_axis", 1)
            M[:, 0] = np.where(p[:, ax] < tmeta["seam_x"], tmeta["m_left"], tmeta["m_right"])
        return M

    out = dict(B_meas=B_meas, B_true=B_true, meas=meas, M_fine=M_fine,
               fine_pts=fine_pts, fine_vol=fine_vol, sigma=sigma,
               truth_meta=tmeta, truth_fn=truth_fn, magnet_box=box,
               noise_frac=noise_frac, meas_mode=meas_mode, seed=seed)
    if verbose:
        print(f"[truth] {kind} / {meas_mode}: 真值细网格 {len(fine_pts)} 点, "
              f"测点 {len(meas)} 个, |B|rms={np.sqrt((B_true**2).mean()):.3e} T, "
              f"注入噪声 σ={sigma:.3e} T ({noise_frac*100:.1f}%)")
    return out


def cache_path(name):
    os.makedirs(CACHE_DIR, exist_ok=True)
    return os.path.join(CACHE_DIR, name)
