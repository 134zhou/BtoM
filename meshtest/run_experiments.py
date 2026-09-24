# -*- coding: utf-8 -*-
"""meshtest 实验驱动：E0–E6。

    python run_experiments.py --quick          # 粗网格、少量 λ，几分钟
    python run_experiments.py --full           # 更细网格、完整 λ 扫描（默认）
    python run_experiments.py --only E2        # 只跑某个实验

结论写入 results.md，图写入 figures/。

为什么所有 λ 都写成 λ = λ_bal·10^k
-----------------------------------
λ 的绝对大小同时依赖网格、长度单位、权重约定和数据标度，直接扫 "1e-10, 1e-9, ..." 会
整段落在过正则化区（什么都看不出来，已踩过）。这里统一用**力的平衡**标定：

    λ_bal = ‖∇E_data(M_ref)‖ / ‖∇E_reg(M_ref)‖|_(λ=1)

即"数据拉力"与"Huber/TV 摩擦"在该参考模型处相等的强度。于是 k<0 欠正则化、
k>0 过正则化、k≈0 两项同量级；跨变体、跨网格、跨单位都可比。
生产代码的 λ=1e-10 会换算成"相当于本变体标度下的 λ，等于 λ_bal×10^几"标注出来。
"""

from __future__ import annotations

import argparse
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

# 以脚本方式运行（输出重定向到文件）时 Windows 控制台默认 GBK，会把中文/上标打印炸掉
try:
    sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    sys.stderr.reconfigure(encoding="utf-8", errors="replace")
except Exception:
    pass

from scipy.sparse.linalg import LinearOperator, lsqr

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei"]
plt.rcParams["axes.unicode_minus"] = False

from b2m_core import (REG_VARIANTS, build_mesh, build_regularizer, calibrate_lambda, chi,
                      check_gradient, data_term, magnet_cells, magnet_levels, make_forward,
                      make_reference_model, metrics, objective_and_gradient, operator_coef,
                      production_objective, reg_term, solve, to_m, variable_scale,
                      volume_weighted_norm)
import truth as T

FIG_DIR = os.path.join(HERE, "figures")
RES_PATH = os.path.join(HERE, "results.md")
LAMBDA_PROD = 1e-10        # 生产代码的 lambda_reg（volume_raw / mm³ 约定下）
EPS = 4e3                  # huber_epsilon

RESULTS = []


def log(exp, what, ok, detail=""):
    RESULTS.append(dict(exp=exp, what=what, ok=bool(ok), detail=detail))
    print(f"  [{'PASS' if ok else 'FAIL'}] {exp} {what}  {detail}")


# ======================================================================================
def default_geom():
    """域 64×64×16，磁体块 32×32×16 居中。

    关键：盒边界 (16,48,0,16) 必须落在**所有**细化档的格边界上（基础网格 4mm/2mm、
    粗档 8mm/4mm 都能整除 16），否则 coarse 单元会跨出盒子，按"单元中心在盒内"选出的
    磁体体积会大于盒体积 —— 这个坑踩过一次：体积偏大 39%，导致幅值整体系统性偏小、
    还原误差 45%~150%。setup() 里有断言，别再犯。
    """
    return (64.0, 64.0, 16.0), ((16.0, 16.0, 0.0), (48.0, 48.0, 16.0))


def setup(domain, magnet_box, style, nbase, variant, meas_mode, kind,
          step=4.0, noise_frac=0.01, seed=0, mag_levels=1, split_axis=0, z_off=24.0,
          verbose=False):
    mesh = build_mesh(domain, nbase, style=style, magnet_box=magnet_box,
                      split_axis=split_axis, split_frac=0.5, mag_levels=mag_levels)
    midx = magnet_cells(mesh, magnet_box)
    cc = mesh.cell_centers[midx]
    vol = mesh.cell_volumes[midx]
    is_mag = np.zeros(mesh.nC, bool)
    is_mag[midx] = True
    reg = build_regularizer(mesh, is_mag, variant=variant)

    data = T.synthesize(kind, domain, magnet_box, meas_mode=meas_mode, step=step,
                        noise_frac=noise_frac, seed=seed, z_off=z_off, verbose=False)
    mv, rv = make_forward(data["meas"] / 1000.0, cc / 1000.0, vol / 1e9)
    st = dict(mesh=mesh, midx=midx, cc=cc, vol=vol, reg=reg, data=data, mv=mv, rv=rv,
              b=data["B_meas"].ravel(), M_true=data["truth_fn"](cc),
              level=magnet_levels(vol), variant=variant,
              tag=f"{style}|{variant}|{meas_mode}|{kind}")
    v_box = float(np.prod(np.array(magnet_box[1]) - np.array(magnet_box[0])))
    v_sel = float(vol.sum())
    if abs(v_sel - v_box) / v_box > 1e-6:
        raise RuntimeError(
            f"几何未对齐：磁体单元体积和 {v_sel:.6g} ≠ 盒体积 {v_box:.6g}"
            f"（相对差 {(v_sel-v_box)/v_box*100:+.2f}%）。请让磁体盒边界落在所有细化档的格边界上。")
    st["M_ref"] = make_reference_model(cc, M_true=st["M_true"])
    st["lam_bal"] = calibrate_lambda(st["M_ref"].ravel(), st["b"], st["mv"], st["rv"],
                                     reg["G_sub"], reg["w"], EPS)
    if verbose:
        print(f"  [{st['tag']}] nC={mesh.nC} 磁体单元={len(midx)} "
              f"V档={np.unique(np.round(vol, 3))} λ_bal={st['lam_bal']:.3e}")
    return st


def invert(st, lam, maxiter=60, x0=None):
    M, res = solve(st["b"], st["mv"], st["rv"], st["reg"]["G_sub"], st["reg"]["w"],
                   EPS, lam, len(st["midx"]), maxiter=maxiter, x0=x0)
    return M, res


def evaluate(st, M):
    m = metrics(M, st["M_true"], st["vol"], st["b"], st["mv"], level=st["level"])
    m["res_rms"] = float(np.sqrt(np.mean((st["mv"](M.ravel()) - st["b"]) ** 2)))
    m["bias_max"] = max(abs(v - 1.0) for v in m["bias_level"].values())
    return m


def lambda_equiv_table(st):
    """同一物理强度下各约定所需的 λ（以 volume_raw 为 1）。

    volume_norm 权重 V_f/W、volume_m3 权重 V_f/1e9、none 权重 1 ⇒
        λ_norm = λ_raw·W,  λ_m3 = λ_raw·1e-9,  λ_none = λ_raw·W/nf
    """
    W = st["reg"]["W"]
    nf = st["reg"]["n_faces"]
    return {"volume_raw": 1.0, "volume_norm": W, "volume_m3": 1e-9, "none": W / nf,
            "_W": W, "_nf": nf, "_mean_w": W / nf}


def prod_lambda_equiv(st, variant):
    """生产 λ=1e-10（volume_raw 标度）换算到当前 variant 标度下的等效 λ。"""
    return LAMBDA_PROD * lambda_equiv_table(st)[variant]


# ======================================================================================
# E0 一致性锚定
# ======================================================================================
def E0(quick):
    print("\n=== E0 harness 与生产公式一致性 ===")
    domain, mbox = default_geom()
    nbase = (16, 16, 4) if quick else (32, 32, 8)
    st = setup(domain, mbox, "split", nbase, "volume_raw", "sparse", "uniform")
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1e4, 3 * len(st["midx"]))
    l1, g1 = objective_and_gradient(x, st["b"], st["mv"], st["rv"],
                                    st["reg"]["G_sub"], st["reg"]["w"], EPS, LAMBDA_PROD)
    l2, g2 = production_objective(x, st["b"], st["mv"], st["rv"],
                                  st["reg"]["G_sub"], st["reg"]["w_face"], EPS, LAMBDA_PROD)
    dl = abs(l1 - l2) / max(abs(l2), 1e-30)
    dg = float(np.max(np.abs(g1 - g2)) / max(np.max(np.abs(g2)), 1e-30))
    log("E0", "loss 相对差 < 1e-12", dl < 1e-12, f"Δloss/loss={dl:.2e}")
    log("E0", "gradient 相对差 < 1e-12", dg < 1e-12, f"Δg/g={dg:.2e}")
    return dict(dl=dl, dg=dg)


# ======================================================================================
# E1 离散化保真度
# ======================================================================================
def E1(quick):
    print("\n=== E1 体积权重是否为正确离散化（加密收敛）===")
    L = 4.0

    def Mfun(X, Y, Z):
        r2 = X ** 2 + Y ** 2 + Z ** 2
        w = np.where(r2 < 9.0, (1.0 - r2 / 9.0) ** 2, 0.0)
        return np.stack([w, w * 0.5, w * 0.25], axis=-1)

    def dM2(X, Y, Z):
        r2 = X ** 2 + Y ** 2 + Z ** 2
        f = lambda c: np.where(r2 < 9.0, -4.0 * c * (1.0 - r2 / 9.0) / 9.0, 0.0)
        return 1.3125 * (f(X) ** 2 + f(Y) ** 2 + f(Z) ** 2)

    from discretize import TreeMesh
    rows = []
    for nbase in ([16, 32] if quick else [16, 32, 64]):
        h = 2 * L / nbase
        mesh = TreeMesh([[(h, nbase)], [(h, nbase)], [(h, nbase)]],
                        origin=[-L, -L, -L], diagonal_balance=False)
        ml = mesh.max_level
        mesh.refine_bounding_box(np.array([[-L, -L, -L], [L, L, L]]),
                                 level=ml - 1, finalize=False)
        mesh.refine(lambda c: ml if c.center[0] < 0 else ml - 1, finalize=False)
        mesh.finalize(); mesh.number()
        cc = mesh.cell_centers
        V = mesh.cell_volumes
        G = mesh.cell_gradient.tocsr()
        w_all = mesh.get_face_inner_product().diagonal()
        Mv = Mfun(cc[:, 0], cc[:, 1], cc[:, 2])
        ref = float(np.sum(V * dM2(cc[:, 0], cc[:, 1], cc[:, 2])))
        g = G @ Mv
        interior = np.diff(G.indptr) >= 2
        Ew = float(np.sum(w_all[interior, None] * g[interior] ** 2))
        Eu = float(np.sum(g[interior] ** 2))
        rows.append(dict(nbase=nbase, nC=mesh.nC, ref=ref, E_w=Ew, E_u=Eu,
                         err_w=abs(Ew - ref) / ref, err_u=abs(Eu - ref) / ref))
        print(f"  nbase={nbase:3d} nC={mesh.nC:7d} 解析={ref:8.3f} "
              f"体积加权={Ew:8.3f}({abs(Ew-ref)/ref*100:6.2f}%) "
              f"不加权={Eu:11.3f}({abs(Eu-ref)/ref*100:9.2f}%)")
    order = np.log(rows[0]["err_w"] / rows[-1]["err_w"]) / np.log(rows[-1]["nbase"] / rows[0]["nbase"])
    log("E1", "体积加权按 O(h^1.2+) 收敛", order >= 1.2,
        f"观测阶数={order:.2f}（{rows[0]['err_w']*100:.2f}% -> {rows[-1]['err_w']*100:.2f}%）")
    log("E1", "最细网格体积加权误差 < 2%", rows[-1]["err_w"] < 0.02,
        f"{rows[-1]['err_w']*100:.2f}%")
    log("E1", "不加权明显错误 (>50%)", rows[-1]["err_u"] > 0.5, f"{rows[-1]['err_u']*100:.0f}%")

    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    nb = [r["nbase"] for r in rows]
    ax.loglog(nb, [r["err_w"] for r in rows], "o-", label="体积加权 w=S·h")
    ax.loglog(nb, [r["err_u"] for r in rows], "s-", label="不加权（旧写法）")
    ax.loglog(nb, [rows[0]["err_w"] * (nb[0] / n) ** 2 for n in nb], "k--", lw=1,
              label="O(h²) 参考斜率")
    ax.set_xlabel("每方向基础网格数 nbase")
    ax.set_ylabel("相对误差 |离散化−解析| / 解析")
    ax.set_title("E1 梯度能量离散化保真度（非均匀网格）")
    ax.grid(alpha=0.3, which="both"); ax.legend()
    fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, "E1_discretization.png"), dpi=180)
    plt.close(fig)
    return dict(rows=rows, order=order)


# ======================================================================================
# E2a 标定 λ 下的还原度与 λ 换算；E2b 欠正则化指纹
# ======================================================================================
def E2(quick):
    print("\n=== E2a 标定 λ 下的还原度 + 各约定 λ 换算表（split/dense/均匀块）===")
    domain, mbox = default_geom()
    nbase = (16, 16, 4) if quick else (32, 32, 8)
    ks = [-4, -3, -2, -1, 0, 1] if quick else [-4, -3, -2, -1, 0, 1, 2, 3]
    table, sts = {}, {}
    for variant in REG_VARIANTS:
        st = setup(domain, mbox, "split", nbase, variant, "dense", "uniform")
        sts[variant] = st
        eq = lambda_equiv_table(st)
        lam_prod = prod_lambda_equiv(st, variant)
        k_prod = float(np.log10(max(lam_prod, 1e-300) / st["lam_bal"]))
        print(f"  [{variant:11s}] λ_bal={st['lam_bal']:.3e}  W=Σw_f={eq['_W']:.4g}  "
              f"W/nf={eq['_mean_w']:.4g}  生产 λ=1e-10 等效={lam_prod:.2e} "
              f"(=λ_bal×10^{k_prod:+.1f})")
        rows = []
        for k in ks:
            lam = st["lam_bal"] * 10.0 ** k
            M, res = invert(st, lam)
            m = evaluate(st, M)
            rows.append(dict(k=k, lam=lam, rel_err=m["rel_err"], maxM=m["maxM_ratio"],
                             bias=m["bias_level"], bias_max=m["bias_max"],
                             res_rms=m["res_rms"], nit=res.nit, k_prod=k_prod))
            bs = " ".join(f"档{a}:{b:.2f}" for a, b in sorted(m["bias_level"].items()))
            print(f"      k={k:+d} λ={lam:9.2e} 还原误差={m['rel_err']*100:6.2f}% "
                  f"max|M|比={m['maxM_ratio']:5.2f} 残差={m['res_rms']:.2e}  {bs}")
        table[variant] = rows

    tol = 0.05 if quick else 0.02
    # 安全窗口取 k∈[0,+1]（λ 在 λ_bal 的 1~10 倍）：这是粗网格实测安全区 [-2,+1] 与
    # 细网格实测安全区 [0,+3] 的**交集**。窗口是单边的 —— λ 低于 λ_bal 一侧很快变差
    # （细网格 k=-1 就到 9%，k=-2 到 46%），所以"宁大勿小"，别往欠正则化方向推。
    win = [r for r in table["volume_norm"] if 0 <= r["k"] <= 1]
    worst_err = max(r["rel_err"] for r in win)
    worst_bias = max(r["bias_max"] for r in win)
    safe = [r["k"] for r in table["volume_norm"] if r["rel_err"] <= tol and r["bias_max"] <= 0.05]
    log("E2a", f"volume_norm 在 k∈[0,+1] 窗口还原误差 ≤ {tol*100:.0f}%", worst_err <= tol,
        f"最差={worst_err*100:.2f}%；实测安全 k={safe}")
    log("E2a", "volume_norm 在 k∈[0,+1] 窗口各档偏置 ≤ 0.05", worst_bias <= 0.05,
        f"最大|bias-1|={worst_bias:.3f}")

    st_n, st_r = sts["volume_norm"], sts["volume_raw"]
    lam_n = st_n["lam_bal"]
    lam_r = lam_n / lambda_equiv_table(st_r)["volume_norm"]
    Mn, _ = invert(st_n, lam_n)
    Mr, _ = invert(st_r, lam_r)
    fn, _ = objective_and_gradient(Mn.ravel(), st_n["b"], st_n["mv"], st_n["rv"],
                                   st_n["reg"]["G_sub"], st_n["reg"]["w"], EPS, lam_n)
    fr, _ = objective_and_gradient(Mn.ravel(), st_r["b"], st_r["mv"], st_r["rv"],
                                   st_r["reg"]["G_sub"], st_r["reg"]["w"], EPS, lam_r)
    dloss = abs(fn - fr) / max(abs(fr), 1e-30)
    dn = float(np.max(np.abs(Mn - Mr)) / max(np.max(np.abs(Mn)), 1e-30))
    log("E2a", "等效 λ 下两约定目标函数相同 (<1e-9)", dloss < 1e-9, f"Δf/f={dloss:.2e}")
    log("E2a", "等效 λ 下两约定解一致 (<1e-2, 仅求解路径差)", dn < 1e-2, f"max相对差={dn:.2e}")

    st_m3 = sts["volume_m3"]
    k_m3 = float(np.log10(max(prod_lambda_equiv(st_m3, "volume_m3"), 1e-300) / st_m3["lam_bal"]))
    log("E2a", "体积写成 m³ 又沿用生产 λ → 落到欠正则化区 (k ≤ -1)", k_m3 <= -1.0,
        f"k={k_m3:+.1f}（即 λ 只有该有的 10^{k_m3:.0f}）")

    print("\n=== E2b 欠正则化指纹：不加权 vs 体积加权（split/sparse/均匀块）===")
    sts_b, fp = {}, {}
    ks_b = [-10, -8, -6, -4, -2, 0] if quick else [-11, -10, -9, -8, -7, -6, -4, -2, 0]
    for variant in ("none", "volume_norm"):
        st = setup(domain, mbox, "split", nbase, variant, "sparse", "uniform", noise_frac=0.03)
        sts_b[variant] = st
        rows = []
        for k in ks_b:
            M, _ = invert(st, st["lam_bal"] * 10.0 ** k, maxiter=400)
            m = evaluate(st, M)
            rows.append(dict(k=k, rel_err=m["rel_err"], bias=m["bias_level"],
                             bias_max=m["bias_max"]))
            print(f"  {variant:11s} k={k:+3d} 还原误差={m['rel_err']*100:7.2f}% "
                  f"max|bias-1|={m['bias_max']:.3f}  "
                  + " ".join(f"档{a}:{b:.2f}" for a, b in sorted(m["bias_level"].items())))
        fp[variant] = rows
    # 取"指纹刚出现"的那个 λ 作对比：即不加权写法偏置首次 ≥ 0.1 的**最大** k
    # （在更深的欠正则化区两者都被最小范数解主导、正则化已不起作用，比在那里没意义）
    cand = [r for r in fp["none"] if r["bias_max"] >= 0.1]
    k_worst = max(r["k"] for r in cand) if cand else max(fp["none"], key=lambda r: r["bias_max"])["k"]
    r_none = next(r for r in fp["none"] if r["k"] == k_worst)
    r_norm = next(r for r in fp["volume_norm"] if r["k"] == k_worst)
    deep_none = max(fp["none"], key=lambda r: r["bias_max"])
    deep_norm = next(r for r in fp["volume_norm"] if r["k"] == deep_none["k"])
    log("E2b", "欠正则化时出现「细档偏小、粗档偏大」幅值指纹 (|bias-1| ≥ 0.1)",
        r_none["bias_max"] >= 0.1,
        f"指纹起始 @k={k_worst:+d}: " + " ".join(f"档{a}:{b:.2f}" for a, b in sorted(r_none["bias"].items())))
    log("E2b", "同一物理强度下体积加权把起始处偏置压到 ≤ 0.7×",
        r_norm["bias_max"] <= 0.7 * r_none["bias_max"],
        f"k={k_worst:+d}: {r_norm['bias_max']:.3f} vs {r_none['bias_max']:.3f} "
        f"(比值 {r_norm['bias_max']/max(r_none['bias_max'],1e-9):.2f})")
    print(f"  [注] 更深的欠正则化区（k={deep_none['k']:+d}）两者都被最小范数解主导："
          f"偏置 narrow={deep_none['bias_max']:.2f} vs volume={deep_norm['bias_max']:.2f}，"
          f"此时正则化已不起作用，属于「λ 太小」本身的问题")

    fig, axes = plt.subplots(1, 2, figsize=(12.6, 4.6))
    for variant, rows in table.items():
        ls = "-" if variant == "volume_norm" else "--"
        axes[0].semilogy([r["k"] for r in rows], [max(r["rel_err"], 1e-5) for r in rows],
                         "o" + ls, label=variant)
    axes[0].axvspan(-2, 2, color="tab:green", alpha=0.12, label="安全窗口 k∈[-2,2]")
    axes[0].set_xlabel("log10(λ/λ_bal)"); axes[0].set_ylabel("体积加权还原误差")
    axes[0].set_title("E2a 还原度 vs 正则化强度（split/dense）")
    axes[0].grid(alpha=0.3); axes[0].legend(fontsize=8)
    lv_names = ["细档", "次细档", "粗档", "更粗档"]
    for variant, rows in fp.items():
        for j, lv in enumerate(sorted(rows[0]["bias"])):
            axes[1].plot([r["k"] for r in rows], [r["bias"].get(lv, np.nan) for r in rows],
                         ("o-" if j % 2 == 0 else "s--"),
                         label=f"{variant} {lv_names[min(lv, 3)]}")
    axes[1].axhline(1.0, color="k", lw=0.8)
    axes[1].axvspan(-2, 2, color="tab:green", alpha=0.12)
    axes[1].set_xlabel("log10(λ/λ_bal)"); axes[1].set_ylabel("各档 |M| / 真值")
    axes[1].set_title("E2b 幅值指纹：欠正则化 → 细小粗大")
    axes[1].grid(alpha=0.3); axes[1].legend(fontsize=7)
    fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, "E2_lambda_sweep.png"), dpi=180)
    plt.close(fig)
    return dict(table=table, sts=sts, sts_b=sts_b, fp=fp, k_worst=k_worst,
                r_none=r_none, r_norm=r_norm, k_m3=k_m3,
                equiv=lambda_equiv_table(sts["volume_raw"]))


# ======================================================================================
# E3 带接缝块
# ======================================================================================
def _width_cells(st, M, y_seam):
    """有符号 Mx 剖面沿 y 从 10% 到 90% 的过渡长度（以 y 方向单元尺寸为单位）。"""
    order = np.argsort(st["cc"][:, 1])
    y, mx = st["cc"][order, 1], np.asarray(M)[order, 0]
    v_hi, v_lo = np.percentile(mx, 95), np.percentile(mx, 5)
    if v_hi - v_lo <= 0:
        return float("nan")
    l90, l10 = v_lo + 0.9 * (v_hi - v_lo), v_lo + 0.1 * (v_hi - v_lo)
    a_hi, a_lo = np.where(mx >= l90)[0], np.where(mx <= l10)[0]
    ah, al = a_hi[y[a_hi] < y_seam], a_lo[y[a_lo] > y_seam]
    if not (len(ah) and len(al)):
        return float("nan")
    cell = float(np.median(np.diff(np.unique(np.round(y, 3)))))
    return abs(y[al].min() - y[ah].max()) / max(cell, 1e-9)


def E3(quick, e2):
    print("\n=== E3 带接缝块（M 反向拼接：对比度/边缘宽度/按档公平性）===")
    domain, mbox = default_geom()
    nbase = (16, 16, 4) if quick else (32, 32, 8)
    y_seam = 0.5 * (mbox[0][1] + mbox[1][1])
    ks = [0, -2, -4] if quick else [0, -1, -2, -3, -4]
    out = {}
    for variant in ("volume_norm", "volume_raw", "none"):
        st = setup(domain, mbox, "split", nbase, variant, "dense", "seam")
        lo, hi = st["cc"][:, 1] < y_seam, st["cc"][:, 1] >= y_seam
        jump_true = st["data"]["truth_meta"]["m_left"] - st["data"]["truth_meta"]["m_right"]

        def jump(mx, m_lo, m_hi):
            return float(mx[m_lo].mean() - mx[m_hi].mean())

        Md, _ = invert(st, st["lam_bal"] * 1e-12, maxiter=400)   # 纯数据拟合（几乎无正则化）
        md = evaluate(st, Md)
        j_data = jump(Md[:, 0], lo, hi)
        width_data = _width_cells(st, Md, y_seam)

        best = None
        for k in ks:
            M, _ = invert(st, st["lam_bal"] * 10.0 ** k, maxiter=300)
            m = evaluate(st, M)
            per_level = {}
            for L in np.unique(st["level"]):
                sel = st["level"] == L
                per_level[int(L)] = (jump(M[sel, 0], lo[sel], hi[sel])
                                     if (sel & lo).any() and (sel & hi).any() else np.nan)
            bias = m["bias_level"]
            vals = [bias[a] for a in sorted(bias)]
            bias_spread = (max(vals) - min(vals)) if len(vals) > 1 else 0.0
            width_cells = _width_cells(st, M, y_seam)
            rec = dict(k=k, rel_err=m["rel_err"], bias=bias, bias_spread=bias_spread,
                       per_level=per_level, jump_rec=jump(M[:, 0], lo, hi),
                       width_cells=width_cells, maxM=m["maxM_ratio"])
            print(f"  {variant:11s} k={k:+d} 还原误差={m['rel_err']*100:5.1f}% "
                  f"(纯数据 {md['rel_err']*100:5.1f}%)  跳变={rec['jump_rec']/1e5:+.2f}/"
                  f"{jump_true/1e5:+.2f} ×1e5  各档偏置={'/'.join(f'{b:.2f}' for b in vals)} "
                  f"档间差={bias_spread:.3f}  边缘宽度={width_cells:.1f} 单元")
            if best is None or rec["rel_err"] < best["rel_err"]:
                best = rec
        best.update(jump_true=jump_true, jump_data=j_data, rel_err_data=md["rel_err"],
                    width_data=width_data, lam_bal=st["lam_bal"], ks=ks)
        out[variant] = best
    r, rn = out["volume_norm"], out["none"]
    print(f"  → 取各自最优 λ：volume_norm k={r['k']:+d}, none k={rn['k']:+d}")
    log("E3", "volume_norm 各档偏置之差 ≤ 0.05（体积加权按体积公平）",
        r["bias_spread"] <= 0.05, f"档间差={r['bias_spread']:.3f}（none={rn['bias_spread']:.3f}）")
    log("E3", "volume_norm 按档不公 ≤ 0.5× 不加权",
        r["bias_spread"] <= 0.5 * rn["bias_spread"],
        f"{r['bias_spread']:.3f} vs {rn['bias_spread']:.3f}")
    log("E3", "最优 λ 下 volume_norm 还原误差 ≤ 纯数据拟合", r["rel_err"] <= r["rel_err_data"],
        f"{r['rel_err']*100:.1f}% vs {r['rel_err_data']*100:.1f}%")
    log("E3", "最优 λ 下 volume_norm 不差于不加权", r["rel_err"] <= rn["rel_err"] + 0.01,
        f"{r['rel_err']*100:.1f}% vs {rn['rel_err']*100:.1f}%")
    log("E3", "接缝跳变被捕捉到（≥40% 真值）", r["jump_rec"] >= 0.4 * r["jump_true"],
        f"{r['jump_rec']/1e5:+.2f} / {r['jump_true']/1e5:+.2f} ×1e5")
    lv = [v for v in r["per_level"].values() if np.isfinite(v)]
    spread = (max(lv) - min(lv)) / abs(r["jump_rec"]) if len(lv) > 1 and r["jump_rec"] else np.nan
    log("E3", "各档独立估出的跳变离散 ≤ 20%", np.isfinite(spread) and spread <= 0.20,
        f"档间离散={spread*100:.1f}%")
    # 绝对宽度受"幅值整体收缩"影响很大（缩得越狠、剖面越像缓坡），因此判据用**相对**比较：
    # 正则化后的过渡宽度不应比纯数据拟合更宽（这才叫"边缘保持"）。
    log("E3", "正则化没把边缘抹得比纯数据拟合更宽（≤1.2×）",
        np.isfinite(r["width_cells"]) and np.isfinite(r["width_data"])
        and r["width_cells"] <= 1.2 * r["width_data"],
        f"正则化 {r['width_cells']:.0f} 单元 vs 纯数据 {r['width_data']:.0f} 单元"
        f"（none: {rn['width_cells']:.0f} / {rn['width_data']:.0f}）")
    return out


# ======================================================================================
# E4 网格无关性
# ======================================================================================
def E4(quick, e2):
    print("\n=== E4 网格无关性（同一数据、同一个 λ、三张网格）===")
    domain, mbox = default_geom()
    meshes = [("uniform_fine", "uniform", (32, 32, 8)),
              ("uniform_base", "uniform", (16, 16, 4)),
              ("split", "split", (16, 16, 4))]
    if quick:
        meshes = meshes[1:]
    out = {}
    for variant in ("volume_norm", "none"):
        st_ref = setup(domain, mbox, "split", (16, 16, 4), variant, "dense", "uniform")
        lam = st_ref["lam_bal"]
        rows = {}
        for name, style, nbase in meshes:
            st = setup(domain, mbox, style, nbase, variant, "dense", "uniform")
            M, _ = invert(st, lam)
            m = evaluate(st, M)
            rows[name] = dict(rel_err=m["rel_err"], bias=m["bias_level"],
                              bias_max=m["bias_max"], nC=st["mesh"].nC,
                              ncell=len(st["midx"]), lam_bal_own=st["lam_bal"])
            print(f"  {variant:11s} {name:14s} nC={st['mesh'].nC:6d} 磁体单元={len(st['midx']):5d} "
                  f"λ_bal(本网格)={st['lam_bal']:.2e} 还原误差={m['rel_err']*100:5.2f}%  "
                  + " ".join(f"档{a}:{b:.2f}" for a, b in sorted(m["bias_level"].items())))
        out[variant] = dict(rows=rows, lam=lam, lam_bal_ref=st_ref["lam_bal"])
    for variant, d in out.items():
        es = [r["rel_err"] for r in d["rows"].values()]
        spread = max(es) / max(min(es), 1e-12)
        bmax = max(r["bias_max"] for r in d["rows"].values())
        if variant == "volume_norm":
            log("E4", "volume_norm 跨网格误差变化 ≤ 1.5x（同一 λ 未重调）", spread <= 1.5,
                f"spread={spread:.2f} err={[f'{e*100:.2f}%' for e in es]}")
            log("E4", "volume_norm 跨网格各档偏置 ≤ 0.05", bmax <= 0.05,
                f"最大|bias-1|={bmax:.3f}")
        else:
            log("E4", "不加权写法跨网格变化更明显 (≥1.5x)", spread >= 1.5,
                f"spread={spread:.2f} err={[f'{e*100:.2f}%' for e in es]}")
    return out


# ======================================================================================
# E5 λ 选择规则
# ======================================================================================
def E5(quick, e2):
    print("\n=== E5 λ 选择：discrepancy principle vs 上帝视角最优 ===")
    rows = sorted(e2["table"]["volume_norm"], key=lambda r: r["k"])
    st = e2["sts"]["volume_norm"]
    sig = st["data"]["sigma"]
    errs = np.array([r["rel_err"] for r in rows])
    i_oracle = int(np.argmin(errs))
    # discrepancy：**满足残差≈σ 的最大 λ**（最正则化但仍与噪声相符的解）
    ok = np.array([r["res_rms"] <= 1.05 * sig for r in rows])
    i_disc = int(np.where(ok)[0].max()) if ok.any() else int(np.argmin(np.abs(
        np.array([r["res_rms"] for r in rows]) - sig)))
    lam_o, lam_d = rows[i_oracle]["lam"], rows[i_disc]["lam"]
    rl = max(lam_d / lam_o, lam_o / lam_d)
    re = errs[i_disc] / max(errs[i_oracle], 1e-12)
    log("E5", "自动 λ 落在最优 λ 的 3 倍以内", rl <= 3.0,
        f"k*={rows[i_oracle]['k']:+d} vs k_auto={rows[i_disc]['k']:+d} (x{rl:.2f})")
    log("E5", "自动 λ 的误差 ≤ 最优误差的 1.3 倍", re <= 1.3,
        f"{errs[i_disc]*100:.2f}% vs {errs[i_oracle]*100:.2f}%")
    fig, ax = plt.subplots(figsize=(6.8, 4.6))
    ks = [r["k"] for r in rows]
    ax.semilogy(ks, errs, "o-", label="还原误差")
    ax.semilogy(ks, [r["res_rms"] / sig for r in rows], "s--", label="数据残差 / 噪声 σ")
    ax.axhline(1.05, color="tab:green", ls=":", lw=1, label="残差≈σ 判据 1.05")
    ax.axvline(rows[i_oracle]["k"], color="tab:red", ls="--", label=f"最优 k={rows[i_oracle]['k']:+d}")
    ax.axvline(rows[i_disc]["k"], color="tab:green", ls=":", label=f"自动 k={rows[i_disc]['k']:+d}")
    ax.set_xlabel("log10(λ/λ_bal)"); ax.set_ylabel("相对误差 / 残差比")
    ax.set_title("E5 λ 选择规则（split/dense, volume_norm）")
    ax.grid(alpha=0.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(FIG_DIR, "E5_lambda_rule.png"), dpi=180)
    plt.close(fig)
    return dict(k_oracle=rows[i_oracle]["k"], k_disc=rows[i_disc]["k"],
                err_oracle=errs[i_oracle], err_disc=errs[i_disc], sigma=sig,
                lam_disc=lam_d, lam_bal=st["lam_bal"], equiv=lambda_equiv_table(st),
                k_ratio=rl)


# ======================================================================================
# E6 求解器收敛性
# ======================================================================================
def E6(quick, e2):
    print("\n=== E6 求解器收敛性（maxiter 60 vs 400）===")
    st = e2["sts"]["volume_norm"]
    lam = st["lam_bal"]
    M60, r60 = invert(st, lam, maxiter=60)
    M400, r400 = invert(st, lam, maxiter=400, x0=r60.x)
    m60, m400 = evaluate(st, M60), evaluate(st, M400)
    d = abs(m60["rel_err"] - m400["rel_err"]) / max(m400["rel_err"], 1e-12)
    log("E6", "maxiter=60 与 400 误差差 < 1%", d < 0.01,
        f"nit={r60.nit}->{r400.nit}, 误差 {m60['rel_err']*100:.2f}% vs {m400['rel_err']*100:.2f}%")
    log("E6", "目标函数相对改善 < 1%", abs(r60.fun - r400.fun) / max(abs(r400.fun), 1e-30) < 0.01,
        f"f={r60.fun:.6e} -> {r400.fun:.6e}")
    return dict(nit60=r60.nit, nit400=r400.nit, err60=m60["rel_err"], err400=m400["rel_err"],
                fun60=r60.fun, fun400=r400.fun)


# ======================================================================================
# E7 粗细网格互换不变性（唯一判据：消除粗/细网格的区别，且对调粗细不改变结果）
# ======================================================================================
SWAP_TOL = 0.02          # 互换后同一区域均值的相对变化上限
RATIO_BAND = (0.98, 1.02)  # 每张网格上"左区均值/右区均值"的允许区间


def E7(quick):
    """同一批数据（均匀磁化真值、无噪声、λ=0 完全不带正则化），

    两张**互为粗细互换**的网格（A：左半细；B：右半细），
    比较左/右半区的体积加权平均 |m|：
      * 判据：√v·m 变量下，互换两网格后同一区域均值变化 ≤2%，
              且每张网格上"左/右"比值落在 [0.98,1.02]；
      * 对照：同样条件下 m 变量（只打印，不作判据）。
    """
    print("\n=== E7 粗细网格互换不变性（λ=0，无正则化）===")
    domain, mbox = default_geom()
    nbase = (16, 16, 4)
    x_split = 0.5 * domain[0]

    data = T.synthesize("uniform", domain, mbox, meas_mode="dense", step=4.0,
                        noise_frac=0.0, seed=0, verbose=False)
    b = data["B_meas"].ravel()
    print(f"  真值：整块均匀磁化 {T.M0_DEFAULT:.0e} A/m；测点 {len(data['meas'])} 个、"
          f"无噪声；λ=0（目标函数里没有正则化项）")

    # 三组"对四张求解完全一致"的设置
    settings = [("LSQR@50", "lsqr", 50), ("LSQR@200", "lsqr", 200), ("L-BFGS-B@50", "lbfgs", 50)]

    cases = {}       # (setting, variable) -> {"A":..., "B":...}
    for sname, skind, itlim in settings:
        for variable in ("m", "sqrt_v_m"):
            for side in ("low", "high"):
                mesh = build_mesh(domain, nbase, style="split", magnet_box=mbox,
                                  split_axis=0, split_frac=0.5, mag_levels=1, fine_side=side)
                midx = magnet_cells(mesh, mbox)
                cc = mesh.cell_centers[midx]
                vol = mesh.cell_volumes[midx]
                vol_SI = vol * 1e-9
                scale = variable_scale(vol_SI, variable)
                coef = operator_coef(vol_SI, variable)
                mv, rv = make_forward(data["meas"] / 1000.0, cc / 1000.0, coef)
                if skind == "lsqr":
                    A = LinearOperator((len(b), 3 * len(cc)), matvec=mv, rmatvec=rv,
                                       dtype=float)
                    x = lsqr(A, b, atol=0.0, btol=0.0, iter_lim=itlim)[0]
                else:
                    reg = build_regularizer(mesh, np.isin(np.arange(mesh.nC), midx),
                                            variant="volume_raw")
                    x = solve(b, mv, rv, reg["G_sub"], np.zeros(reg["n_faces"]), EPS, 0.0,
                              len(cc), maxiter=itlim, scale=scale)[0].ravel()
                mag = np.linalg.norm(to_m(x, scale), axis=1)
                left = cc[:, 0] < x_split
                mean_L = float(np.sum(vol[left] * mag[left]) / np.sum(vol[left]))
                mean_R = float(np.sum(vol[~left] * mag[~left]) / np.sum(vol[~left]))
                cases.setdefault((sname, variable), {})[side] = dict(mean_L=mean_L, mean_R=mean_R)

    # 打印 + 判据（只判 √v·m）
    table = []
    ok_all = True
    for sname, _, _ in settings:
        for variable in ("m", "sqrt_v_m"):
            A_ = cases[(sname, variable)]["low"]     # 左半细
            B_ = cases[(sname, variable)]["high"]    # 右半细
            inv_L = A_["mean_L"] / B_["mean_L"]
            inv_R = A_["mean_R"] / B_["mean_R"]
            swap_err = max(abs(inv_L - 1.0), abs(inv_R - 1.0))
            ratio_A = A_["mean_L"] / A_["mean_R"]
            ratio_B = B_["mean_L"] / B_["mean_R"]
            row = dict(setting=sname, variable=variable, mean_L_A=A_["mean_L"],
                       mean_R_A=A_["mean_R"], mean_L_B=B_["mean_L"], mean_R_B=B_["mean_R"],
                       ratio_A=ratio_A, ratio_B=ratio_B, swap_err=swap_err)
            table.append(row)
            flag = "√v·m" if variable == "sqrt_v_m" else "m   "
            print(f"  [{sname:11s}] {flag}: 网格A(左细) 左={A_['mean_L']:9.0f} 右={A_['mean_R']:9.0f}"
                  f" 左/右={ratio_A:5.2f} | 网格B(右细) 左={B_['mean_L']:9.0f} 右={B_['mean_R']:9.0f}"
                  f" 左/右={ratio_B:5.2f} | 互换不变性误差={swap_err*100:5.1f}%")
            if variable == "sqrt_v_m":
                ok = (swap_err <= SWAP_TOL
                      and RATIO_BAND[0] <= ratio_A <= RATIO_BAND[1]
                      and RATIO_BAND[0] <= ratio_B <= RATIO_BAND[1])
                ok_all = ok_all and ok
                log("E7", f"[{sname}] √v·m 互换不变性 ≤{SWAP_TOL*100:.0f}% 且左右比值 ∈ "
                          f"[{RATIO_BAND[0]},{RATIO_BAND[1]}]", ok,
                    f"互换误差={swap_err*100:.1f}%，比值 {ratio_A:.2f}/{ratio_B:.2f}")

    # 对照行（不作判据）：m 变量的最差互换误差
    worst_m = max(r["swap_err"] for r in table if r["variable"] == "m")
    print(f"  [对照] m 变量的互换不变性误差 = {worst_m*100:.1f}%（不作判据，仅记录）")

    # 内部健全性：梯度链式法则（有限差分）
    st_probe = None
    for side in ("low",):
        mesh = build_mesh(domain, nbase, style="split", magnet_box=mbox,
                          split_axis=0, split_frac=0.5, mag_levels=1, fine_side=side)
        midx = magnet_cells(mesh, mbox)
        cc = mesh.cell_centers[midx]; vol = mesh.cell_volumes[midx]
        reg = build_regularizer(mesh, np.isin(np.arange(mesh.nC), midx), variant="volume_raw")
        vol_SI = vol * 1e-9
        for variable in ("m", "sqrt_v_m"):
            scale = variable_scale(vol_SI, variable)
            mv, rv = make_forward(data["meas"] / 1000.0, cc / 1000.0,
                                  operator_coef(vol_SI, variable))
            # 用非零探测点：x=0 处 λ 项导数为 0，测不出数据项被多除 scale 的 bug
            xp = np.random.default_rng(0).normal(0.0, float(np.mean(scale)), 3 * len(cc))
            err = check_gradient(xp, b, mv, rv, reg["G_sub"], np.zeros(reg["n_faces"]),
                                 EPS, 0.0, scale=scale)
            print(f"  [健全性] 变量 {variable:9s} 梯度有限差分相对范数误差 = {err:.2e}")
    del st_probe
    return dict(table=table, worst_m=worst_m, ok=ok_all)


# ======================================================================================
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--full", action="store_true")
    ap.add_argument("--only", default=None)
    args = ap.parse_args(argv)
    quick = bool(args.quick)
    only = args.only

    os.makedirs(FIG_DIR, exist_ok=True)
    t0 = time.time()
    print(f"meshtest（{'quick' if quick else 'full'}）：λ 用力的平衡标定 λ_bal，"
          f"生产 λ={LAMBDA_PROD:g}(volume_raw/mm³)，ε={EPS:g}")

    if only in (None, "E0"):
        E0(quick)
    if only in (None, "E1"):
        E1(quick)

    e7 = None
    if only in (None, "E7"):
        e7 = E7(quick)

    e2 = None
    if only in (None, "E2", "E3", "E4", "E5", "E6"):
        e2 = E2(quick)
    if only in (None, "E3") and e2:
        E3(quick, e2)
    if only in (None, "E4") and e2:
        E4(quick, e2)
    if only in (None, "E5") and e2:
        E5(quick, e2)
    if only in (None, "E6") and e2:
        E6(quick, e2)

    dt = time.time() - t0
    write_results(dt, quick, e2, e7)
    n_fail = sum(1 for r in RESULTS if not r["ok"])
    print(f"\n完成：{len(RESULTS)} 条判据，{len(RESULTS)-n_fail} 通过，{n_fail} 失败，"
          f"耗时 {dt/60:.1f} min -> {RES_PATH}")
    return 0 if n_fail == 0 else 1


def write_results(dt, quick, e2, e7=None):
    lines = ["# meshtest 实验结果", "",
             f"- 生成时间：{time.strftime('%Y-%m-%d %H:%M:%S')}（{'quick' if quick else 'full'} 模式）",
             "- λ 标定：λ_bal = ‖∇E_data(M_ref)‖ / ‖∇E_reg(M_ref)‖；扫描 λ = λ_bal·10^k",
             f"- 生产代码：λ={LAMBDA_PROD:g}（volume_raw，w=V_f mm³），ε={EPS:g}",
             f"- 耗时：{dt/60:.1f} min", ""]
    if e2 is not None:
        eq = e2["equiv"]
        lines += ["## 同一物理强度下各约定所需的 λ（相对 volume_raw = 1）", "",
                  "| 约定 | λ 倍数 | 说明 |", "|---|---|---|",
                  "| volume_raw | 1 | 当前生产：w=V_f（mm³） |",
                  f"| volume_norm | {eq['volume_norm']:.4g} | 推荐：w=V_f/W，加密与单位都无关 |",
                  f"| volume_m3 | {eq['volume_m3']:.1e} | w=V_f（m³）：与 raw 差 10⁹ |",
                  f"| none | {eq['none']:.4g} | 旧写法：每面权重 1 |", "",
                  f"（本配置 W=Σw_f={eq['_W']:.4g} mm³，内部面 {eq['_nf']} 个，"
                  f"平均面权重 W/nf={eq['_mean_w']:.4g} mm³）", "",
                  "→ 加体积权重后，**同一个物理强度对应的 λ 要乘 W/nf**"
                  f"（本配置 ≈ {eq['_mean_w']:.0f} 倍；若把体积换成 m³ 则是 10⁹ 倍）。"
                  "这就是「加了体积项结果就非物理」的直接原因：λ 没跟着换算，"
                  "等效正则化强度掉到欠正则化区，磁化强度就会堆到大单元（粗网格）上。", ""]
    if e7 is not None:
        lines += ["## E7 粗细网格互换不变性（λ=0，无正则化；均匀真值 ⇒ 物理上左右均值应相等）", "",
                  "| 设置 | 变量 | 网格A(左细) 左/右 | 网格B(右细) 左/右 | 互换不变性误差 |",
                  "|---|---|---|---|---|"]
        for r in e7["table"]:
            lines.append(f"| {r['setting']} | {'√v·m' if r['variable']=='sqrt_v_m' else 'm'} | "
                         f"{r['ratio_A']:.2f} | {r['ratio_B']:.2f} | {r['swap_err']*100:.1f}% |")
        lines += [f"", f"对照：m 变量的最差互换误差 = {e7['worst_m']*100:.1f}%（不参与判据）。", ""]
    lines += ["## 判据汇总", "", "| 实验 | 判据 | 结果 | 数据 |", "|---|---|---|---|"]
    for r in RESULTS:
        lines.append(f"| {r['exp']} | {r['what']} | {'✅ PASS' if r['ok'] else '❌ FAIL'} | {r['detail']} |")
    n_fail = sum(1 for r in RESULTS if not r["ok"])
    lines += ["", f"**{len(RESULTS)-n_fail}/{len(RESULTS)} 条通过。**", ""]
    with open(RES_PATH, "w", encoding="utf-8") as f:
        f.write("\n".join(lines))
    print(f"已写入 {RES_PATH}")


if __name__ == "__main__":
    sys.exit(main())
