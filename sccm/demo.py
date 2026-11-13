"""无窗口演示：把界面上的两步流程脚本化。

演示用的是一段**已知的三维半环**，把它投到两幅真实像上，再当作"人工点选"
的结果逐点走界面同一条链路（pick_a -> epipolar -> pick_b -> triangulate）。
这样可以同时验收三件事：投影像上的位置对不对、极线是否穿过对映点、
三角化是否把三维点原样收回来。它不是对这颗暗条的测量结论，报告里会写明。

真实数据上的自动配对（``auto_pairs``）保留作辅助，判据是"极线在另一视角
轮廓上的最近点"，只用于快速预填，不能替代人工确认。
"""

from __future__ import annotations

import os

import numpy as np

from . import geometry as G
from . import render
from .app import Session

__all__ = ["run_demo", "synthetic_loop", "auto_pairs"]


def synthetic_loop(n=11, r0=1.07, lon0=-8.0, lat0=-26.0,
                   ra=0.028, rb=0.055, t0=0.0, t1=np.pi):
    """在 (lon0, lat0, r0) 的局部切平面上建一段半椭圆，返回 HEEQ 直角坐标 (N,3)，km。"""
    lon, lat = np.radians(lon0), np.radians(lat0)
    er = np.array([np.cos(lat) * np.cos(lon), np.cos(lat) * np.sin(lon), np.sin(lat)])
    elon = np.array([-np.sin(lon), np.cos(lon), 0.0])
    elat = np.cross(er, elon)
    c = r0 * G.RSUN_KM * er
    t = np.linspace(t0, t1, n)
    pts = (c[None, :]
           + ra * G.RSUN_KM * np.cos(t)[:, None] * elon[None, :]
           + rb * G.RSUN_KM * np.sin(t)[:, None] * elat[None, :])
    return pts


def _project_points(m, pts, frame):
    from astropy.coordinates import CartesianRepresentation, SkyCoord
    import astropy.units as u
    sky = SkyCoord(CartesianRepresentation(*[c * u.km for c in pts.T]), frame=frame)
    return G.project_to(m, sky)


def _window_around(a, b, xa, ya, xb, yb, pad=140.0):
    from .views import Window
    ca = (0.5 * (np.nanmin(xa) + np.nanmax(xa)), 0.5 * (np.nanmin(ya) + np.nanmax(ya)))
    half = max(np.nanmax(xa) - np.nanmin(xa), np.nanmax(ya) - np.nanmin(ya)) / 2 + pad
    return Window(ca[0], ca[1], half)


def run_demo(obs1, obs2, outdir=None, n=11, verbose=True, real=False,
             center=(200.0, -300.0), half=330.0):
    outdir = outdir or os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "out")
    s = Session(obs1, obs2, outdir=outdir, verbose=verbose)

    if real:
        s.zoom("a", (center[0] - half, center[0] + half,
                     center[1] - half, center[1] + half), link=True)
        msg = s.identify()
        if verbose:
            print("[demo] 识别:", msg)
        acc, rej = auto_pairs(s, n_max=n)
        if verbose:
            print(f"[demo] 自动配对：接受 {len(acc)}，拒绝 {len(rej)}")
        for pa, pb, r, dist, gap in acc:
            s.pending = dict(x=pa[0], y=pa[1], epi=G.epipolar(
                s.a.m, s.b.m, pa[0], pa[1], s.frame, r_lo=1.0, r_hi=2.0, n=161))
            s.pick_b(pb[0], pb[1])
        s.pending = None
    else:
        pts = synthetic_loop(n=n)
        xa, ya = _project_points(s.a.m, pts, s.frame)
        xb, yb = _project_points(s.b.m, pts, s.frame)
        ok = np.isfinite(xa) & np.isfinite(ya) & np.isfinite(xb) & np.isfinite(yb)
        pts, xa, ya, xb, yb = pts[ok], xa[ok], ya[ok], xb[ok], yb[ok]
        if len(pts) < 3:
            raise RuntimeError("合成环投不到两幅像上")
        wa = _window_around(s.a, s.b, xa, ya, xb, yb)
        s.set_window("a", wa, link=True)
        if verbose:
            print(f"[demo] 合成半环 {len(pts)} 点，视场 a {s.win['a']}  b {s.win['b']}")
        truth = pts / G.RSUN_KM
        err = []
        for i in range(len(pts)):
            s.pick_a(float(xa[i]), float(ya[i]))
            p = s.pick_b(float(xb[i]), float(yb[i]))
            if p is not None:
                err.append(np.linalg.norm(p.xyz - truth[i]))
        if verbose and err:
            print(f"[demo] 回真值最大偏差 {max(err):.2e} R_sun（{max(err)*G.RSUN_KM:.2f} km）")

    if verbose:
        sm = s.summary()
        print("[demo] 小结:", {k: (round(v, 4) if isinstance(v, float) else v)
                               for k, v in sm.items()})
    paths = s.save(make_figure=True, dpi=300)
    if verbose:
        print("[demo] 已保存:\n  " + "\n  ".join(paths))
    return s, paths


def auto_pairs(session, n_max=10, tol=12.0, r_lo=1.0, r_hi=1.3,
               gap_max_km=4000.0):
    """真实数据上的快速预填：沿参考脊线取点，找极线在另一视角轮廓上的最近点。"""
    if not session.contours["a"] or not session.contours["b"]:
        return [], []
    ca = session.contours["a"][0]
    cb = session.contours["b"][0]
    c = ca.get("spine")
    if c is None:
        c = ca["contour"]
    c = np.asarray(c, float)
    target = np.asarray(cb["contour"], float)
    step = max(1, len(c) // 60)
    idx = list(range(0, len(c), step))
    accepted, rejected = [], []
    for k in idx:
        pa = (float(c[k, 0]), float(c[k, 1]))
        epi = G.epipolar(session.a.m, session.b.m, pa[0], pa[1], session.frame,
                         r_lo=r_lo, r_hi=r_hi, n=121)
        d = np.hypot(target[None, :, 0] - epi["x"][:, None],
                     target[None, :, 1] - epi["y"][:, None])
        j, i = np.unravel_index(np.argmin(d), d.shape)
        dist = float(d[j, i])
        if dist > tol:
            rejected.append((pa, tuple(target[i]), float(epi["r"][j]), dist, np.nan))
            continue
        pb = (float(target[i, 0]), float(target[i, 1]))
        res = G.triangulate(session.a.m, session.b.m, pa, pb, session.frame)
        gap = float(res["gap_km"])
        row = (pa, pb, float(res["r"]), dist, gap)
        if gap <= gap_max_km:
            accepted.append(row)
        else:
            rejected.append(row)
    if len(accepted) > n_max:
        sel = np.linspace(0, len(accepted) - 1, n_max).astype(int)
        accepted = [accepted[i] for i in sel]
    return accepted, rejected
