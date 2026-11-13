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


def run_demo(obs1, obs2, outdir=None, n=11, verbose=True,
             center=(-335.0, 71.0), half=180.0):
    outdir = outdir or os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "out")
    s = Session(obs1, obs2, outdir=outdir, verbose=verbose)

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
