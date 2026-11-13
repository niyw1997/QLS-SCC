"""数值自检：用已知三维点走一遍投影 → 极线 → 三角化 → PCA。

不需要界面，也不需要窗。默认拿本目录 data/20250302 的真实几何当试验台，
真值是一组人为指定的 HEEQ 点；只要链路正确，三角化必须把它们原样收回来。
"""

from __future__ import annotations

import numpy as np

from . import geometry as G
from .views import load_view

__all__ = ["run_selftest"]

DEFAULT_A = "examples/20250416/data/aia193_20250416_cutout.fits"
DEFAULT_B = "examples/20250416/data/euvi195_20250416_cutout.fts"


def _pass(name, value, tol, unit="", ok=None):
    good = (abs(value) <= tol) if ok is None else ok
    print(f"  [{'PASS' if good else 'FAIL'}] {name}: {value:.3e} {unit} (tol {tol:.1e})")
    return bool(good)


def run_selftest(obs1=None, obs2=None, verbose=True):
    import os
    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    obs1 = obs1 or os.path.join(here, DEFAULT_A)
    obs2 = obs2 or os.path.join(here, DEFAULT_B)
    a = load_view(obs1)
    b = load_view(obs2)
    frame = G.ref_frame(max(a.obstime, b.obstime))
    print(f"[selftest] {a.short} + {b.short}")
    print(f"  视角差 {G.separation_deg(a.m, b.m, frame):.3f}°   "
          f"极线灵敏度 {G.epipolar_sensitivity(a.m, b.m, frame):.1f}\"/0.1 R☉")

    ok = True
    # 1) 三角化回真值：在东西两侧各取几个点，含 r 1.02 / 1.10 / 1.20
    rng = np.random.default_rng(20261113)
    truth = []
    for r in (1.02, 1.08, 1.15, 1.22):
        for lon, lat in ((0.0, 0.0), (20.0, -20.0), (-25.0, 15.0), (10.0, -45.0)):
            truth.append((r, lon, lat))
    x = G.RSUN_KM * np.array([r * np.cos(np.radians(lat)) * np.cos(np.radians(lon))
                              for r, lon, lat in truth])
    y = G.RSUN_KM * np.array([r * np.cos(np.radians(lat)) * np.sin(np.radians(lon))
                              for r, lon, lat in truth])
    z = G.RSUN_KM * np.array([r * np.sin(np.radians(lat)) for r, lon, lat in truth])
    pts = np.column_stack([x, y, z])
    from astropy.coordinates import SkyCoord
    from astropy.coordinates import CartesianRepresentation
    from sunpy.coordinates import HeliographicStonyhurst
    import astropy.units as u
    sky = SkyCoord(CartesianRepresentation(*[c * u.km for c in pts.T]),
                   frame=HeliographicStonyhurst(obstime=frame.obstime))
    xa, ya = G.project_to(a.m, sky)
    xb, yb = G.project_to(b.m, sky)
    dr, dlon, dlat, gaps = [], [], [], []
    for i in range(len(truth)):
        t = G.triangulate(a.m, b.m, (xa[i], ya[i]), (xb[i], yb[i]), frame)
        dr.append(t["r"] - truth[i][0])
        dlon.append((t["lon"] - truth[i][1] + 180) % 360 - 180)
        dlat.append(t["lat"] - truth[i][2])
        gaps.append(t["gap_km"])
    print("  三角化（%d 个合成真值点）:" % len(truth))
    ok &= _pass("max |Δr|", float(np.max(np.abs(dr))), 2e-4, "R☉")
    ok &= _pass("max |Δlon|", float(np.max(np.abs(dlon))), 2e-3, "deg")
    ok &= _pass("max |Δlat|", float(np.max(np.abs(dlat))), 2e-3, "deg")
    ok &= _pass("max ray gap", float(np.max(gaps)), 5.0, "km")

    # 2) 极线一致性：真值投到 b 的点，必须落在 a 点发出的极线上
    dev = []
    for i in range(len(truth)):
        e = G.epipolar(a.m, b.m, xa[i], ya[i], frame, r_lo=truth[i][0] - 1e-4,
                       r_hi=truth[i][0] + 1e-4, n=3)
        dev.append(np.hypot(e["x"][1] - xb[i], e["y"][1] - yb[i]))
    print("  极线:")
    ok &= _pass("max 偏离", float(np.max(dev)), 2e-2, "arcsec")

    # 3) 视场联动：a 的四角 r=1.02 投影必须落在联动框内
    box = (-200.0, 400.0, -500.0, -100.0)
    lb = G.link_box(a.m, b.m, box, frame, r_lo=1.0, r_hi=1.10)
    corners = np.array([[-200, -500], [400, -500], [400, -100], [-200, -100]], float)
    ex, ey = G.project_to(b.m, G.common_point(a.m, corners[:, 0], corners[:, 1],
                                              1.02, frame))
    inside = ((ex >= lb[0]) & (ex <= lb[1]) & (ey >= lb[2]) & (ey <= lb[3]))
    print(f"  联动框: {tuple(round(v, 1) for v in lb)}")
    ok &= _pass("四角入框", float(np.sum(~inside)), 0.5, "个角", ok=inside.all())

    # 4) PCA：已知椭圆，点均匀采样
    ra_t, rb_t = 0.04, 0.09
    # 开弧：从一脚到另一脚，首末中点才是椭圆中心（闭合环会把中点拉到环边）
    th = np.linspace(-0.5 * np.pi, 0.5 * np.pi, 181)
    n = np.array([0.3, 0.5, 0.81])
    n /= np.linalg.norm(n)
    e1 = np.cross(n, [0, 0, 1.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(n, e1)
    ell = (ra_t * np.cos(th)[:, None] * e1 + rb_t * np.sin(th)[:, None] * e2) * G.RSUN_KM
    fit = G.pca_fit(ell)
    print("  PCA:")
    ok &= _pass("ra 相对误差", (fit["ra"] / G.RSUN_KM - ra_t) / ra_t, 1e-2, "",
                ok=abs(fit["ra"] / G.RSUN_KM - ra_t) / ra_t < 1e-2)
    ok &= _pass("rb 相对误差", (fit["rb"] / G.RSUN_KM - rb_t) / rb_t, 1e-2, "",
                ok=abs(fit["rb"] / G.RSUN_KM - rb_t) / rb_t < 1e-2)
    ok &= _pass("法向对齐", float(1 - abs(fit["evector"][:, 0] @ n)), 1e-6, "",
                ok=abs(abs(fit["evector"][:, 0] @ n) - 1) < 1e-6)

    print(f"[selftest] {'全部通过' if ok else '有失败项'}")
    return int(not ok)
