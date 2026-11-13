"""暗结构识别与脊线。

以种子点为准做局部生长。AIA 304 的暗区常常连成一片，全局阈值会把整个活动区
吞进来；给了种子之后，阈值取种子暗度的一个比例，并限制在种子周围一个局部
窗口内，结果既可预期，也便于在界面上解释。

界面里的 Outline 现在在 AIA 与 EUVI 上各用一个种子独立识别，互不引导。
走廊相关函数保留给离线分析，主流程不再调用。
"""

from __future__ import annotations

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy import ndimage as ndi
from sunpy.coordinates import Heliocentric

from . import geometry as G

__all__ = ["contour_from_mask", "spine_from_mask", "detect", "corridor_points",
           "corridor_mask", "darkest_in_corridor"]


def contour_from_mask(mask):
    """marching squares 取最长闭合轮廓，返回 (col, row) 数组。"""
    from matplotlib.figure import Figure

    if not np.any(mask):
        return None
    fig = Figure()
    ax = fig.add_subplot(111)
    try:
        cs = ax.contour(np.asarray(mask, float), levels=[0.5])
        segs = cs.allsegs[0] if cs.allsegs else []
        if not segs:
            return None
        c = max(segs, key=len)
        return np.asarray(c, float).reshape(-1, 2)
    finally:
        fig.clf()


def spine_from_mask(mask, view, n=18):
    """用主成分方向把掩膜压成一条有序脊线，返回 HPC 角秒 (M,2)。

    这不是精确的骨架，但对"沿暗条布点"足够：把像素投到主轴，分箱取中位
    位置，箱的顺序就是沿结构的顺序。
    """
    ys, xs = np.where(mask)
    if len(xs) < 30:
        return None
    p = np.column_stack([xs, ys]).astype(float)
    c = p.mean(axis=0)
    q = p - c
    cov = q.T @ q / len(p)
    ev, evec = np.linalg.eigh(cov)
    u = evec[:, -1]
    t = q @ u
    lo, hi = np.percentile(t, [2, 98])
    if hi - lo < 4:
        return None
    edges = np.linspace(lo, hi, n + 1)
    out = []
    for i in range(n):
        sel = (t >= edges[i]) & (t < edges[i + 1])
        if sel.sum() < 6:
            continue
        out.append(p[sel].mean(axis=0))
    if len(out) < 3:
        return None
    out = np.asarray(out)
    x, y = view.to_world(out[:, 0], out[:, 1])
    return np.column_stack([x, y])


def _component(mask, lab, seed_ij=None):
    """有种子时只认种子所在的那一块；种子被阈值排除就返回 0。"""
    if seed_ij is not None:
        j, i = seed_ij
        return int(lab[j, i])
    sizes = ndi.sum(mask, lab, range(1, lab.max() + 1))
    if len(sizes) == 0:
        return 0
    return int(np.argmax(sizes)) + 1


def _grow(d, region, thr, seed_ij, min_area, bridge):
    """按阈值取连通域，桥接断口，返回布尔掩膜或 None。"""
    m = (d > thr) & region
    if m.sum() < min_area:
        return None
    m = ndi.binary_closing(m, np.ones((5, 5)), iterations=1)
    lab, n = ndi.label(m)
    if n == 0:
        return None
    k = _component(m, lab, seed_ij)
    if k == 0:
        return None
    sel = lab == k
    if sel.sum() < min_area:
        return None
    if bridge:
        keep = np.zeros_like(m)
        sizes = ndi.sum(m, lab, range(1, n + 1))
        for i2, s in enumerate(sizes, 1):
            if s >= max(min_area, 0.03 * sel.sum()):
                keep |= lab == i2
        br = ndi.binary_dilation(keep, ndi.generate_binary_structure(2, 2),
                                 iterations=bridge)
        lab2, n2 = ndi.label(br)
        k2 = _component(br, lab2, seed_ij)
        if k2 > 0:
            sel = ndi.binary_erosion(lab2 == k2, ndi.generate_binary_structure(2, 2),
                                     iterations=bridge)
            sel = ndi.binary_closing(sel, np.ones((7, 7)), iterations=1)
            if sel.sum() < min_area:
                sel = lab == k
    return sel


def detect(view, win, seed=None, thr=None, local_radius=180.0, sigma=45.0,
           percentile=92.0, min_area=120, bridge=6, fill=True, restrict=None,
           seed_frac=0.95, mean_frac=0.90, max_iter=5):
    """在窗口中找一块暗结构，围绕 ``seed`` 生长。

    除了按种子暗度的比例定阈值，还要求整块区域的**平均暗度不低于种子暗度的
    ``mean_frac`` 倍**。不满足就把阈值逐步抬高再取一次，直到只留下更暗的核心。
    这条规则让识别更严格，避免把一片中等暗度的大区域当成暗条。

    参数
    ----
    seed     : (Tx,Ty) 角秒。给定时取局部窗口内、以该点所属连通域为准的一块。
    thr      : 绝对暗化阈值。不给时按种子暗度的 ``seed_frac`` 倍（下限 0.08）。
    mean_frac: 区域平均暗度相对种子暗度的下限，默认 0.90。
    local_radius : 有种子时只在这个半径内找。

    返回 dict 或 None：``contour``/``spine``（HPC 角秒）、``area_px``、
    ``mean_dark``、``seed_dark``、``thr``、``mask``、``centroid``。
    """
    d = np.clip(np.nan_to_num(view.darkening(sigma), nan=0.0), 0.0, None)
    region = view.region_mask(win)
    if restrict is not None:
        region = region & np.asarray(restrict, bool)
    seed_ij = None
    seed_dark = None
    if seed is not None:
        px, py = view.to_pixel(seed[0], seed[1])
        i, j = int(round(float(px[0]))), int(round(float(py[0])))
        ny, nx = d.shape
        i, j = int(np.clip(i, 0, nx - 1)), int(np.clip(j, 0, ny - 1))
        seed_ij = (j, i)
        seed_dark = float(d[j, i])
        rr = max(local_radius / max(view.scale, 1e-6), 4.0)
        yy, xx = np.mgrid[0:ny, 0:nx]
        region = region & ((xx - i) ** 2 + (yy - j) ** 2 <= rr * rr)
        if thr is None:
            thr = max(seed_frac * seed_dark, 0.06)
    elif thr is None:
        pool = d[region]
        if pool.size < 64:
            return None
        thr = max(float(np.percentile(pool, percentile)), 0.06)

    best = None
    n_iter = max_iter if seed is not None else 1
    for _ in range(n_iter):
        sel = _grow(d, region, thr, seed_ij, min_area, bridge)
        if sel is None:
            break
        mean = float(np.nanmean(d[sel]))
        best = (sel, mean, int(sel.sum()), float(thr))
        if seed is None or seed_dark is None or mean >= mean_frac * seed_dark:
            break
        thr *= 1.12
    if best is None:
        return None
    sel, mean, area, thr = best
    if fill:
        sel = ndi.binary_fill_holes(sel)
    c = contour_from_mask(sel)
    if c is None or len(c) < 12:
        return None
    x, y = view.to_world(c[:, 0], c[:, 1])
    contour = np.column_stack([x, y])
    good = np.isfinite(contour).all(axis=1)
    contour = contour[good]
    spine = spine_from_mask(sel, view)
    return dict(contour=contour, spine=spine, area_px=area,
                mean_dark=mean, seed_dark=seed_dark, thr=thr, mask=sel,
                centroid=(float(contour[:, 0].mean()), float(contour[:, 1].mean())))


def corridor_points(m_src, m_dst, contour_hpc, frame, r_lo=1.0, r_hi=1.5,
                    n_r=5, max_pts=240):
    """把参考像上的轮廓按三维体投到另一视角，得到走廊点集。"""
    c = np.asarray(contour_hpc, float)
    if len(c) > max_pts:
        idx = np.linspace(0, len(c) - 1, max_pts).astype(int)
        c = c[idx]
    rs = np.linspace(r_lo, r_hi, n_r)
    X = np.repeat(c[:, 0], n_r)
    Y = np.repeat(c[:, 1], n_r)
    R = np.tile(rs, len(c))
    D = m_src.observer_coordinate.radius.to_value(u.km)
    v = G.los_point(X, Y, R, D)
    hc = Heliocentric(x=v[:, 0] * u.km, y=v[:, 1] * u.km, z=v[:, 2] * u.km,
                      observer=m_src.observer_coordinate, obstime=m_src.date)
    sky = SkyCoord(hc).transform_to(frame)
    q = sky.transform_to(m_dst.coordinate_frame)
    x = np.asarray(q.Tx.to_value(u.arcsec), float)
    y = np.asarray(q.Ty.to_value(u.arcsec), float)
    ok = np.isfinite(x) & np.isfinite(y)
    return np.column_stack([x[ok], y[ok]])


def corridor_mask(view, pts_hpc, grow=14):
    """把走廊点集光栅化成原生像素掩膜。"""
    ny, nx = view.m.data.shape
    pts = np.asarray(pts_hpc, float)
    px, py = view.to_pixel(pts[:, 0], pts[:, 1])
    m = np.zeros((ny, nx), bool)
    ii = np.round(px).astype(int)
    jj = np.round(py).astype(int)
    ok = (ii >= 0) & (ii < nx) & (jj >= 0) & (jj < ny)
    m[jj[ok], ii[ok]] = True
    if not m.any():
        return m
    m = ndi.binary_dilation(m, ndi.generate_binary_structure(2, 2), iterations=int(grow))
    m = ndi.binary_closing(m, ndi.generate_binary_structure(2, 2), iterations=3)
    return m


def darkest_in_corridor(view, pts_hpc, sigma=45.0):
    """走廊内暗度最大的点，返回 (x, y, dark)。"""
    d = np.clip(np.nan_to_num(view.darkening(sigma), nan=0.0), 0.0, None)
    ny, nx = d.shape
    px, py = view.to_pixel(pts_hpc[:, 0], pts_hpc[:, 1])
    ii = np.round(px).astype(int)
    jj = np.round(py).astype(int)
    ok = (ii >= 0) & (ii < nx) & (jj >= 0) & (jj < ny)
    if not ok.any():
        return None
    vals = d[jj[ok], ii[ok]]
    k = int(np.argmax(vals))
    ii, jj = ii[ok], jj[ok]
    x, y = view.to_world([ii[k]], [jj[k]])
    return float(x[0]), float(y[0]), float(vals[k])
