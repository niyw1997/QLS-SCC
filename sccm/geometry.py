"""双视角几何内核。

四件事
------
1. 视线：从观测者出发，用真三维直线（DSUN_OBS + 1001 式解析求交）；
2. 极线：把一条视线上的三维点逐个转到另一观测者，再投成 HPC 角秒；
3. 三角化：两条视线的最近点取中点，给出 HEEQ 坐标与两视线的失配量；
4. PCA：点云拟合椭圆环，给长短半径、环面法向与拟合残差。

坐标一律交给 SunPy。所有 `Heliocentric(observer=...)` 的点先转到
`HeliographicStonyhurst`（HEEQ 的定义：z 沿自转轴，x 指向地球中心子午面与
赤道面的交线），之后的 lon/lat 就是 HEEQ 经纬度。
"""

from __future__ import annotations

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import Heliocentric, HeliographicStonyhurst

__all__ = [
    "RSUN_KM", "RSUN_MM", "ref_frame", "observer_origin", "los_point",
    "common_point", "project_to", "epipolar", "triangulate", "pca_fit",
    "link_box", "separation_deg", "height_mm",
]

RSUN_KM = 695700.0          # IAU 2015 名义半径，km
RSUN_MM = 696.0             # 1 R_sun 折合多少 Mm，只用于显示高度


def ref_frame(obstime):
    """公共坐标架。HEEQ 的 lon 以地球中心子午面为零点。"""
    return HeliographicStonyhurst(obstime=obstime)


def _xyz(sky):
    c = sky.cartesian
    return np.array([c.x.to_value(u.km), c.y.to_value(u.km),
                     c.z.to_value(u.km)], float)


def height_mm(r):
    """r（R_sun）折合到日面之上的高度，单位 Mm。"""
    return (np.asarray(r, float) - 1.0) * RSUN_MM


# --------------------------------------------------------------------------- #
# 视线
# --------------------------------------------------------------------------- #
def los_point(tx_arcsec, ty_arcsec, r, D_km):
    """观测者本地系里的视线点，单位 km。

    观测者在 (0, 0, D)，太阳在原点；HPC (Tx, Ty) 给出方向
    ``d = (tan Tx, tan Ty, -1) / |...|``，沿 d 走到日心距 ``r * R_sun`` 的近根。
    视线点就是 ``D*(0,0,1) + s*d``。

    参数可以是标量也可以是数组（广播），返回形状 ``(..., 3)``。
    """
    tx = np.atleast_1d(np.asarray(tx_arcsec, float))
    ty = np.atleast_1d(np.asarray(ty_arcsec, float))
    rs = np.atleast_1d(np.asarray(r, float))
    tx, ty, rs = np.broadcast_arrays(tx, ty, rs)
    Tx = np.radians(tx / 3600.0)
    Ty = np.radians(ty / 3600.0)
    d = np.stack([np.tan(Tx), np.tan(Ty), -np.ones_like(Tx)], axis=-1)
    d /= np.linalg.norm(d, axis=-1, keepdims=True)
    D = float(D_km)
    R = rs * RSUN_KM
    b = 2.0 * D * d[..., 2]
    c = D * D - R * R
    disc = b * b - 4.0 * c
    disc = np.maximum(disc, 0.0)
    s = (-b - np.sqrt(disc)) / 2.0
    out = d * s[..., None]
    out[..., 2] += D
    return out


def observer_origin(m, frame):
    """观测者在公共坐标架里的位置，km。"""
    D = m.observer_coordinate.radius.to_value(u.km)
    hc = Heliocentric(x=0.0 * u.km, y=0.0 * u.km, z=D * u.km,
                      observer=m.observer_coordinate, obstime=m.date)
    return _xyz(SkyCoord(hc).transform_to(frame))


def common_point(m, tx_arcsec, ty_arcsec, r, frame):
    """把观测者本地系里的视线点转到公共坐标架。"""
    D = m.observer_coordinate.radius.to_value(u.km)
    v = los_point(tx_arcsec, ty_arcsec, r, D)
    v = np.atleast_2d(v)
    hc = Heliocentric(x=v[:, 0] * u.km, y=v[:, 1] * u.km, z=v[:, 2] * u.km,
                      observer=m.observer_coordinate, obstime=m.date)
    return SkyCoord(hc).transform_to(frame)


def project_to(m, sky):
    """公共坐标架的点投到 ``m`` 的 HPC，返回 (Tx, Ty) 角秒。

    投不上（在观测者背后、坐标变换失败）返回 ``(nan, nan)``。
    """
    try:
        q = SkyCoord(sky).transform_to(m.coordinate_frame)
        return (np.asarray(q.Tx.to_value(u.arcsec), float),
                np.asarray(q.Ty.to_value(u.arcsec), float))
    except Exception:
        n = np.shape(np.atleast_1d(sky))[0]
        return np.full(n, np.nan), np.full(n, np.nan)


def epipolar(m_a, m_b, tx_arcsec, ty_arcsec, frame, r_lo=1.0, r_hi=2.0, n=161):
    """``m_a`` 上一点 (Tx, Ty) 的视线投到 ``m_b``，按 r 采样。

    返回 dict：

    =========  ==================================================
    ``x, y``   ``m_b`` 的 HPC 角秒
    ``r``      对应的日心距（R_sun）
    ``sph``    公共架里的三维单位方向（供画图）
    =========  ==================================================
    """
    rs = np.linspace(r_lo, r_hi, int(n))
    D = m_a.observer_coordinate.radius.to_value(u.km)
    v = los_point(np.full_like(rs, tx_arcsec), np.full_like(rs, ty_arcsec), rs, D)
    hc = Heliocentric(x=v[:, 0] * u.km, y=v[:, 1] * u.km, z=v[:, 2] * u.km,
                      observer=m_a.observer_coordinate, obstime=m_a.date)
    sky = SkyCoord(hc).transform_to(frame)
    q = SkyCoord(sky).transform_to(m_b.coordinate_frame)
    x = np.asarray(q.Tx.to_value(u.arcsec), float)
    y = np.asarray(q.Ty.to_value(u.arcsec), float)
    return dict(x=x, y=y, r=rs, sph=sky)


# --------------------------------------------------------------------------- #
# 三角化
# --------------------------------------------------------------------------- #
def _ray(m, tx, ty, frame):
    """一条视线的 (起点, 单位方向)，都在公共坐标架里，km。"""
    o = observer_origin(m, frame)
    p = common_point(m, tx, ty, 1.5, frame)
    p = _xyz(p)[:, 0]
    d = p - o
    d /= np.linalg.norm(d)
    return o, d


def triangulate(m_a, m_b, pa, pb, frame):
    """两条视线的最近点。

    ``pa``/``pb`` 是两幅图上的 (Tx, Ty) 角秒。返回 dict：

    ``xyz`` 最短连线的中点（km，公共架），``r`` 日心距（R_sun），
    ``lon``/``lat`` HEEQ 经纬度（deg），``height_mm`` 高度，
    ``gap_km`` 两条视线在最近处的距离（判别配点是否自洽）。
    """
    o1, u1 = _ray(m_a, pa[0], pa[1], frame)
    o2, u2 = _ray(m_b, pb[0], pb[1], frame)
    w0 = o1 - o2
    b = float(u1 @ u2)
    d = float(u1 @ w0)
    e = float(u2 @ w0)
    den = 1.0 - b * b
    if abs(den) < 1e-12:                       # 两条视线平行，几何退化
        return dict(xyz=np.full(3, np.nan), r=np.nan, lon=np.nan, lat=np.nan,
                    height_mm=np.nan, gap_km=np.nan, degenerate=True)
    s = (b * e - d) / den
    t = b * s + e
    p1 = o1 + s * u1
    p2 = o2 + t * u2
    p = 0.5 * (p1 + p2)
    # frame 是 HGS（z 沿自转轴、x 指向地球中心子午面），它的直角坐标就是 HEEQ
    r = float(np.linalg.norm(p) / RSUN_KM)
    x, y, z = p
    lon = float(np.degrees(np.arctan2(y, x)))
    lat = float(np.degrees(np.arcsin(np.clip(z / (r * RSUN_KM), -1, 1))))
    return dict(xyz=p, xyz_rs=p / RSUN_KM, r=r, lon=lon, lat=lat,
                height_mm=float(height_mm(r)),
                gap_km=float(np.linalg.norm(p1 - p2)), degenerate=False)


# --------------------------------------------------------------------------- #
# PCA
# --------------------------------------------------------------------------- #
def pca_fit(points_km):
    """把一串有序三维点拟合成椭圆环。

    ``points_km`` 形状 (N, 3)，顺序必须从一脚到另一脚（不能闭合）。
    中点取首末两点的中点，与论文 ``r_M = (p_1 + p_N)/2`` 一致。

    返回 dict：``midpoint``、``evalue``（升序）、``evector``（列 = n, a, b）、
    ``ra``/``rb``（半轴，km）、``rn``（法向厚度）、``residual``（平面残差 km）、
    ``sigma_ra``/``sigma_rb``（半轴不确定度 km）、``local``（N×3 局部坐标）。
    """
    p = np.asarray(points_km, float)
    n = len(p)
    if n < 3:
        return None
    mid = 0.5 * (p[0] + p[-1])
    d = p - mid
    cov = (d.T @ d) / n
    ev, evec = np.linalg.eigh(cov)             # 升序
    order = np.argsort(ev)
    ev, evec = ev[order], evec[:, order]
    if np.linalg.det(evec) < 0:                # 右手 triad
        evec[:, [0, 2]] = -evec[:, [0, 2]]
    ra = float(np.sqrt(2.0 * max(ev[1], 0.0)))
    rb = float(np.sqrt(2.0 * max(ev[2], 0.0)))
    rn = float(np.sqrt(2.0 * max(ev[0], 0.0)))
    local = d @ evec                          # (N,3)，列依次为 n, a, b
    # 椭圆半轴的不确定度：先归一化相位，再算残差均方根
    x = local[:, 1]
    y = local[:, 2]
    rr = np.hypot(x, y)
    ok = rr > 0
    cos_t = np.where(ok, x / np.where(ok, rr, 1), 0.0)
    sin_t = np.where(ok, y / np.where(ok, rr, 1), 0.0)
    sig_a = float(np.sqrt(np.mean((x - ra * cos_t) ** 2)))
    sig_b = float(np.sqrt(np.mean((y - rb * sin_t) ** 2)))
    resid = float(np.sqrt(np.mean(local[:, 0] ** 2)))
    return dict(midpoint=mid, evalue=ev, evector=evec, ra=ra, rb=rb, rn=rn,
                residual=resid, sigma_ra=sig_a, sigma_rb=sig_b,
                local=local, cos_t=cos_t, sin_t=sin_t)


def ellipse_curve(fit, n=361):
    """按 (ra, rb) 画一条闭合椭圆，返回到公共架的三维点 (n,3)。"""
    th = np.linspace(0.0, 2.0 * np.pi, n)
    loc = np.column_stack([np.zeros(n), fit["ra"] * np.cos(th),
                           fit["rb"] * np.sin(th)])
    return fit["midpoint"] + loc @ fit["evector"].T


# --------------------------------------------------------------------------- #
# 视场联动
# --------------------------------------------------------------------------- #
def link_box(m_src, m_dst, box, frame, r_lo=1.0, r_hi=1.5,
             n_edge=36, n_r=5, pad=20.0):
    """把 ``m_src`` 的一块视场当作三维体，投到 ``m_dst`` 上的包围盒。

    视场矩形的四条边各取 ``n_edge`` 点，每条视线在 ``r_lo`` 到 ``r_hi``
    之间取 ``n_r`` 个半径，全部转到 ``m_dst`` 的 HPC，取包围盒外扩 ``pad``
    角秒。返回 ``(x0, x1, y0, y1)``，投不上返回 None。

    这个方法就是 pcp3 的 Link FoV：视场不再按像素平移，而是按三维体投影。
    """
    x0, x1, y0, y1 = box
    gx = np.linspace(x0, x1, n_edge)
    gy = np.linspace(y0, y1, n_edge)
    ex = np.r_[gx, gx, np.full(n_edge, x0), np.full(n_edge, x1)]
    ey = np.r_[np.full(n_edge, y0), np.full(n_edge, y1), gy, gy]
    rs = np.linspace(r_lo, r_hi, n_r)
    X = np.repeat(ex, n_r)
    Y = np.repeat(ey, n_r)
    R = np.tile(rs, len(ex))
    D = m_src.observer_coordinate.radius.to_value(u.km)
    v = los_point(X, Y, R, D)
    hc = Heliocentric(x=v[:, 0] * u.km, y=v[:, 1] * u.km, z=v[:, 2] * u.km,
                      observer=m_src.observer_coordinate, obstime=m_src.date)
    sky = SkyCoord(hc).transform_to(frame)
    q = SkyCoord(sky).transform_to(m_dst.coordinate_frame)
    tx = np.asarray(q.Tx.to_value(u.arcsec), float)
    ty = np.asarray(q.Ty.to_value(u.arcsec), float)
    ok = np.isfinite(tx) & np.isfinite(ty)
    if ok.sum() < 8:
        return None
    return (float(tx[ok].min() - pad), float(tx[ok].max() + pad),
            float(ty[ok].min() - pad), float(ty[ok].max() + pad))


def separation_deg(m_a, m_b, frame):
    """两观测者的日心经度差（度）。"""
    a = m_a.observer_coordinate.transform_to(frame)
    b = m_b.observer_coordinate.transform_to(frame)
    return float(abs(a.lon.deg - b.lon.deg))


def epipolar_sensitivity(m_a, m_b, frame, tx=0.0, ty=0.0,
                         r0=1.1, dr=0.1):
    """沿极线每 ``dr`` R_sun 在 ``m_b`` 上移动多少角秒。"""
    e1 = epipolar(m_a, m_b, tx, ty, frame, r_lo=r0, r_hi=r0 + dr, n=2)
    if not np.isfinite(e1["x"]).all():
        return float("nan")
    return float(np.hypot(e1["x"][1] - e1["x"][0], e1["y"][1] - e1["y"][0]))
