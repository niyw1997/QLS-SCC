"""四格排版与出图。

(a) 参考视角（AIA）  (b) 另一视角（EUVI）  (c) HEEQ 三维  (d) PCA 局部平面

沿用工作区的多子图约定：不画总标题；左上角 (a)(b)(c)(d)；子图内时刻用白底
小框放在右下角；同区域的 (a)(b) 用同一套极线与点色标。
"""

from __future__ import annotations

import numpy as np
import matplotlib
from matplotlib.collections import LineCollection
from matplotlib.colors import LogNorm, Normalize
from matplotlib.figure import Figure
from matplotlib.cm import ScalarMappable

from . import geometry as G

matplotlib.rcParams["font.sans-serif"] = [
    "PingFang SC", "Heiti SC", "Hiragino Sans GB", "Arial Unicode MS", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

__all__ = ["make_figure", "save"]

BG = "#0e0e0e"
FG = "#e8e8e8"
GRID = "#4a4a4a"
CYAN = "#00e5ff"
YELLOW = "#ffe14d"
MAGENTA = "#ff5cf0"

CMAP_R = matplotlib.colormaps["turbo"]
NORM_R = Normalize(vmin=1.0, vmax=1.6)


def _cmap(name):
    try:
        return matplotlib.colormaps[name]
    except Exception:
        return matplotlib.colormaps["inferno"]


def _letter(ax, ch):
    kw = dict(transform=ax.transAxes, ha="left", va="top", fontsize=13,
              color="white", zorder=20,
              bbox=dict(fc="black", ec="none", alpha=0.55, pad=2.2))
    if hasattr(ax, "text2D"):
        ax.text2D(0.018, 0.975, f"({ch})", **kw)
    else:
        ax.text(0.018, 0.975, f"({ch})", **kw)


def _stamp(ax, text):
    ax.text(0.975, 0.028, text, transform=ax.transAxes, ha="right", va="bottom",
            fontsize=7.6, color="black", zorder=20,
            bbox=dict(fc="white", ec="none", alpha=0.9, pad=2.0))


def _tag(ax, text):
    ax.text(0.025, 0.045, text, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.6, color=FG, zorder=20,
            bbox=dict(fc="black", ec="none", alpha=0.5, pad=2.2))


def _style_2d(ax):
    ax.set_facecolor(BG)
    for s in ax.spines.values():
        s.set_color("#666666")
        s.set_linewidth(0.7)
    ax.tick_params(colors="#bbbbbb", labelsize=7.5, direction="in",
                   top=True, right=True, width=0.6, length=3)
    ax.xaxis.label.set_color("#cccccc")
    ax.yaxis.label.set_color("#cccccc")


def _panel_image(ax, view, win, levels=None, max_px=760):
    sub, ext = view.subimage(win, max_px=max_px)
    lo, hi = levels or view.levels(win)
    if view.norm_kind == "log":
        lo = max(lo, 1e-6)
        hi = max(hi, lo * 1.0001)
        norm = LogNorm(vmin=lo, vmax=hi)
    elif view.norm_kind == "power":
        from matplotlib.colors import PowerNorm
        norm = PowerNorm(gamma=getattr(view, "gamma", 0.4), vmin=lo,
                         vmax=max(hi, lo + 1e-6))
    else:
        norm = Normalize(vmin=lo, vmax=max(hi, lo + 1e-6))
    cmap = _cmap(view.cmap).copy()
    cmap.set_bad("black")
    ax.imshow(sub, origin="lower", extent=ext, cmap=cmap, norm=norm,
              interpolation="nearest", rasterized=True, zorder=1)
    ax.set_xlim(win.x0, win.x1)
    ax.set_ylim(win.y0, win.y1)
    _style_2d(ax)
    ax.set_aspect("equal")


def _scale_bar(ax, win, length=200.0):
    x0, x1, y0, y1 = win.as_tuple()
    xa = x1 - 0.06 * (x1 - x0) - length
    ya = y1 - 0.06 * (y1 - y0)
    ax.plot([xa, xa + length], [ya, ya], color="white", lw=2.0, zorder=15,
            solid_capstyle="butt")
    ax.text(xa + 0.5 * length, ya + 0.02 * (y1 - y0), f"{length:.0f}\"",
            color="white", fontsize=7.5, ha="center", va="bottom", zorder=15)


def _draw_points(ax, pts, which, count_r=True):
    for k, p in enumerate(pts):
        q = p.pa if which == "a" else p.pb
        col = CMAP_R(NORM_R(p.r))
        ax.plot([q[0]], [q[1]], "o", ms=7.2, mfc=col, mec="white", mew=0.8,
                zorder=12)
        ax.text(q[0] + 10, q[1] + 8, str(p.i), color="white", fontsize=6.5,
                zorder=13, clip_on=True,
                bbox=dict(fc="black", ec="none", alpha=0.45, pad=1.2))


def _epipolar(ax, epi):
    x, y, r = epi["x"], epi["y"], epi["r"]
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 2:
        return
    seg = np.stack([np.column_stack([x, y])[:-1], np.column_stack([x, y])[1:]], axis=1)
    lc = LineCollection(seg, cmap=CMAP_R, norm=NORM_R, linewidths=2.2, zorder=9)
    lc.set_array(0.5 * (r[:-1] + r[1:]))
    ax.add_collection(lc)
    for rv in (1.0, 1.2, 1.4, 1.6, 1.8, 2.0):
        k = int(np.argmin(abs(r - rv)))
        if not (np.isfinite(x[k]) and np.isfinite(y[k])):
            continue
        ax.plot([x[k]], [y[k]], "|", ms=6, color="white", mew=1.0, zorder=10)
        ax.text(x[k] + 8, y[k], f"{rv:.1f}", color="white", fontsize=6.0,
                va="center", zorder=10, clip_on=True,
                bbox=dict(fc="black", ec="none", alpha=0.4, pad=1.0))


CMAP_H = matplotlib.colormaps["viridis"]


def _norm_h(H):
    lo, hi = float(np.min(H)), float(np.max(H))
    if hi - lo < 0.5:            # 有时两点几乎同高，给色标留 0.5 Mm 的最小跨度
        mid = 0.5 * (lo + hi)
        lo, hi = mid - 0.25, mid + 0.25
    return Normalize(vmin=lo, vmax=hi)


def _draw_height3d(ax, session):
    """(c) 经度、纬度、高度 H 的三维直角坐标，底面固定在 H=0。"""
    ax.set_facecolor(BG)
    if not session.points:
        lon0, lon1, lat0, lat1 = -1.0, 1.0, -1.0, 1.0
        for x in np.linspace(lon0, lon1, 5):
            ax.plot([x, x], [lat0, lat1], [0, 0], color="#2f2f2f", lw=0.5)
        for y in np.linspace(lat0, lat1, 5):
            ax.plot([lon0, lon1], [y, y], [0, 0], color="#2f2f2f", lw=0.5)
        ax.text2D(0.5, 0.90, "至少 2 个点才能画三维结构",
                  transform=ax.transAxes, ha="center", va="center",
                  color="#b0b0b0", fontsize=9, zorder=20,
                  bbox=dict(fc=BG, ec="#333333", alpha=0.75, pad=4))
        ax.set_xlim(lon0, lon1)
        ax.set_ylim(lat0, lat1)
        ax.set_zlim(0, 1)
        ax.set_xlabel("经度 (deg)", fontsize=7.2, labelpad=2, color="white")
        ax.set_ylabel("纬度 (deg)", fontsize=7.2, labelpad=2, color="white")
        ax.set_zlabel("H (Mm)", fontsize=7.2, labelpad=2, color="white")
        ax.set_box_aspect((1, 1, 0.7), zoom=0.85)
        ax.view_init(elev=24, azim=-60)
        ax.tick_params(labelsize=6.2, colors="#dddddd", pad=-2)
        ax.locator_params(nbins=4)
        for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
            axis.pane.fill = False
            axis.pane.set_edgecolor("#333333")
        ax.grid(color="#2e2e2e", lw=0.4)
        return

    lon = np.array([p.res["lon"] for p in session.points], float)
    lat = np.array([p.res["lat"] for p in session.points], float)
    H = np.array([p.res["height_mm"] for p in session.points], float)
    n = len(H)
    dlon = max(lon.max() - lon.min(), 0.4)
    dlat = max(lat.max() - lat.min(), 0.4)
    lon0, lon1 = lon.min() - 0.18 * dlon, lon.max() + 0.18 * dlon
    lat0, lat1 = lat.min() - 0.18 * dlat, lat.max() + 0.18 * dlat
    h1 = max(float(H.max()) * 1.18, 10.0)

    # H=0 的底平面
    LX, LY = np.meshgrid([lon0, lon1], [lat0, lat1])
    ax.plot_surface(LX, LY, np.zeros_like(LX), color="#555555", alpha=0.16,
                    linewidth=0, zorder=0)
    for x in np.linspace(lon0, lon1, 5):
        ax.plot([x, x], [lat0, lat1], [0, 0], color="#3a3a3a", lw=0.5, zorder=0)
    for y in np.linspace(lat0, lat1, 5):
        ax.plot([lon0, lon1], [y, y], [0, 0], color="#3a3a3a", lw=0.5, zorder=0)

    # 每点到 H=0 的垂线
    for i in range(n):
        ax.plot([lon[i], lon[i]], [lat[i], lat[i]], [0, H[i]], ls=":",
                color="#777777", lw=0.6, zorder=1)

    # 平滑三维曲线
    if n >= 2:
        from scipy.interpolate import CubicSpline
        t = np.linspace(0.0, 1.0, n)
        tt = np.linspace(0.0, 1.0, 240)
        if n == 2:
            X, Y, Z = (np.interp(tt, t, lon), np.interp(tt, t, lat),
                       np.interp(tt, t, H))
        else:
            X = CubicSpline(t, lon, bc_type="natural")(tt)
            Y = CubicSpline(t, lat, bc_type="natural")(tt)
            Z = CubicSpline(t, H, bc_type="natural")(tt)
        ax.plot(X, Y, Z, color=YELLOW, lw=1.4, zorder=3)

    sc = ax.scatter(lon, lat, H, c=H, cmap=CMAP_H, norm=_norm_h(H), s=42,
                    depthshade=False, edgecolors="white", linewidths=0.6, zorder=5)

    ax.set_xlim(lon0, lon1)
    ax.set_ylim(lat0, lat1)
    ax.set_zlim(0, h1)
    ax.set_box_aspect((1, 1, 0.85), zoom=1.0)
    ax.view_init(elev=24, azim=-60)
    ax.set_xlabel("经度 (deg)", fontsize=7.2, labelpad=2, color="white")
    ax.set_ylabel("纬度 (deg)", fontsize=7.2, labelpad=2, color="white")
    ax.set_zlabel("H (Mm)", fontsize=7.2, labelpad=2, color="white")
    ax.tick_params(labelsize=6.2, colors="#dddddd", pad=-2)
    ax.locator_params(nbins=4)
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.fill = False
        axis.pane.set_edgecolor("#333333")
    ax.grid(color="#2e2e2e", lw=0.4)

    cb = ax.figure.colorbar(sc, ax=ax, shrink=0.55, pad=0.10, aspect=24)
    cb.set_label("H (Mm)", color="white", fontsize=7.2)
    cb.ax.tick_params(labelsize=6.2, colors="#dddddd")
    cb.ax.yaxis.get_major_formatter().set_useOffset(False)
    cb.outline.set_edgecolor("#555555")


def _draw_height_profile(ax, session):
    """(d) 沿 AIA 像面 xy 弧长的高度 H 分布。"""
    ax.set_facecolor(BG)
    _style_2d(ax)
    ax.set_xlabel("沿 AIA 像面的弧长 (arcsec)", fontsize=7.5, color="#cccccc")
    ax.set_ylabel("H (Mm)", fontsize=7.5, color="#cccccc")
    if len(session.points) < 2:
        ax.text(0.5, 0.5, "至少 2 个点才能画 H 分布",
                transform=ax.transAxes, ha="center", va="center",
                color="#8a8a8a", fontsize=9.5)
        ax.set_xticks([])
        ax.set_yticks([])
        return

    s = session.arc_length_a()
    H = np.array([p.res["height_mm"] for p in session.points], float)
    if len(H) >= 3:
        from scipy.interpolate import CubicSpline
        ss = np.linspace(s[0], s[-1], 300)
        HH = CubicSpline(s, H, bc_type="natural")(ss)
    else:
        ss, HH = s, H
    ax.plot(ss, HH, color=YELLOW, lw=1.4, zorder=3)
    ax.scatter(s, H, c=H, cmap=_cmap("viridis"), norm=_norm_h(H),
               s=40, edgecolors="white", linewidths=0.6, zorder=5)
    pad_x = 0.04 * max(s[-1] - s[0], 1.0)
    ax.set_xlim(s[0] - pad_x, s[-1] + pad_x)
    mid = 0.5 * (float(H.min()) + float(H.max()))
    span = max(float(H.max()) - float(H.min()), 0.5)
    lo, hi = mid - span / 2, mid + span / 2
    pad_y = 0.12 * span
    ax.set_ylim(lo - pad_y, hi + pad_y)
    ax.yaxis.get_major_formatter().set_useOffset(False)
    if lo > 0:
        ax.axhline(0.0, color="#555555", lw=0.6, ls=":")   # 在视野外，仅作提示
    ax.locator_params(nbins=5)
    txt = (f"弧长 {s[-1]:.0f}″\n"
           f"H {H.min():.1f} 到 {H.max():.1f} Mm\n"
           f"中位 {np.median(H):.1f} Mm")
    ax.text(0.03, 0.90, txt, transform=ax.transAxes, va="top", ha="left",
            color="#dddddd", fontsize=7.2,
            bbox=dict(fc="#181818", ec="#333333", pad=4))


def make_figure(session, size=(12.4, 9.0), dpi=150, levels=None, fig=None):
    """画四格总图。GUI 与保存共用同一套代码。

    传入 ``fig`` 时在它上面重画（Tk 的画布要复用同一个 Figure 对象）。
    """
    views = session.views
    win = session.win
    if fig is None:
        fig = Figure(figsize=size, dpi=dpi, facecolor=BG)
    else:
        fig.clf()
        fig.set_facecolor(BG)
    gs = fig.add_gridspec(2, 2, left=0.055, right=0.985, bottom=0.10, top=0.97,
                          wspace=0.09, hspace=0.16)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0], projection="3d")
    # 没有点时 3D 里是一个整球，摆成正方形免得被格子裁掉；
    # 有点时视野会围住点云，保持整格更饱满。
    b = ax_c.get_position()
    fw, fh = fig.get_figwidth(), fig.get_figheight()
    if not session.points:
        side = min(b.width * fw, b.height * fh) * 0.90
        ax_c.set_position([b.x0 + (b.width - side / fw) / 2,
                           b.y0 + (b.height - side / fh) / 2 + 0.015,
                           side / fw, side / fh])

    ax_d = fig.add_subplot(gs[1, 1])
    fig._sccm_axes = dict(a=ax_a, b=ax_b, c=ax_c, d=ax_d)

    # (a)(b) 图像
    for ax, key, ch in ((ax_a, "a", "a"), (ax_b, "b", "b")):
        v = views[key]
        lev = (levels or {}).get(key)
        _panel_image(ax, v, win[key], levels=lev)
        _letter(ax, ch)
        _tag(ax, f"{v.label}   {v.obstime.strftime('%H:%M:%S')} UT")
        _stamp(ax, v.obstime.strftime("%H:%M:%S UT"))
        _scale_bar(ax, win[key], 200.0)
        ax.set_xlabel("HPC x (arcsec)", fontsize=7.5)
        ax.set_ylabel("HPC y (arcsec)", fontsize=7.5)

    if len(session.points) >= 2:
        cv = session.curve_B(300)
        if cv is not None:
            ax_b.plot(cv["x"], cv["y"], "-", color=YELLOW, lw=1.1, alpha=0.8,
                      zorder=4)
        ca = session.curve_A(300)
        if ca is not None:
            ax_a.plot(ca["x"], ca["y"], "-", color=CYAN, lw=1.0, alpha=0.75,
                      zorder=4)

    # 待确认的起点与极线，起点可以在任一格
    if session.pending is not None:
        x, y = session.pending["x"], session.pending["y"]
        src = session.pending.get("src", "a")
        ax_src = ax_a if src == "a" else ax_b
        ax_dst = ax_b if src == "a" else ax_a
        ax_src.plot([x], [y], "+", ms=16, color=MAGENTA, mew=1.8, zorder=14)
        ax_src.plot([x], [y], "o", ms=18, mfc="none", mec=MAGENTA, mew=1.1, zorder=14)
        _epipolar(ax_dst, session.pending["epi"])
        if session.guide_idx is not None:
            ax_src.text(x + 6, y + 6, str(session.guide_idx + 1),
                        color=MAGENTA, fontsize=8, zorder=15)

    # 预置点已取过的对映点
    if session.guides and session.guide_idx is not None:
        g = session.guides[session.guide_idx]
        dst = "b" if session.guide_src == "a" else "a"
        ax_dst = ax_b if dst == "b" else ax_a
        if g.dst is not None:
            ax_dst.plot([g.dst[0]], [g.dst[1]], "o", ms=17, mfc="none",
                        mec="#7cffb2", mew=1.6, zorder=14)
            ax_dst.annotate(f"点 {session.guide_idx + 1}",
                            xy=g.dst, xytext=(7, 7), textcoords="offset points",
                            color="#7cffb2", fontsize=7.5, zorder=15)

    # 已测点
    _draw_points(ax_a, session.points, "a")
    _draw_points(ax_b, session.points, "b")

    # (c)(d)
    _draw_height3d(ax_c, session)
    _letter(ax_c, "c")
    _draw_height_profile(ax_d, session)
    _letter(ax_d, "d")

    # 底部小结（不抢子图版面）
    s = session.summary()
    if s.get("n"):
        txt = (f"n = {s['n']}    r = {s['r_median']:.4f} R_sun    "
               f"height = {s['height_median_mm']:.0f} Mm    "
               f"median ray gap = {s['gap_median_km']:.0f} km")
        fig.text(0.055, 0.012, txt, color="#cfcfcf", fontsize=8.6, ha="left")
    return fig


def save(session, outdir=None, dpi=300, size=(12.4, 9.0)):
    outdir = outdir or session.outdir
    return session.save(make_figure=True, dpi=dpi, size=size)
