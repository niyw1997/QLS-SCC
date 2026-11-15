"""测量会话：把界面、演示与命令行共用的状态与动作放在一起。

GUI 只负责画和接鼠标，真正的"发生了什么"都在这里：
视场联动、取点、三角化、PCA、保存。
"""

from __future__ import annotations

import csv
import json
import os
from dataclasses import dataclass

import numpy as np

from . import __version__, geometry as G
from .views import View, Window, load_view

__all__ = ["Session", "Point", "Guide"]

LINK_R_HI_DEFAULT = 1.10      # 视场联动用 1.0-1.1 的薄壳来定中心（pcp3 的做法）
EPI_R = (1.0, 2.0)            # 极线画到 2 R_sun，便于看高度刻度


@dataclass
class Point:
    i: int
    pa: tuple            # (Tx,Ty) 角秒
    pb: tuple
    res: dict
    seq: int | None = None      # 预置点列表里的序号，自由取点为 None
    note: str = ""

    @property
    def r(self):
        return self.res["r"]

    @property
    def xyz(self):
        return self.res["xyz_rs"]

    def as_dict(self):
        return dict(i=self.i, pa=list(map(float, self.pa)), pb=list(map(float, self.pb)),
                    heeq=[float(v) for v in self.xyz],
                    r=float(self.res["r"]), lon=float(self.res["lon"]),
                    lat=float(self.res["lat"]), height_mm=float(self.res["height_mm"]),
                    gap_km=float(self.res["gap_km"]),
                    guide=None if self.seq is None else self.seq + 1,
                    note=self.note)


@dataclass
class Guide:
    i: int
    src: tuple            # 预置点在源图上的 (x,y) 角秒
    dst: tuple | None = None
    point: Point | None = None


class Session:
    def __init__(self, obs1, obs2, outdir=".", r_link_hi=LINK_R_HI_DEFAULT,
                 verbose=True):
        self.a = load_view(obs1)
        self.b = load_view(obs2)
        self.views = {"a": self.a, "b": self.b}
        self.frame = G.ref_frame(max(self.a.obstime, self.b.obstime))
        self.outdir = os.path.abspath(outdir or ".")
        os.makedirs(self.outdir, exist_ok=True)
        self.r_link_hi = float(r_link_hi)
        self.win = {"a": self.a.full_window(), "b": self.b.full_window()}
        self.points: list[Point] = []
        self.pending = None          # dict(src,x,y,epi) 在某一格起过点、等另一格确认
        self.guide_src = None        # "a" 或 "b"，预置点在哪一格
        self.guides: list[Guide] = []
        self.guide_idx = None
        self.message = "在任一格点第一个点，另一格会画极线。"
        self.separation = G.separation_deg(self.a.m, self.b.m, self.frame)
        self.sensitivity = G.epipolar_sensitivity(self.a.m, self.b.m, self.frame)
        if verbose:
            dt = abs((self.a.obstime - self.b.obstime).to_value("s"))
            print(f"[sccm] {self.a.short} {self.a.obstime} + {self.b.short} {self.b.obstime}")
            print(f"[sccm] 视角差 {self.separation:.2f}°   时刻差 {dt:.0f} s   "
                  f"极线灵敏度 {self.sensitivity:.1f}\"/0.1 R_sun")
            if self.separation < 5.0:
                print("[sccm] 警告：视角差小于 5°，高度几乎没有分辨力。")

    # ------------------------------------------------------------------ #
    # 视场
    # ------------------------------------------------------------------ #
    def set_window(self, key, win, link=True):
        self.win[key] = win
        other = "b" if key == "a" else "a"
        if not link:
            return
        src, dst = self.views[key], self.views[other]
        box = G.link_box(src.m, dst.m, win.as_tuple(), self.frame,
                         r_lo=1.0, r_hi=self.r_link_hi)
        if box is None:
            return
        cx = 0.5 * (box[0] + box[1])
        cy = 0.5 * (box[2] + box[3])
        # 用与源相同的视场边长，包围盒只用来定中心（pcp3 的 center_r_hi 做法）
        self.win[other] = Window(cx, cy, win.half).clipped(dst.extent())

    def zoom(self, key, box, link=True):
        w = Window.from_box(box).clipped(self.views[key].extent())
        self.set_window(key, w, link=link)
        return w

    def pan(self, key, dx, dy, link=True):
        w = Window(self.win[key].cx + dx, self.win[key].cy + dy, self.win[key].half)
        w = w.clipped(self.views[key].extent())
        self.set_window(key, w, link=link)
        return w

    def reset(self):
        self.win = {"a": self.a.full_window(), "b": self.b.full_window()}
        self.pending = None
        self.message = "已回到全盘。"

    # ------------------------------------------------------------------ #
    # 预置点序列
    # ------------------------------------------------------------------ #
    def load_guides(self, panel, points, windows=None):
        """载入一组预置源点，并停在第 1 个点。

        panel 是源图，"a" 或 "b"。points 是源图上的 (x,y) 角秒。
        windows 可选，给 {"a": (x0,x1,y0,y1), "b": (...)} 或 aia/euvi 键。
        """
        panel = str(panel).lower()
        if panel not in ("a", "b"):
            raise ValueError("panel 只能是 a 或 b")
        pts = []
        for q in points or []:
            if isinstance(q, dict):
                x = q.get("x", q.get("aia_x_arcsec", q.get("euvi_x_arcsec")))
                y = q.get("y", q.get("aia_y_arcsec", q.get("euvi_y_arcsec")))
            else:
                x, y = q[0], q[1]
            if x is None or y is None:
                continue
            pts.append((float(x), float(y)))
        if not pts:
            raise ValueError("预置点列表是空的")

        if windows:
            wa = windows.get("a", windows.get("aia"))
            wb = windows.get("b", windows.get("euvi"))
            if wa:
                self.win["a"] = Window.from_box(wa)
            if wb:
                self.win["b"] = Window.from_box(wb)

        self.guide_src = panel
        self.guide_idx = 0
        self.guides = [Guide(i=i, src=p) for i, p in enumerate(pts)]
        self.points = []
        self.pending = None
        self.goto_guide(0)
        return self.guides

    def goto_guide(self, i):
        if not self.guides:
            self.message = "还没有载入预置点列表。"
            return None
        i = max(0, min(int(i), len(self.guides) - 1))
        self.guide_idx = i
        g = self.guides[i]
        self.start_point(self.guide_src, g.src[0], g.src[1], guide=i)
        src_name = self.a.short if self.guide_src == "a" else self.b.short
        dst = "b" if self.guide_src == "a" else "a"
        dst_name = self.b.short if dst == "b" else self.a.short
        self.message = (f"引导点 {i + 1}/{len(self.guides)}，{src_name} 上已标出，"
                        f"到 {dst_name} 的极线上点对映点。")
        if g.point is not None:
            self.message += " 已取过，可再点一次覆盖。"
        return g

    def prev_guide(self):
        if not self.guides:
            self.message = "还没有载入预置点列表。"
            return None
        if self.guide_idx is None or self.guide_idx <= 0:
            self.message = "已经是第一个预置点。"
            return None
        return self.goto_guide(self.guide_idx - 1)

    def next_guide(self):
        if not self.guides:
            self.message = "还没有载入预置点列表。"
            return None
        if self.guide_idx is None or self.guide_idx >= len(self.guides) - 1:
            self.message = "已经是最后一个预置点。"
            return None
        return self.goto_guide(self.guide_idx + 1)

    def _sort_points_by_seq(self):
        seq = [p for p in self.points if p.seq is not None]
        rest = [p for p in self.points if p.seq is None]
        seq.sort(key=lambda p: p.seq)
        self.points = seq + rest
        for j, p in enumerate(self.points, start=1):
            p.i = j

    def curve_of(self, which="b", n=400):
        """把已测点在某一幅图上的位置拟合成一条平滑曲线。

        以弦长为参数，2 个点用直线，3 个点及以上用自然三次样条。
        返回 dict(x, y, t, knots)，点不够返回 None。
        """
        if len(self.points) < 2:
            return None
        p = np.array([q.pb if which == "b" else q.pa for q in self.points], float)
        keep = [0]
        for i in range(1, len(p)):
            if np.hypot(*(p[i] - p[keep[-1]])) > 1e-6:
                keep.append(i)
        p = p[keep]
        if len(p) < 2:
            return None
        d = np.r_[0.0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]
        if d[-1] <= 0:
            return None
        t = d / d[-1]
        tt = np.linspace(0.0, 1.0, int(n))
        if len(p) == 2:
            x = np.interp(tt, t, p[:, 0])
            y = np.interp(tt, t, p[:, 1])
        else:
            from scipy.interpolate import CubicSpline
            x = CubicSpline(t, p[:, 0], bc_type="natural")(tt)
            y = CubicSpline(t, p[:, 1], bc_type="natural")(tt)
        return dict(x=x, y=y, t=tt, knots=p)

    def curve_B(self, n=400):
        return self.curve_of("b", n)

    def curve_A(self, n=400):
        return self.curve_of("a", n)

    def arc_length_a(self):
        """已测点按 AIA 像面 xy 的累计弧长，角秒。"""
        if len(self.points) < 1:
            return np.zeros(0)
        p = np.array([q.pa for q in self.points], float)
        return np.r_[0.0, np.cumsum(np.hypot(*np.diff(p, axis=0).T))]

    def auto_snap(self, tol=4.0, gap_max_km=1500.0):
        """用已测点的平滑曲线吸附新的对映点。

        起点可以在任一格。取另一格已有的点拟合成曲线，与当前极线求最近点，
        把那一格的对映点吸到曲线上。至少要有 2 个已测点。
        """
        if len(self.points) < 2:
            self.message = "至少要有 2 个已测点，才能拟合曲线。"
            return None
        if self.pending is None:
            self.message = "先在任一格点一个起点，再按 M 让曲线吸附。"
            return None
        src = self.pending.get("src", "a")
        dst = "b" if src == "a" else "a"
        c = self.curve_of(dst)
        if c is None:
            self.message = "点重合太多，拟合不出曲线。"
            return None
        epi = self.pending["epi"]
        if epi is None:
            return None
        dx = c["x"][:, None] - epi["x"][None, :]
        dy = c["y"][:, None] - epi["y"][None, :]
        d = np.hypot(dx, dy)
        j, i = np.unravel_index(int(np.argmin(d)), d.shape)
        dmin = float(d[j, i])
        bx, by = float(c["x"][j]), float(c["y"][j])
        if not self.win[dst].contains(bx, by):
            self.message = f"吸附点落在 ({dst}) 视场外，先缩小视场再试。"
            return None
        if dmin > tol:
            self.message = (f"拟合曲线离极线最近 {dmin:.1f}″，超过 {tol:.0f}″，"
                            f"没有吸附。先手点一个，或再多测几个点把曲线延长。")
            return None
        px, py = self.pending["x"], self.pending["y"]
        if src == "a":
            pa, pb = (px, py), (bx, by)
        else:
            pa, pb = (bx, by), (px, py)
        res = G.triangulate(self.a.m, self.b.m, pa, pb, self.frame)
        if res["gap_km"] > gap_max_km:
            self.message = (f"曲线与极线并不真正相交，吸附后失配 "
                            f"{res['gap_km']:.0f} km，已放弃。")
            return None
        p = self.complete_point(dst, bx, by)
        if p is not None:
            p.note = "curve snap (smooth curve near epipolar)"
            self.message = (f"曲线吸附，离极线 {dmin:.1f}″，"
                            f"r = {p.r:.4f} R_sun，高度 {p.res['height_mm']:.1f} Mm")
        return p

    # ------------------------------------------------------------------ #
    # 取点
    # ------------------------------------------------------------------ #
    def start_point(self, key, x, y, guide=None):
        """从 key 这一格起一个新点，并把极线画到另一格。"""
        src = key
        dst = "b" if key == "a" else "a"
        m_src, m_dst = ((self.a.m, self.b.m) if key == "a"
                        else (self.b.m, self.a.m))
        self.pending = dict(src=src, x=float(x), y=float(y), guide=guide,
                            epi=G.epipolar(m_src, m_dst, x, y, self.frame,
                                           r_lo=EPI_R[0], r_hi=EPI_R[1], n=161))
        self.message = f"已画极线。到 ({dst}) 上点同一条结构的对映点。"
        return self.pending

    def complete_point(self, key, x, y):
        """在 key 这一格点对映点，与待确认的另一点三角化。"""
        if self.pending is None:
            self.message = "先在某一格点一个起点。"
            return None
        src = self.pending.get("src", "a")
        if src == key:
            self.message = "同一格不能配对，请到另一格点对映点。"
            return None
        px, py = self.pending["x"], self.pending["y"]
        if src == "a":
            pa, pb = (px, py), (float(x), float(y))
        else:
            pa, pb = (float(x), float(y)), (px, py)
        res = G.triangulate(self.a.m, self.b.m, pa, pb, self.frame)
        gi = self.pending.get("guide")
        if gi is not None and 0 <= gi < len(self.guides):
            g = self.guides[gi]
            if g.point is not None and g.point in self.points:
                p = g.point
                p.pa, p.pb, p.res = pa, pb, res
            else:
                p = Point(i=len(self.points) + 1, pa=pa, pb=pb, res=res,
                          seq=gi)
                self.points.append(p)
                g.point = p
            g.dst = (float(x), float(y))
            self._sort_points_by_seq()
            self.message = (f"引导点 {gi + 1}/{len(self.guides)} 已记录，"
                            f"r = {res['r']:.4f} R_sun，"
                            f"高度 {res['height_mm']:.1f} Mm，"
                            f"失配 {res['gap_km']:.1f} km。"
                            f"再点一次可覆盖。")
            return p
        p = Point(i=len(self.points) + 1, pa=pa, pb=pb, res=res)
        self.points.append(p)
        self.pending = None
        self.message = (f"点 {p.i}: r = {res['r']:.4f} R_sun，"
                        f"高度 {res['height_mm']:.1f} Mm，"
                        f"失配 {res['gap_km']:.1f} km")
        return p

    def pick_a(self, x, y):
        if self.pending is not None and self.pending.get("src", "a") == "b":
            return self.complete_point("a", x, y)
        return self.start_point("a", x, y)

    def pick_b(self, x, y):
        if self.pending is not None and self.pending.get("src", "a") == "a":
            return self.complete_point("b", x, y)
        return self.start_point("b", x, y)

    def delete_at(self, key, x, y, tol=14.0):
        """删掉离 (x,y) 最近且在该图上的点。"""
        best, bd = None, tol
        for k, p in enumerate(self.points):
            q = p.pa if key == "a" else p.pb
            d = float(np.hypot(q[0] - x, q[1] - y))
            if d < bd:
                best, bd = k, d
        if best is not None:
            p = self.points.pop(best)
            if p.seq is not None and 0 <= p.seq < len(self.guides):
                self.guides[p.seq].point = None
                self.guides[p.seq].dst = None
            for j, q in enumerate(self.points, start=1):
                q.i = j
            self.message = f"已删除点 {p.i}（r = {p.r:.4f}）。"
            return True
        return False

    def clear_points(self):
        self.points = []
        self.pending = None
        for g in self.guides:
            g.point = None
            g.dst = None
        if self.guides:
            self.goto_guide(self.guide_idx or 0)
            self.message = "点已清空，预置点列表还在。"
        else:
            self.message = "点已清空。"

    # ------------------------------------------------------------------ #
    # 拟合与统计
    # ------------------------------------------------------------------ #
    def pca(self):
        if len(self.points) < 3:
            return None
        xyz = np.array([p.res["xyz"] for p in self.points])
        return G.pca_fit(xyz)

    def summary(self):
        if not self.points:
            return dict(n=0)
        rs = np.array([p.r for p in self.points])
        gaps = np.array([p.res["gap_km"] for p in self.points])
        arc = self.arc_length_a()
        return dict(n=len(self.points), r_median=float(np.median(rs)),
                    r_min=float(rs.min()), r_max=float(rs.max()),
                    height_median_mm=float(np.median(G.height_mm(rs))),
                    height_min_mm=float(np.min(G.height_mm(rs))),
                    height_max_mm=float(np.max(G.height_mm(rs))),
                    gap_median_km=float(np.median(gaps)),
                    gap_max_km=float(gaps.max()),
                    arc_length_arcsec=float(arc[-1]))

    # ------------------------------------------------------------------ #
    # 保存
    # ------------------------------------------------------------------ #
    def stem(self):
        t = self.a.obstime
        return t.strftime("%Y-%m-%d T %H:%M:%S") + " UT SCC"

    def save(self, make_figure=True, dpi=300, size=(12.4, 9.0)):
        from . import render

        os.makedirs(self.outdir, exist_ok=True)
        stem = self.stem()
        paths = []
        if make_figure:
            fig = render.make_figure(self, size=size, dpi=dpi)
            for ext in ("png", "pdf"):
                p = os.path.join(self.outdir, f"{stem}.{ext}")
                fig.savefig(p, dpi=dpi, facecolor=fig.get_facecolor())
                paths.append(p)
        js = self.to_json()
        p = os.path.join(self.outdir, f"{stem}.json")
        with open(p, "w") as fh:
            json.dump(js, fh, indent=1, ensure_ascii=False)
        paths.append(p)
        p = os.path.join(self.outdir, f"{stem}.csv")
        self.to_csv(p)
        paths.append(p)
        return paths

    def to_json(self):
        d = dict(
            program="SCC Measure", version=__version__,
            obs1=dict(path=self.a.path, instrument=self.a.label,
                      time=self.a.obstime.isot, scale_arcsec=self.a.scale,
                      observer_lon=self.a.observer.transform_to(self.frame).lon.deg,
                      observer_lat=self.a.observer.transform_to(self.frame).lat.deg),
            obs2=dict(path=self.b.path, instrument=self.b.label,
                      time=self.b.obstime.isot, scale_arcsec=self.b.scale,
                      observer_lon=self.b.observer.transform_to(self.frame).lon.deg,
                      observer_lat=self.b.observer.transform_to(self.frame).lat.deg),
            separation_deg=self.separation,
            epipolar_sensitivity_arcsec_per_0p1Rsun=self.sensitivity,
            window_a=self.win["a"].as_tuple(), window_b=self.win["b"].as_tuple(),
            points=[p.as_dict() for p in self.points],
            summary=self.summary(),
        )
        return d

    def to_csv(self, path):
        with open(path, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["i", "aia_x_arcsec", "aia_y_arcsec", "euvi_x_arcsec", "euvi_y_arcsec",
                        "heeq_x_rsun", "heeq_y_rsun", "heeq_z_rsun", "r_rsun", "lon_deg",
                        "lat_deg", "height_mm", "ray_gap_km", "guide"])
            for p in self.points:
                d = p.as_dict()
                w.writerow([d["i"], *d["pa"], *d["pb"], *d["heeq"], d["r"], d["lon"],
                            d["lat"], d["height_mm"], d["gap_km"],
                            "" if d["guide"] is None else d["guide"]])

    # ------------------------------------------------------------------ #
    def status_text(self):
        s = [f"{self.a.short} {self.a.obstime.strftime('%H:%M:%S')}",
             f"{self.b.short} {self.b.obstime.strftime('%H:%M:%S')}",
             f"视角差 {self.separation:.2f}°   Δt {abs((self.a.obstime-self.b.obstime).to_value('s')):.0f} s",
             f"极线 {self.sensitivity:.0f}\"/0.1 R☉"]
        s.append(self.message)
        if self.guides:
            done = sum(1 for g in self.guides if g.point is not None)
            s.append(f"引导点 {(self.guide_idx or 0) + 1}/{len(self.guides)}，已取 {done}")
        if self.points:
            sm = self.summary()
            s += [f"点数 {sm['n']}",
                  f"r 中位 {sm['r_median']:.4f}",
                  f"高度中位 {sm['height_median_mm']:.0f} Mm",
                  f"失配中位 {sm['gap_median_km']:.0f} km"]
            s += [f"像面弧长 {sm['arc_length_arcsec']:.0f}\""]
        return "\n".join(s)
