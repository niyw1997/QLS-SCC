"""SCC Measure 的 Tk 界面。

左边一列控件，右边四格：两幅像、HEEQ 三维、PCA 局部平面。
交互照 pcp3 的两步流程，但把"看哪一块"和"点哪个点"分开：
缩放联动（Z）先对好视场，识别（I）给出候选轮廓，然后左键按 A→B 取点。
"""

from __future__ import annotations

import math
import os
import tkinter as tk
from tkinter import filedialog, messagebox

import matplotlib
matplotlib.use("TkAgg")
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
from matplotlib.widgets import RectangleSelector

from . import render
from .app import Session

__all__ = ["SCCApp", "run"]

BG = "#161616"
FG = "#cfcfcf"

SEL_PROPS = dict(facecolor="#ffe14d", edgecolor="#ffffff", alpha=0.16, linewidth=1.0)


class _Button(tk.Label):
    """深色按钮（macOS 原生 tk.Button 不认 bg/fg）。"""

    def __init__(self, master, text, command, height=None, **kw):
        super().__init__(master, text=text, bg="#2d2d2d", fg="#f2f2f2",
                         font=("Helvetica", 11), padx=10,
                         pady=6 if height is None else height,
                         anchor="w", cursor="hand2", **kw)
        self._cmd = command
        self._on = False
        self._enabled = True
        self.bind("<Button-1>", lambda e: self._cmd() if self._enabled else None)
        self.bind("<Enter>", lambda e: self.configure(bg="#454545" if not self._on else "#8a6a00"))
        self.bind("<Leave>", lambda e: self.configure(bg=self._bg()))

    def _bg(self):
        return "#7a5c00" if self._on else "#2d2d2d"

    def set_on(self, on):
        self._on = bool(on)
        self.configure(bg=self._bg())

    def set_enabled(self, on):
        self._enabled = bool(on)
        self.configure(fg="#f2f2f2" if on else "#6a6a6a",
                       bg=self._bg() if on else "#222222",
                       cursor="hand2" if on else "arrow")


class _Range(tk.Canvas):
    """双柄对比度滑杆，数值直接写在两端。"""

    H = 46
    C_LO = "#3fd4ff"
    C_HI = "#ffe14d"

    def __init__(self, master, command, width=246, height=H):
        super().__init__(master, width=width, height=height, bg=BG,
                         highlightthickness=0, bd=0, cursor="hand2")
        self.command = command
        self.lo, self.hi, self.log = 1.0, 10.0, True
        self.t0, self.t1 = 0.0, 1.0
        self._drag = None
        self._pad = 12
        self._y0, self._y1 = 26, 38
        self.bind("<Button-1>", self._press)
        self.bind("<B1-Motion>", self._motion)
        self.bind("<ButtonRelease-1>", self._release)
        self.bind("<Configure>", lambda e: self._draw())
        self._draw()

    def v_of(self, t):
        t = max(0.0, min(1.0, t))
        if self.log:
            return math.exp(math.log(self.lo) + t * (math.log(self.hi) - math.log(self.lo)))
        return self.lo + t * (self.hi - self.lo)

    def t_of(self, v):
        v = max(self.lo, min(self.hi, float(v)))
        if self.log:
            return (math.log(v) - math.log(self.lo)) / (math.log(self.hi) - math.log(self.lo))
        return (v - self.lo) / (self.hi - self.lo)

    def set_span(self, lo, hi, log):
        self.lo, self.hi, self.log = float(lo), float(hi), bool(log)

    def set_values(self, vmin, vmax):
        self.t0, self.t1 = self.t_of(vmin), self.t_of(vmax)
        self._draw()

    def values(self):
        return self.v_of(self.t0), self.v_of(self.t1)

    def _x(self, t):
        w = int(self.winfo_width() or self["width"])
        return self._pad + t * (w - 2 * self._pad)

    def _draw(self):
        self.delete("all")
        w = int(self.winfo_width() or self["width"])
        x0, x1 = self._pad, w - self._pad
        y0, y1 = self._y0, self._y1
        xa, xb = self._x(self.t0), self._x(self.t1)
        self.create_rectangle(x0, y0, x1, y1, fill="#2f2f2f", outline="#5a5a5a")
        self.create_rectangle(xa, y0, xb, y1, fill="#6b6320", outline="")
        for f in (0.25, 0.5, 0.75):
            xx = x0 + f * (x1 - x0)
            self.create_line(xx, y0, xx, y0 + 3, fill="#555555")
        for x, col in ((xa, self.C_LO), (xb, self.C_HI)):
            self.create_rectangle(x - 5, y0 - 7, x + 5, y1 + 7,
                                  fill=col, outline="#ffffff", width=1)
        self.create_text(x0, 5, anchor="w", fill=self.C_LO,
                         text=f"{self.v_of(self.t0):.4g}", font=("Menlo", 8))
        self.create_text(x1, 5, anchor="e", fill=self.C_HI,
                         text=f"{self.v_of(self.t1):.4g}", font=("Menlo", 8))

    def _t_at(self, x):
        w = int(self.winfo_width() or self["width"])
        return max(0.0, min(1.0, (x - self._pad) / max(1, w - 2 * self._pad)))

    def _press(self, ev):
        t = self._t_at(ev.x)
        xa, xb = self._x(self.t0), self._x(self.t1)
        self._drag = "min" if abs(ev.x - xa) <= abs(ev.x - xb) else "max"
        self._move(t)

    def _motion(self, ev):
        if self._drag:
            self._move(self._t_at(ev.x))

    def _release(self, ev):
        self._drag = None

    def _move(self, t):
        gap = 0.004
        if self._drag == "min":
            self.t0 = min(t, self.t1 - gap)
        else:
            self.t1 = max(t, self.t0 + gap)
        self._draw()
        if self.command:
            self.command()


class SCCApp:
    def __init__(self, session: Session, dialogs=True):
        self.s = session
        self.dialogs = dialogs
        self.zoom_mode = False
        self.selectors = {}
        self.levels = {}
        self.manual = {}
        self._draw_pending = False
        self._pan = None
        self._slider = {}
        self._guard = False

        self.root = tk.Tk()
        self.root.title(f"SCC Measure · {self.s.a.obstime.strftime('%Y-%m-%d %H:%M:%S')} UT")
        self.root.configure(bg=BG)
        sw = self.root.winfo_screenwidth()
        sh = self.root.winfo_screenheight()
        w = min(1320, sw - 40)
        h = min(940, sh - 80)
        self.root.geometry(f"{w}x{h}+20+30")

        self._build_sidebar()
        self._build_canvas()
        self.root.protocol("WM_DELETE_WINDOW", self.close)
        self.root.lift()
        self.root.attributes("-topmost", True)
        self.root.after(400, lambda: self.root.attributes("-topmost", False))

    # ------------------------------------------------------------------ #
    def _build_sidebar(self):
        side = tk.Frame(self.root, bg=BG, width=292)
        side.pack(side="left", fill="y")
        side.pack_propagate(False)

        tk.Label(side, text="SCC Measure", bg=BG, fg="#ffe14d",
                 font=("Helvetica", 22, "bold")).pack(anchor="w", padx=14, pady=(6, 0))
        tk.Label(side, text="双视角几何三角化 · Nisticò 2023", bg=BG, fg="#aaaaaa",
                 font=("Helvetica", 9)).pack(anchor="w", padx=14, pady=(0, 4))

        tk.Label(side, text=self.s.a.obstime.strftime("%Y-%m-%d"), bg=BG, fg="white",
                 font=("Helvetica", 12, "bold")).pack(anchor="w", padx=14)
        tk.Label(side, text=f"{self.s.a.short} {self.s.a.obstime.strftime('%H:%M:%S')}",
                 bg=BG, fg="#dddddd", font=("Menlo", 9)).pack(anchor="w", padx=14)
        tk.Label(side, text=f"{self.s.b.short} {self.s.b.obstime.strftime('%H:%M:%S')}",
                 bg=BG, fg="#dddddd", font=("Menlo", 9)).pack(anchor="w", padx=14,
                                                              pady=(0, 3))

        self.btn_zoom = _Button(side, "联动缩放 Zoom / Link  (Z)", self.toggle_zoom)
        self.btn_zoom.pack(fill="x", padx=14, pady=2)
        self.btn_snap = _Button(side, "曲线吸附 Curve  (M)", self.auto_snap)
        self.btn_snap.pack(fill="x", padx=14, pady=2)
        _Button(side, "删最后 Undo  (U)", self.undo).pack(fill="x", padx=14, pady=2)
        _Button(side, "清空 Clear  (X)", self.clear).pack(fill="x", padx=14, pady=2)
        _Button(side, "保存 Save  (S)", self.save).pack(fill="x", padx=14, pady=2)
        _Button(side, "全盘 Reset  (R)", self.reset).pack(fill="x", padx=14, pady=2)

        rrow = tk.Frame(side, bg=BG)
        rrow.pack(fill="x", padx=14, pady=(4, 0))
        tk.Label(rrow, text="联动壳层 r_hi（1.0 到 r_hi）", bg=BG, fg=FG,
                 font=("Helvetica", 9)).pack(anchor="w")
        self.r_var = tk.DoubleVar(value=self.s.r_link_hi)
        self.r_scale = tk.Scale(rrow, from_=1.05, to=1.40, resolution=0.05,
                                orient="horizontal", variable=self.r_var,
                                bg=BG, fg="white", troughcolor="#3a3a3a",
                                highlightthickness=0, bd=0, length=230,
                                sliderlength=14, width=8, sliderrelief="flat",
                                font=("Helvetica", 7), command=self._on_r)
        self.r_scale.pack(fill="x")

        tk.Label(side, text="对比度 min / max", bg=BG, fg="#9a9a9a",
                 font=("Helvetica", 9)).pack(anchor="w", padx=14, pady=(2, 1))
        self.contrast = tk.Frame(side, bg=BG)
        self.contrast.pack(fill="x", padx=14)
        _Button(side, "对比度自动  (A)", self.auto_levels, height=3).pack(
            fill="x", padx=14, pady=(2, 2))

        self.info = tk.Text(side, height=9, width=32, bg="#0b0b0b", fg="#9fe8b0",
                            relief="flat", font=("Menlo", 8), wrap="word",
                            padx=6, pady=6, borderwidth=0)
        self.info.pack(fill="both", expand=True, padx=14, pady=(2, 4))
        self.info.configure(state="disabled")

        for k, fn in (("z", self.toggle_zoom), ("m", self.auto_snap),
                      ("u", self.undo),
                      ("x", self.clear), ("s", self.save), ("r", self.reset),
                      ("a", self.auto_levels), ("q", self.close)):
            self.root.bind(f"<KeyPress-{k}>", lambda e, f=fn: f())
        self.root.bind("<Escape>", lambda e: self.toggle_zoom(False))

    def _build_canvas(self):
        frame = tk.Frame(self.root, bg="black")
        frame.pack(side="left", fill="both", expand=True)
        self.fig = Figure(figsize=(10.0, 8.4), dpi=100, facecolor=render.BG)
        self.canvas = FigureCanvasTkAgg(self.fig, master=frame)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.canvas.mpl_connect("button_press_event", self.on_press)
        self.canvas.mpl_connect("button_release_event", self.on_release)
        self.canvas.mpl_connect("motion_notify_event", self.on_motion)
        self.canvas.mpl_connect("scroll_event", self.on_scroll)
        self._build_contrast()
        self.redraw(rebuild_sliders=False)

    # ------------------------------------------------------------------ #
    def _build_contrast(self):
        for w in self.contrast.winfo_children():
            w.destroy()
        self._slider = {}
        for key in ("a", "b"):
            v = self.s.views[key]
            tk.Label(self.contrast, text=f"({key}) {v.short}", bg=BG, fg="#e8e8e8",
                     font=("Helvetica", 9), anchor="w").pack(fill="x", pady=(2, 0))
            rg = _Range(self.contrast, command=lambda k=key: self._on_level(k))
            rg.pack(fill="x", pady=(0, 3))
            self._slider[key] = rg
        self._sync_contrast()

    def _span(self, key, lo=None, hi=None):
        """滑杆量程。有手动值时把它也包进来，免得柄被夹在两端。"""
        v = self.s.views[key]
        a_lo, a_hi = v.levels(self.s.win[key])
        lo = a_lo if lo is None else min(a_lo, lo)
        hi = a_hi if hi is None else max(a_hi, hi)
        if v.norm_kind == "log":
            return max(lo / 10.0, 1e-3), hi * 10.0, True
        return 0.0, hi * 1.4, False

    def _sync_contrast(self):
        """把手柄对齐到当前生效的 vmin/vmax。没手动过就跟自动值。"""
        if not self._slider:
            return
        self._guard = True
        try:
            for key, rg in self._slider.items():
                auto = self.s.views[key].levels(self.s.win[key])
                if self.manual.get(key) and key in self.levels:
                    lo, hi = self.levels[key]
                else:
                    lo, hi = auto
                rg.set_span(*self._span(key, lo, hi))
                rg.set_values(lo, hi)
        finally:
            self._guard = False

    def _manual_levels(self):
        return {k: v for k, v in self.levels.items() if self.manual.get(k)}

    def _on_level(self, key):
        if self._guard:
            return
        vmin, vmax = self._slider[key].values()
        self.levels[key] = (vmin, vmax)
        self.manual[key] = True
        self._request_draw()

    def auto_levels(self):
        self.levels = {}
        self.manual = {}
        self._sync_contrast()
        self._request_draw()

    # ------------------------------------------------------------------ #
    def redraw(self, rebuild_sliders=True):
        render.make_figure(self.s, fig=self.fig, levels=self._manual_levels())
        self._rebuild_selectors()
        self.canvas.draw()
        if rebuild_sliders:
            self._sync_contrast()
        self._set_info(self.s.status_text())
        self._update_buttons()

    def _update_buttons(self):
        if hasattr(self, "btn_snap"):
            self.btn_snap.set_enabled(len(self.s.points) >= 2)

    def _request_draw(self):
        if self._draw_pending:
            return
        self._draw_pending = True
        self.root.after(35, self._flush)

    def _flush(self):
        self._draw_pending = False
        self.redraw()

    def _rebuild_selectors(self):
        for sel in self.selectors.values():
            try:
                sel.disconnect_events()
            except Exception:
                pass
        self.selectors = {}
        axes = getattr(self.fig, "_sccm_axes", {})
        for ax in (axes.get("a"), axes.get("b")):
            if ax is None:
                continue
            sel = RectangleSelector(
                ax, lambda e0, e1, a=ax: self._on_zoom_box(a, e0, e1),
                useblit=False, button=[1], minspanx=4, minspany=4,
                spancoords="pixels", interactive=False, props=SEL_PROPS)
            sel.set_active(self.zoom_mode)
            self.selectors[ax] = sel

    # ------------------------------------------------------------------ #
    def _ax_key(self, ax):
        axes = getattr(self.fig, "_sccm_axes", {})
        for name in ("a", "b"):
            if ax is axes.get(name):
                return name
        return None

    def _on_zoom_box(self, ax, e0, e1):
        if not self.zoom_mode:
            return
        key = self._ax_key(ax)
        if key is None:
            return
        if None in (e0.xdata, e0.ydata, e1.xdata, e1.ydata):
            return
        x0, x1 = sorted([e0.xdata, e1.xdata])
        y0, y1 = sorted([e0.ydata, e1.ydata])
        if (x1 - x0) < 3 or (y1 - y0) < 3:
            return
        self.s.zoom(key, (x0, x1, y0, y1), link=True)
        self.zoom_mode = False
        self.btn_zoom.set_on(False)
        for sel in self.selectors.values():
            sel.set_active(False)
        self.redraw()

    def toggle_zoom(self, force=None):
        self.zoom_mode = (not self.zoom_mode) if force is None else bool(force)
        for sel in self.selectors.values():
            sel.set_active(self.zoom_mode)
        self.btn_zoom.set_on(self.zoom_mode)
        self.canvas.get_tk_widget().configure(
            cursor="crosshair" if self.zoom_mode else "")
        self.s.message = ("在任一格里拖框，另一格按三维体联动。"
                          if self.zoom_mode else "缩放关。")
        self._set_info(self.s.status_text())

    # ------------------------------------------------------------------ #
    def on_press(self, ev):
        if ev.inaxes is None:
            return
        key = self._ax_key(ev.inaxes)
        if key is None:
            return
        if ev.button == 3 or (ev.button == 1 and ev.key == "shift"):
            if self.s.delete_at(key, ev.xdata, ev.ydata):
                self.redraw()
            return
        if ev.button == 2:
            self._pan = dict(key=key, x=ev.xdata, y=ev.ydata)
            return
        if ev.button != 1:
            return
        if self.zoom_mode:
            return                    # 交给 RectangleSelector
        x, y = ev.xdata, ev.ydata
        if x is None or y is None:
            return
        if self.s.pending is not None:
            if key == self.s.pending.get("src", "a"):
                self.s.start_point(key, x, y)
                self.s.message = "已移动起点，极线已更新。"
            else:
                self.s.complete_point(key, x, y)
            self.redraw()
            return
        # 没有待确认点时，从当前这一格起新点，另一格自动画极线
        self.s.start_point(key, x, y)
        self.redraw()

    def on_release(self, ev):
        if self._pan is None:
            return
        key = self._pan["key"]
        if ev.xdata is not None:
            dx = self._pan["x"] - ev.xdata
            dy = self._pan["y"] - ev.ydata
            if abs(dx) + abs(dy) > 2:
                self.s.pan(key, dx, dy, link=True)
                self.redraw()
        self._pan = None

    def on_motion(self, ev):
        if self._pan is not None and ev.inaxes is not None:
            ax = ev.inaxes
            if None not in (ev.xdata, ev.ydata):
                dx = self._pan["x"] - ev.xdata
                dy = self._pan["y"] - ev.ydata
                ax.set_xlim(ax.get_xlim()[0] - dx, ax.get_xlim()[1] - dx)
                ax.set_ylim(ax.get_ylim()[0] - dy, ax.get_ylim()[1] - dy)
                self.canvas.draw_idle()

    def on_scroll(self, ev):
        if ev.inaxes is None:
            return
        key = self._ax_key(ev.inaxes)
        if key is None:
            return
        if ev.button == "up":
            f = 1 / 1.3
        else:
            f = 1.3
        w = self.s.win[key]
        cx = ev.xdata if ev.xdata is not None else w.cx
        cy = ev.ydata if ev.ydata is not None else w.cy
        from .views import Window
        nw = Window(cx + (w.cx - cx) * f, cy + (w.cy - cy) * f, w.half * f)
        nw = nw.clipped(self.s.views[key].extent())
        self.s.set_window(key, nw, link=True)
        self.redraw()

    # ------------------------------------------------------------------ #
    def auto_snap(self):
        self.s.auto_snap()
        self.redraw()

    def undo(self):
        if self.s.points:
            self.s.points.pop()
            for j, q in enumerate(self.s.points, start=1):
                q.i = j
            self.s.message = "已撤销最后一个点。"
            self.redraw()

    def clear(self):
        self.s.clear_points()
        self.redraw()

    def reset(self):
        self.s.reset()
        self.redraw()

    def _on_r(self, _=None):
        self.s.r_link_hi = float(self.r_var.get())
        self.redraw()

    # ------------------------------------------------------------------ #
    def save(self):
        outdir = self.s.outdir
        if self.dialogs:
            chosen = filedialog.askdirectory(title="选择保存目录", initialdir=outdir,
                                             mustexist=True)
            if not chosen:
                return
            outdir = chosen
            self.s.outdir = chosen
        try:
            paths = self.s.save(make_figure=True, dpi=300)
        except Exception as exc:
            if self.dialogs:
                messagebox.showerror("SCC Measure", f"保存失败：{exc}")
            print("[sccm] 保存失败", repr(exc))
            return
        self.s.message = "已保存\n" + "\n".join(os.path.basename(p) for p in paths)
        self._set_info(self.s.status_text())
        if self.dialogs:
            messagebox.showinfo("SCC Measure", "已保存\n" + "\n".join(paths))
        print("[sccm] 已保存:\n  " + "\n  ".join(paths))
        self.last_saved = paths

    def _set_info(self, text):
        self.info.configure(state="normal")
        self.info.delete("1.0", "end")
        self.info.insert("1.0", text)
        self.info.configure(state="disabled")

    def close(self):
        self.root.destroy()

    def run(self):
        self.root.mainloop()


def run(session, dialogs=True):
    app = SCCApp(session, dialogs=dialogs)
    app.run()
    return app
