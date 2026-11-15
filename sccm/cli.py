"""命令行入口。"""

from __future__ import annotations

import argparse
import csv
import json
import os
import sys

import numpy as np

from . import __version__, geometry as G
from .app import Point, Session

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DEF_A = os.path.join(HERE, "examples/20250416/data/aia193_20250416_cutout.fits")
DEF_B = os.path.join(HERE, "examples/20250416/data/euvi195_20250416_cutout.fts")
DEF_FOV = os.path.join(HERE, "examples/20250416/fov.json")


def _box(s):
    if not s:
        return None
    v = [float(x) for x in s.replace(" ", "").split(",")]
    if len(v) != 4:
        raise argparse.ArgumentTypeError("视场要写成 x0,x1,y0,y1（角秒）")
    return tuple(v)


def _cli_info(args):
    a = Session(args.obs1 or DEF_A, args.obs2 or DEF_B, outdir=args.out,
                verbose=False)
    print(f"SCC Measure {__version__}")
    print(f"  (a) {a.a.label}  {a.a.obstime}  尺度 {a.a.scale:.3f}\"/px  "
          f"DSUN {a.a.m.fits_header.get('DSUN_OBS', float('nan'))/1e9:.3f} Gm")
    print(f"  (b) {a.b.label}  {a.b.obstime}  尺度 {a.b.scale:.3f}\"/px  "
          f"DSUN {a.b.m.fits_header.get('DSUN_OBS', float('nan'))/1e9:.3f} Gm")
    dt = abs((a.a.obstime - a.b.obstime).to_value("s"))
    print(f"  视角差 {a.separation:.3f}°   Δt {dt:.0f} s   "
          f"极线 {a.sensitivity:.1f}\"/0.1 R☉")
    print(f"  联动窗口 a {a.win['a']}  b {a.win['b']}")


def _cli_demo(args):
    from .demo import run_demo
    center = args.center
    run_demo(args.obs1 or DEF_A, args.obs2 or DEF_B, outdir=args.out,
             center=center, half=args.half, n=args.n)


def _apply_fov(s, path):
    """读 QLS 保存的 JSON，把 aia 与 euvi 两个窗口原样套到 SCC 上。"""
    from .views import Window
    with open(path, encoding="utf-8") as fh:
        j = json.load(fh)
    w = j.get("mapping", {}).get("windows", {})
    if "aia" not in w or "euvi" not in w:
        raise SystemExit(f"{path} 里没有 aia/euvi 窗口")
    s.win["a"] = Window.from_box(w["aia"])
    s.win["b"] = Window.from_box(w["euvi"])
    s.message = (f"已载入 FoV {os.path.basename(path)}\n"
                 f"AIA {s.win['a']} / EUVI {s.win['b']}\n"
                 f"在任一格点第一个点，另一格画极线。")


def _cli_show(args):
    from . import gui
    s = Session(args.obs1 or DEF_A, args.obs2 or DEF_B, outdir=args.out,
                r_link_hi=args.r_hi)
    if not args.fov_json and not args.zoom and args.obs1 is None and args.obs2 is None:
        args.fov_json = DEF_FOV
    if args.fov_json:
        _apply_fov(s, args.fov_json)
    elif args.zoom:
        s.zoom("a", _box(args.zoom), link=True)
    if args.guide:
        from .guides import read_guides
        try:
            g = read_guides(args.guide, args.guide_panel)
            s.load_guides(g["panel"], g["points"], windows=g.get("windows"))
        except Exception as exc:
            raise SystemExit(f"载入预置点失败: {exc}")
    gui.run(s, dialogs=not args.no_dialogs)


def _cli_selftest(args):
    from .selftest import run_selftest
    return run_selftest(args.obs1, args.obs2)


def _cli_points(args):
    s = Session(args.obs1 or DEF_A, args.obs2 or DEF_B, outdir=args.out, verbose=True)
    rows = []
    with open(args.pairs) as fh:
        for row in csv.reader(fh):
            if not row or row[0].startswith("#"):
                continue
            v = [float(x) for x in row[:4]]
            rows.append(v)
    print(f"{'#':>3} {'AIA x':>8} {'AIA y':>8} {'EUVI x':>8} {'EUVI y':>8} "
          f"{'r':>7} {'height':>8} {'gap':>8}")
    for i, v in enumerate(rows, 1):
        t = G.triangulate(s.a.m, s.b.m, (v[0], v[1]), (v[2], v[3]), s.frame)
        print(f"{i:3d} {v[0]:8.1f} {v[1]:8.1f} {v[2]:8.1f} {v[3]:8.1f} "
              f"{t['r']:7.4f} {t['height_mm']:8.1f} {t['gap_km']:8.1f}")
        s.points.append(Point(i=i, pa=(v[0], v[1]), pb=(v[2], v[3]), res=t))
    paths = s.save(make_figure=True, dpi=300)
    print("已保存:\n  " + "\n  ".join(paths))


def build_parser():
    p = argparse.ArgumentParser(
        prog="sccm", description="SCC Measure：双视角几何三角化工作台")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp):
        sp.add_argument("--obs1", default=None, help=f"参考视角 FITS（默认 {DEF_A}）")
        sp.add_argument("--obs2", default=None, help=f"另一视角 FITS（默认 {DEF_B}）")
        sp.add_argument("--out", default=os.path.join(HERE, "out"), help="输出目录")

    sp = sub.add_parser("show", help="打开交互界面")
    common(sp)
    sp.add_argument("--zoom", default=None, help="启动视场 x0,x1,y0,y1（角秒）")
    sp.add_argument("--fov-json", default=None,
                    help="QLS 保存的 JSON，直接套用里面的 aia/euvi 窗口")
    sp.add_argument("--guide", default=None,
                    help="预置点列表 JSON/CSV，打开后可逐点取对映点")
    sp.add_argument("--guide-panel", default=None,
                    choices=("a", "b", "aia", "euvi"),
                    help="预置点在 a/aia 还是 b/euvi")
    sp.add_argument("--r-hi", type=float, default=1.10, help="联动壳层上限，默认 1.10")
    sp.add_argument("--no-dialogs", action="store_true", help="不弹保存/报错对话框")
    sp.set_defaults(func=_cli_show)

    sp = sub.add_parser("demo", help="无窗口跑一遍并出图")
    common(sp)
    sp.add_argument("--center", type=float, nargs=2, default=(-335.0, 71.0),
                    metavar=("X", "Y"), help="演示视场中心（角秒）")
    sp.add_argument("--half", type=float, default=180.0, help="演示视场半宽（角秒）")
    sp.add_argument("-n", type=int, default=11, help="合成环取几个点")
    sp.set_defaults(func=_cli_demo)

    sp = sub.add_parser("info", help="打印两幅像的几何信息")
    common(sp)
    sp.set_defaults(func=_cli_info)

    sp = sub.add_parser("selftest", help="数值自检")
    sp.add_argument("--obs1", default=None)
    sp.add_argument("--obs2", default=None)
    sp.set_defaults(func=_cli_selftest)

    sp = sub.add_parser("points", help="批量三角化 CSV 里已配好的点")
    common(sp)
    sp.add_argument("pairs", help="CSV，每行 aia_x,aia_y,euvi_x,euvi_y")
    sp.set_defaults(func=_cli_points)
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
