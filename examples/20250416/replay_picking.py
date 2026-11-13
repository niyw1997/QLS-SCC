#!/usr/bin/env python3
"""Replay a recorded measurement and stop in the middle of point picking.

The CSV is the one exported by SCC Measure. Each row has the AIA and EUVI
coordinates of one measured point. This script adds the first N points, starts
point N+1 on AIA, and draws its epipolar line on EUVI, so the state is exactly
a frame from the picking process.
"""
import argparse
import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from sccm.app import Session          # noqa: E402
from sccm.views import Window        # noqa: E402


def read_points(path):
    pts = []
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            pa = (float(row["aia_x_arcsec"]), float(row["aia_y_arcsec"]))
            pb = (float(row["euvi_x_arcsec"]), float(row["euvi_y_arcsec"]))
            pts.append((pa, pb))
    return pts


def build_session(points, csv_path, fov_path, n):
    import json
    data = os.path.join(HERE, "data")
    with open(fov_path, encoding="utf-8") as fh:
        fov = json.load(fh)["mapping"]["windows"]
    s = Session(os.path.join(data, "aia193_20250416_cutout.fits"),
                os.path.join(data, "euvi195_20250416_cutout.fts"),
                outdir=os.path.join(HERE, "output"), verbose=False)
    s.win["a"] = Window.from_box(fov["aia"])
    s.win["b"] = Window.from_box(fov["euvi"])
    for pa, pb in points[:n]:
        s.pick_a(*pa)
        s.pick_b(*pb)
    if n < len(points):
        s.start_point("a", *points[n][0])
        s.message = (f"回放 {os.path.basename(csv_path)} 前 {n} 个点，"
                     f"第 {n + 1} 个点已起，极线已画到另一格。")
    else:
        s.message = f"回放 {os.path.basename(csv_path)} 全部 {n} 个点。"
    return s


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--points", type=int, default=4,
                    help="先放回几个已测点（默认 4），再起下一个点")
    ap.add_argument("--csv", default=os.path.join(HERE, "measurement.csv"))
    ap.add_argument("--fov", default=os.path.join(HERE, "fov.json"))
    ap.add_argument("--save", default=None, help="保存一张四格图，不打开界面")
    args = ap.parse_args()
    points = read_points(args.csv)
    n = max(0, min(args.points, len(points)))
    s = build_session(points, args.csv, args.fov, n)
    if args.save:
        from sccm import render
        fig = render.make_figure(s, size=(10.0, 8.4), dpi=150)
        fig.savefig(args.save, dpi=150)
        print(f"wrote {args.save}")
        return 0
    from sccm import gui
    gui.run(s, dialogs=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
