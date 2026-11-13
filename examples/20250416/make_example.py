#!/usr/bin/env python3
"""Headless example: two tie-points on the 2025-04-16 193/195 pair.

Run from the repository root, or set PYTHONPATH to it. The GUI version is run.sh.
"""
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, ROOT)

from sccm.app import Session          # noqa: E402
from sccm.views import Window        # noqa: E402


def main():
    data = os.path.join(HERE, "data")
    out = os.path.join(HERE, "output")
    with open(os.path.join(HERE, "fov.json"), encoding="utf-8") as fh:
        fov = json.load(fh)["mapping"]["windows"]

    s = Session(os.path.join(data, "aia193_20250416_cutout.fits"),
                os.path.join(data, "euvi195_20250416_cutout.fts"),
                outdir=out, verbose=True)
    s.win["a"] = Window.from_box(fov["aia"])
    s.win["b"] = Window.from_box(fov["euvi"])

    # Two points on the same structure. The second one was obtained by
    # back-projecting a point on the EUVI contour at the same height.
    for pa, pb in (((-400.0, 100.0), (-943.0, 24.3)),
                   ((-459.6677, 130.9985), (-980.7457, 57.4983))):
        s.pick_a(*pa)
        s.pick_b(*pb)

    print("[example] summary:", s.summary())
    paths = s.save(make_figure=True, dpi=300)
    print("[example] wrote:\n  " + "\n  ".join(paths))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
