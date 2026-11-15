"""读预置点列表。

支持 JSON 和 CSV。JSON 可以直接给出 FoV 窗口，也可以把 QLS 的
``mapping.windows`` 原样放进 ``windows`` 或 ``mapping``。

示例 JSON

    {
      "panel": "a",
      "points": [[-379.3, 117.6], [-371.7, 101.4]],
      "windows": {
        "aia": [-513.4, -155.8, -107.4, 250.2],
        "euvi": [-1065.7, -681.6, -184.5, 199.5]
      }
    }

``panel`` 是源图，``a``/``aia`` 或 ``b``/``euvi``。points 是源图上的
HPC 角秒坐标。CSV 的表头可以是 ``x,y``，也可以是
``aia_x_arcsec,aia_y_arcsec,euvi_x_arcsec,euvi_y_arcsec``。
"""
from __future__ import annotations

import csv
import json
import os

__all__ = ["read_guides", "normalize_windows"]


def _panel(value, default=None):
    if value is None:
        return default
    s = str(value).strip().lower()
    if s in ("a", "aia", "sdo", "earth"):
        return "a"
    if s in ("b", "euvi", "stereo", "stereo-a"):
        return "b"
    raise ValueError(f"不认识的源图 {value!r}，请用 a/aia 或 b/euvi")


def normalize_windows(w):
    """把 aia/euvi 或 a/b 的窗口统一成 {'a':..., 'b':...}。"""
    if not w:
        return {}
    if "mapping" in w and isinstance(w["mapping"], dict):
        w = w["mapping"].get("windows", w)
    out = {}
    for key in ("a", "aia"):
        if key in w and w[key] is not None:
            out["a"] = [float(v) for v in w[key]]
            break
    for key in ("b", "euvi"):
        if key in w and w[key] is not None:
            out["b"] = [float(v) for v in w[key]]
            break
    return out


def _points_from_json(raw):
    pts = []
    for q in raw or []:
        if isinstance(q, dict):
            x = q.get("x", q.get("aia_x_arcsec", q.get("euvi_x_arcsec")))
            y = q.get("y", q.get("aia_y_arcsec", q.get("euvi_y_arcsec")))
        else:
            x, y = q[0], q[1]
        if x is None or y is None:
            continue
        pts.append((float(x), float(y)))
    return pts


def _read_csv(path, panel):
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise ValueError("CSV 是空的")
    keys = {k.strip() for k in rows[0]}
    p = panel
    if p is None and "panel" in keys:
        p = _panel(rows[0].get("panel"))
    p = p or "a"
    xkey = "x" if "x" in keys else ("aia_x_arcsec" if p == "a" else "euvi_x_arcsec")
    ykey = "y" if "y" in keys else ("aia_y_arcsec" if p == "a" else "euvi_y_arcsec")
    if xkey not in keys or ykey not in keys:
        raise ValueError("CSV 需要 x,y 或 aia_x_arcsec/euvi_x_arcsec 坐标列")
    pts = [(float(r[xkey]), float(r[ykey])) for r in rows
           if r.get(xkey) not in (None, "") and r.get(ykey) not in (None, "")]
    if not pts:
        raise ValueError("CSV 里没有有效坐标")
    return dict(panel=p, points=pts, windows={}, source=path)


def read_guides(path, panel=None):
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    ext = os.path.splitext(path)[1].lower()
    if ext == ".csv":
        return _read_csv(path, _panel(panel) if panel else None)
    with open(path, encoding="utf-8") as fh:
        j = json.load(fh)
    p = _panel(j.get("panel") or j.get("source") or panel or "a")
    pts = _points_from_json(j.get("points") or j.get("coords"))
    if not pts:
        raise ValueError("JSON 里没有 points 或 coords")
    w = j.get("windows") or j.get("fov")
    if not w and "mapping" in j:
        w = j["mapping"].get("windows")
    return dict(panel=p, points=pts, windows=normalize_windows(w), source=path)
