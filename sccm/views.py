"""读 FITS、按 HPC 网格重采样、拉伸与显示窗口。

与 QLS 的区别只有一处：这里的像可能带 PC 旋转（EUVI 304 的 CROTA 约 3°），
所以不做"裁像素再 imshow"，而是每次按 HPC 角秒网格重采样，
这样轴刻度、叠加的极线与取点都在同一个坐标里，旋转不会偷偷跑掉。
"""

from __future__ import annotations

import os
from dataclasses import dataclass, field

import numpy as np
import astropy.units as u
from astropy.coordinates import SkyCoord
from scipy import ndimage as ndi
from sunpy.map import Map

try:                                     # 注册 sdoaia304 / euvi304 色表
    import sunpy.visualization.colormaps  # noqa: F401
except Exception:                        # pragma: no cover
    pass

__all__ = ["Window", "View", "load_view", "guess_key"]


# --------------------------------------------------------------------------- #
# 视场
# --------------------------------------------------------------------------- #
@dataclass
class Window:
    """方形视场，单位角秒（HPC）。"""

    cx: float
    cy: float
    half: float

    @property
    def x0(self):
        return self.cx - self.half

    @property
    def x1(self):
        return self.cx + self.half

    @property
    def y0(self):
        return self.cy - self.half

    @property
    def y1(self):
        return self.cy + self.half

    @property
    def side(self):
        return 2.0 * self.half

    def as_tuple(self):
        return (self.x0, self.x1, self.y0, self.y1)

    @classmethod
    def from_box(cls, box, square=True):
        x0, x1, y0, y1 = [float(v) for v in box]
        cx, cy = 0.5 * (x0 + x1), 0.5 * (y0 + y1)
        if square:
            half = 0.5 * max(abs(x1 - x0), abs(y1 - y0))
        else:
            half = 0.5 * max(abs(x1 - x0), abs(y1 - y0))
        return cls(cx, cy, max(half, 1.0))

    def contains(self, x, y):
        return self.x0 <= x <= self.x1 and self.y0 <= y <= self.y1

    def clipped(self, lim):
        """把中心夹在 lim=(xmin,xmax,ymin,ymax) 里，保证窗口不跑出像。"""
        if lim[1] - lim[0] <= 2 * self.half or lim[3] - lim[2] <= 2 * self.half:
            return self
        return Window(min(max(self.cx, lim[0] + self.half), lim[1] - self.half),
                      min(max(self.cy, lim[2] + self.half), lim[3] - self.half),
                      self.half)

    def __str__(self):
        return (f"({self.cx:.0f},{self.cy:.0f}) ± {self.half:.0f}\"")


# --------------------------------------------------------------------------- #
# 单幅像
# --------------------------------------------------------------------------- #
_DEFAULT = {
    # key: (label, short, cmap, p_lo, p_hi, norm)
    "aia": ("SDO/AIA 304 Å", "AIA 304", "sdoaia304", 0.5, 99.6, "log"),
    "euvi": ("STEREO-A/EUVI 304 Å", "EUVI 304", "euvi304", 2.0, 99.2, "log"),
    "aia193": ("SDO/AIA 193 Å", "AIA 193", "sdoaia193", 1.0, 99.7, "log"),
    "euvi195": ("STEREO-A/EUVI 195 Å", "EUVI 195", "euvi195", 2.0, 99.2, "log"),
    "other": ("FITS", "FITS", "sdoaia304", 1.0, 99.5, "log"),
}


def guess_key(path):
    """按文件名与头猜仪器波段，只用于配色和标签。"""
    name = os.path.basename(str(path)).lower()
    if "aia" in name:
        if "193" in name:
            return "aia193"
        if "304" in name:
            return "aia"
    elif "euvi" in name or "n5eu" in name or "n4eu" in name:
        if "195" in name:
            return "euvi195"
        if "304" in name:
            return "euvi"
    try:
        h = Map(path).fits_header
        tel = str(h.get("TELESCOP", "")).lower()
        wl = int(float(h.get("WAVELNTH", 0) or 0))
        if "aia" in tel:
            return "aia193" if wl == 193 else "aia"
        if "stereo" in tel or "euvi" in str(h.get("INSTRUME", "")).lower():
            return "euvi195" if wl == 195 else "euvi"
    except Exception:
        pass
    return "other"


@dataclass
class View:
    key: str
    path: str
    m: Map
    label: str
    short: str
    cmap: str = "sdoaia304"
    norm_kind: str = "log"
    gamma: float = 0.4
    p_lo: float = 0.5
    p_hi: float = 99.6
    obstacle: str = ""
    _dark: np.ndarray | None = field(default=None, repr=False)
    _dark_sigma: float | None = field(default=None, repr=False)
    _cache: dict = field(default_factory=dict, repr=False)

    # ---- 基本量 -------------------------------------------------------- #
    @property
    def data(self):
        return self.m.data

    @property
    def obstime(self):
        return self.m.date

    @property
    def shape(self):
        return self.m.data.shape

    @property
    def scale(self):
        """原始像素尺度，角秒/像素。"""
        return float(self.m.scale.axis1.to_value(u.arcsec / u.pix))

    @property
    def rsun_arcsec(self):
        # EUVI 的头用 RSUN（角秒），sunpy 有时读不到，这里显式补一下
        h = self.m.fits_header
        for key in ("RSUN", "RSUN_OBS"):
            if key in h:
                try:
                    v = float(h[key])
                    if v > 0:
                        return v
                except Exception:
                    pass
        return float(getattr(self.m, "rsun_arcsec", 960.0))

    @property
    def observer(self):
        return self.m.observer_coordinate

    # ---- 坐标 ---------------------------------------------------------- #
    def to_world(self, px, py):
        sc = self.m.pixel_to_world(np.atleast_1d(px) * u.pix,
                                   np.atleast_1d(py) * u.pix)
        return (np.asarray(sc.Tx.to_value(u.arcsec), float),
                np.asarray(sc.Ty.to_value(u.arcsec), float))

    def to_pixel(self, x, y):
        sc = SkyCoord(np.atleast_1d(x) * u.arcsec, np.atleast_1d(y) * u.arcsec,
                      frame=self.m.coordinate_frame)
        px, py = self.m.world_to_pixel(sc)
        px = u.Quantity(px).to_value(u.pix) if hasattr(px, "unit") else px
        py = u.Quantity(py).to_value(u.pix) if hasattr(py, "unit") else py
        return np.asarray(px, float), np.asarray(py, float)

    def extent(self):
        """整幅像在 HPC 上的外接矩形，(x0, x1, y0, y1)。"""
        nx, ny = self.m.data.shape[1], self.m.data.shape[0]
        px = np.array([-0.5, nx - 0.5, nx - 0.5, -0.5])
        py = np.array([-0.5, -0.5, ny - 0.5, ny - 0.5])
        tx, ty = self.to_world(px, py)
        return (float(np.nanmin(tx)), float(np.nanmax(tx)),
                float(np.nanmin(ty)), float(np.nanmax(ty)))

    def full_window(self):
        # 两格都按 1.06 倍视半径开窗，太阳看起来一样大，也不把整片天空塞进来
        return Window(0.0, 0.0, 1.06 * self.rsun_arcsec)

    def contains(self, x, y):
        x0, x1, y0, y1 = self.extent()
        return x0 <= x <= x1 and y0 <= y <= y1

    # ---- 显示 ---------------------------------------------------------- #
    def subimage(self, win, max_px=760):
        """按窗口重采样成北朝上的 HPC 网格。返回 (img, extent)。"""
        key = (round(win.cx, 2), round(win.cy, 2), round(win.half, 2), int(max_px))
        if key in self._cache:
            return self._cache[key]
        n = int(np.clip(win.side / max(0.6, 0.8 * self.scale), 160, max_px))
        gx = np.linspace(win.x0, win.x1, n)
        gy = np.linspace(win.y0, win.y1, n)
        GX, GY = np.meshgrid(gx, gy)
        px, py = self.to_pixel(GX.ravel(), GY.ravel())
        img = ndi.map_coordinates(np.asarray(self.m.data, float),
                                  [py, px], order=1, mode="constant", cval=np.nan)
        img = img.reshape(n, n)
        out = (img, (win.x0, win.x1, win.y0, win.y1))
        if len(self._cache) > 12:
            self._cache.clear()
        self._cache[key] = out
        return out

    def levels(self, win=None, n=220):
        """窗口内的 log 拉伸上下限，只在日面内取分位数。"""
        win = win or self.full_window()
        sub, ext = self.subimage(win, max_px=n)
        gx = np.linspace(ext[0], ext[1], sub.shape[1])
        gy = np.linspace(ext[2], ext[3], sub.shape[0])
        GX, GY = np.meshgrid(gx, gy)
        on_disk = np.hypot(GX, GY) < 0.97 * self.rsun_arcsec
        good = np.isfinite(sub) & (sub > 0) & on_disk
        vals = sub[good]
        if vals.size < 32:
            return 1.0, 10.0
        lo, hi = np.nanpercentile(vals, [self.p_lo, self.p_hi])
        if not np.isfinite(lo) or lo <= 0:
            lo = float(np.nanmin(vals))
        if not np.isfinite(hi) or hi <= lo:
            hi = lo * 10.0
        return float(lo), float(hi)

    # ---- 暗化图（识别用） ---------------------------------------------- #
    def darkening(self, sigma_arcsec=60.0, bg_median=0.0):
        """原图减去大尺度高斯背景后的相对暗化，(bg - I)/bg，正值是暗结构。"""
        if self._dark is not None and self._dark_sigma == sigma_arcsec:
            return self._dark
        img = np.asarray(self.m.data, float)
        sig = max(sigma_arcsec / max(self.scale, 1e-6), 3.0)
        bg = ndi.gaussian_filter(img, sig, mode="nearest")
        with np.errstate(invalid="ignore", divide="ignore"):
            d = (bg - img) / np.where(bg > 0, bg, np.nan)
        self._dark = d
        self._dark_sigma = sigma_arcsec
        return d

    def region_mask(self, win, shrink=0.98):
        """窗口对应的原生像素布尔掩膜（含日面限制）。"""
        ny, nx = self.m.data.shape
        gx = np.linspace(win.x0, win.x1, 5)
        gy = np.linspace(win.y0, win.y1, 5)
        GX, GY = np.meshgrid(gx, gy)
        px, py = self.to_pixel(GX.ravel(), GY.ravel())
        i0, i1 = np.clip([px.min(), px.max()], 0, nx - 1).astype(int)
        j0, j1 = np.clip([py.min(), py.max()], 0, ny - 1).astype(int)
        m = np.zeros((ny, nx), bool)
        m[max(j0, 0):min(j1 + 1, ny), max(i0, 0):min(i1 + 1, nx)] = True
        # 日面内：用原生像素坐标算到日心的距离
        xx, yy = np.meshgrid(np.arange(nx), np.arange(ny))
        cx, cy = self.to_pixel(0.0, 0.0)
        r = np.hypot(xx - float(np.ravel(cx)[0]), yy - float(np.ravel(cy)[0]))
        rsun_px = self.rsun_arcsec / self.scale
        m &= r < shrink * rsun_px
        m &= np.isfinite(self.m.data) & (self.m.data > 0)
        return m


# --------------------------------------------------------------------------- #
# 读入
# --------------------------------------------------------------------------- #
def load_view(path, key=None, label=None, cmap=None):
    key = key or guess_key(path)
    lab, short, cmap0, p_lo, p_hi, norm = _DEFAULT.get(key, _DEFAULT["other"])
    m = Map(path)
    v = View(key=key, path=str(path), m=m, label=label or lab, short=short,
             cmap=cmap or cmap0, norm_kind=norm, p_lo=p_lo, p_hi=p_hi)
    if key in ("euvi", "euvi195"):
        v.norm_kind = "power"
        v.gamma = 0.35
        v.p_lo, v.p_hi = 1.0, 99.7
    return v
