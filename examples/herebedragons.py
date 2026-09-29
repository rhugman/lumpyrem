"""Helpers for the lumpyrem example notebooks.

Nothing here is part of lumpyrem.  It builds the synthetic landscape the
spatial notebooks run on, reads the example climate, and draws maps, so the
notebooks themselves can stay about the model.

The landscape is a 10 km x 15 km valley on a 250 m grid, 40 rows by 60
columns.  A river runs down the middle; the valley floor is flat and clayey,
the slopes loamy, the ridges sandy.  Land use follows the terrain: irrigated
crops on the flat floor near the river, pasture on the slopes, native
vegetation on the ridges.  A buried palaeochannel of deep regolith meanders
along the valley.  All of it is made from smooth random fields with a fixed
seed, so every run of every notebook sees the same landscape.
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xarray as xr
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

DATA_D = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
CLIMATE_CSV = os.path.join(DATA_D, "climate.csv")

NROW, NCOL, DX = 40, 60, 250.0

SOIL_CLASSES = ("sand", "loam", "clay")
SOIL_COLORS = ("#e8c872", "#b07d48", "#6b4a3a")
LANDUSE_CLASSES = ("native", "pasture", "crop")
LANDUSE_COLORS = ("#2f6b3a", "#a8c66c", "#e0a030")


# ---------------------------------------------------------------------------
# Climate
# ---------------------------------------------------------------------------

def load_climate(path: str = CLIMATE_CSV) -> pd.DataFrame:
    """The example daily climate as a DataFrame indexed by ``time``."""
    df = pd.read_csv(path, index_col=0, parse_dates=True)
    df.index.name = "time"
    return df


def season_flag(index: pd.DatetimeIndex, months: tuple[int, ...]) -> np.ndarray:
    """1.0 on days whose month is in ``months``, else 0.0."""
    return np.isin(index.month, months).astype(float)


# ---------------------------------------------------------------------------
# Landscape
# ---------------------------------------------------------------------------

def smooth_field(shape: tuple[int, int], length: float, rng: np.random.Generator) -> np.ndarray:
    """Gaussian-smoothed white noise scaled to [0, 1]; ``length`` in cells.

    Filtered on a grid twice the size and cropped, so the field does not wrap
    around at the edges the way a plain FFT filter would.
    """
    ny, nx = 2 * shape[0], 2 * shape[1]
    noise = rng.standard_normal((ny, nx))
    ky = np.fft.fftfreq(ny)[:, None]
    kx = np.fft.fftfreq(nx)[None, :]
    kernel = np.exp(-2.0 * (np.pi * length) ** 2 * (kx**2 + ky**2))
    field = np.real(np.fft.ifft2(np.fft.fft2(noise) * kernel))[: shape[0], : shape[1]]
    return (field - field.min()) / (field.max() - field.min())


def synthetic_landscape(nrow: int = NROW, ncol: int = NCOL, dx: float = DX,
                        seed: int = 42) -> xr.Dataset:
    """The example valley as a Dataset on ``(y, x)``.

    Variables: ``dem`` (m), ``soil`` and ``landuse`` (integer codes into
    :data:`SOIL_CLASSES` and :data:`LANDUSE_CLASSES`), ``deep_regolith``
    (bool) and ``unsat_thickness`` (m, ground to water table).
    """
    rng = np.random.default_rng(seed)
    y = (np.arange(nrow) + 0.5) * dx
    x = (np.arange(ncol) + 0.5) * dx
    yy, xx = np.meshgrid(y, x, indexing="ij")

    # A valley whose axis wanders a little, falling gently to the north.
    axis = x.mean() + 0.08 * np.ptp(x) * np.sin(2 * np.pi * yy / np.ptp(y))
    across = np.abs(xx - axis) / (0.5 * np.ptp(x))
    dem = 150.0 + 90.0 * across**1.5 + 20.0 * yy / np.ptp(y) \
        + 25.0 * smooth_field((nrow, ncol), 4.0, rng)
    height = (dem - dem.min()) / (dem.max() - dem.min())

    texture = height + 0.35 * (smooth_field((nrow, ncol), 3.0, rng) - 0.5)
    soil = (2 - np.digitize(texture, [0.28, 0.55])).astype(np.int64)   # 0 sand .. 2 clay

    use = height + 0.3 * (smooth_field((nrow, ncol), 2.5, rng) - 0.5)
    landuse = (2 - np.digitize(use, [0.25, 0.6])).astype(np.int64)  # 0 native .. 2 crop

    channel = np.abs(xx - axis - 0.12 * np.ptp(x) * np.sin(5 * np.pi * yy / np.ptp(y))) \
        < 0.09 * np.ptp(x)
    unsat = 4.0 + 0.45 * (dem - dem.min())

    return xr.Dataset(
        {
            "dem": (("y", "x"), dem, {"units": "m"}),
            "soil": (("y", "x"), soil, {"classes": ", ".join(SOIL_CLASSES)}),
            "landuse": (("y", "x"), landuse, {"classes": ", ".join(LANDUSE_CLASSES)}),
            "deep_regolith": (("y", "x"), channel),
            "unsat_thickness": (("y", "x"), unsat, {"units": "m"}),
        },
        coords={"y": ("y", y, {"units": "m"}), "x": ("x", x, {"units": "m"})},
    )


# The parameter tables notebook 02 builds step by step, for notebook 03 to reuse.
SOIL_PARAMS = pd.DataFrame(
    {"awc": [0.10, 0.15, 0.18],          # plant-available water, m per m of root depth
     "ks": [0.8, 0.2, 0.03],             # m/day
     "m": [0.6, 0.5, 0.4],
     "mflowmax": [0.0, 0.0, 0.01]},      # macropore flow cap, m/day: cracking clay only
    index=pd.Index(SOIL_CLASSES, name="soil"),
)
LANDUSE_PARAMS = pd.DataFrame(
    {"root_depth": [2.0, 1.0, 1.2],      # m
     "irrigvolfrac": [0.0, 0.0, 0.4]},   # irrigation tops the store up to this fraction
    index=pd.Index(LANDUSE_CLASSES, name="landuse"),
)
RDELAY_PER_M = 3.0                       # days of drainage delay per metre of unsaturated zone


def cell_parameters(land: xr.Dataset, soils: pd.DataFrame = SOIL_PARAMS,
                    landuse: pd.DataFrame = LANDUSE_PARAMS,
                    rdelay_per_m: float = RDELAY_PER_M) -> pd.DataFrame:
    """One row of upper-store parameters per cell, numbered row-major."""
    soil = soils.iloc[land.soil.values.ravel()].reset_index()
    use = landuse.iloc[land.landuse.values.ravel()].reset_index()
    cells = pd.concat([soil, use], axis=1)
    cells["maxvol"] = cells["awc"] * cells["root_depth"]
    cells["rdelay"] = rdelay_per_m * land.unsat_thickness.values.ravel()
    return cells


# ---------------------------------------------------------------------------
# Maps
# ---------------------------------------------------------------------------

def _extent(land: xr.Dataset) -> tuple[float, float, float, float]:
    """imshow extent in km, north up."""
    half = 0.5 * float(land.x[1] - land.x[0])
    return (float(land.x[0] - half) / 1e3, float(land.x[-1] + half) / 1e3,
            float(land.y[0] - half) / 1e3, float(land.y[-1] + half) / 1e3)


def to_map(values, land: xr.Dataset) -> np.ndarray:
    """A flat per-cell array (row-major, as the notebooks number cells) as a map."""
    return np.asarray(values).reshape(land.sizes["y"], land.sizes["x"])


def plot_map(ax, values, land: xr.Dataset, *, title: str = "", label: str = "",
             cmap: str = "viridis", **kwargs):
    """Draw one continuous field, flat or 2-D, with a colour bar."""
    im = ax.imshow(to_map(values, land), origin="lower", extent=_extent(land),
                   cmap=cmap, **kwargs)
    ax.set_title(title)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")
    plt.colorbar(im, ax=ax, label=label, shrink=0.85)
    return im


def plot_classes(ax, codes, land: xr.Dataset, names, colors, *, title: str = ""):
    """Draw a class map with a legend."""
    cmap = ListedColormap(colors)
    norm = BoundaryNorm(np.arange(len(names) + 1) - 0.5, cmap.N)
    ax.imshow(to_map(codes, land), origin="lower", extent=_extent(land),
              cmap=cmap, norm=norm)
    ax.legend(handles=[Patch(color=c, label=n) for n, c in zip(names, colors, strict=True)],
              loc="upper right", fontsize="small", framealpha=0.9)
    ax.set_title(title)
    ax.set_xlabel("x (km)")
    ax.set_ylabel("y (km)")


def mark_cells(ax, cells, land: xr.Dataset, labels=None):
    """Mark flat cell numbers on a map, optionally labelled."""
    ncol = land.sizes["x"]
    dx = float(land.x[1] - land.x[0]) / 1e3
    for k, cell in enumerate(cells):
        row, col = divmod(int(cell), ncol)
        xk, yk = (col + 0.5) * dx, (row + 0.5) * dx
        ax.plot(xk, yk, "o", mfc="white", mec="black", ms=7)
        if labels is not None:
            ax.annotate(labels[k], (xk, yk), xytext=(5, 5), textcoords="offset points",
                        fontsize="small", fontweight="bold",
                        bbox=dict(boxstyle="round,pad=0.15", fc="white", alpha=0.8))
