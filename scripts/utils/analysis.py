"""
analysis.py — time-series extraction and statistical analysis for ATS simulations.

Typical usage
-------------
>>> from utils.analysis import (
...     theil_sen_fit,
...     determine_ice_table_depth,
...     build_ice_table_df,
...     compute_winter_temp,
...     compute_winter_snow,
...     compute_ponded_depth_series,
...     compute_unfrozen_area_series,
... )
"""

from __future__ import annotations

from datetime import datetime, timedelta

import numpy as np
import pandas as pd
from scipy.stats import kendalltau, theilslopes
from tqdm.notebook import tqdm

from .mesh_utils import filter_paths, extract_data_for_paths

# ---------------------------------------------------------------------------
# Physical constants (can be overridden as function arguments where relevant)
# ---------------------------------------------------------------------------
T_FREEZE_DEFAULT: float = 273.15   # freezing point of water [K]
T_TOL_DEFAULT: float = 0.5         # half-width of zero-curtain temperature band [K]


# ---------------------------------------------------------------------------
# Statistical helpers
# ---------------------------------------------------------------------------

def theil_sen_fit(
    yvals: np.ndarray,
    xvals: np.ndarray,
) -> tuple[float, float, np.ndarray, float, float]:
    """Fit a Theil–Sen robust linear regression.

    Parameters
    ----------
    yvals : array-like
        Dependent variable.
    xvals : array-like
        Independent variable (e.g. year index).

    Returns
    -------
    slope : float
        Theil–Sen slope.
    intercept : float
        Theil–Sen intercept.
    fit : numpy.ndarray
        Fitted values ``slope * xvals + intercept``.
    tau : float
        Kendall rank correlation coefficient.
    p : float
        Two-sided p-value for the Kendall test.

    Examples
    --------
    >>> slope, intercept, fit, tau, p = theil_sen_fit(zc_days, years)
    """
    yvals = np.asarray(yvals, dtype=float)
    xvals = np.asarray(xvals, dtype=float)
    slope, intercept, _, _ = theilslopes(yvals, xvals)
    fit = slope * xvals + intercept
    tau, p = kendalltau(xvals, yvals)
    return slope, intercept, fit, tau, p


# ---------------------------------------------------------------------------
# Ice table / permafrost helpers
# ---------------------------------------------------------------------------

def determine_ice_table_depth(
    sat_ice_subset: list | np.ndarray,
    closepoly_y_coords: np.ndarray,
    threshold: float = 0.65,
) -> float | None:
    """Find the depth to the ice table in a single vertical column.

    The ice table is defined as the uppermost cell whose ice saturation
    exceeds *threshold*.  Depth is measured downward from the surface
    (top of the supplied coordinate array).

    Parameters
    ----------
    sat_ice_subset : array-like
        Ice saturation values for cells in a column, ordered top-to-bottom
        (same ordering as *closepoly_y_coords*).
    closepoly_y_coords : array-like
        Vertical (z) centroid coordinates [m] for the same cells, ordered
        top-to-bottom.  ``closepoly_y_coords[0]`` is the surface elevation.
    threshold : float, optional
        Ice saturation threshold above which a cell is considered frozen
        (default 0.65).

    Returns
    -------
    depth : float or None
        Depth to the ice table [m], or ``None`` if no cell exceeds the
        threshold (fully thawed column).

    Examples
    --------
    >>> depth = determine_ice_table_depth(sat_ice_col, z_col, threshold=0.65)
    """
    surface = closepoly_y_coords[0]
    for sat, z in zip(sat_ice_subset, closepoly_y_coords):
        if sat > threshold:
            return float(surface - z)
    return None


def build_ice_table_df(
    vis,
    x_min: float = 90.0,
    dx: float = 1.0092592592592524,
    thresh: float = 0.9,
) -> pd.DataFrame:
    """Build a time series of ice-table depth for a single vertical column.

    The column is selected by horizontal position ``[x_min, x_min + 2*dx]``.

    Parameters
    ----------
    vis : ats_xdmf.VisFile
        ATS VisFile object with a loaded mesh.
    x_min : float, optional
        Left edge of the column selection window [m] (default 90).
    dx : float, optional
        Half-width of the column selection window [m].  The right edge is
        ``x_min + 2 * dx``.  Set this to the actual mesh cell width of your
        simulation.
    thresh : float, optional
        Ice saturation threshold for ``determine_ice_table_depth``
        (default 0.9).

    Returns
    -------
    pandas.DataFrame
        Columns:

        * ``'time'`` — simulation time [days]
        * ``'ice_table_depth'`` — depth to ice table [m], or ``NaN`` when
          the column is fully thawed.

    Examples
    --------
    >>> df = build_ice_table_df(vis_1981_2020, x_min=90, dx=1.009, thresh=0.9)
    >>> df.plot(x='time', y='ice_table_depth')
    """
    x_max = x_min + 2 * dx
    poly = vis.getMeshPolygons(cmap="jet_r", linewidth=0)
    paths = poly.get_paths()

    depths, times_list = [], []
    for n in range(len(vis.times)):
        cycle = vis.cycles[n]
        t = vis.times[n]

        selected_paths = filter_paths(paths, x_min, x_max)
        sat_ice = vis.get("saturation_ice", cycle)
        sat_ice_subset = extract_data_for_paths(sat_ice, paths, x_min, x_max)
        closepoly_y = np.array([p.vertices[-1, 1] for p in selected_paths])

        depths.append(determine_ice_table_depth(sat_ice_subset, closepoly_y, threshold=thresh))
        times_list.append(t)

    return pd.DataFrame({"time": times_list, "ice_table_depth": depths})


# ---------------------------------------------------------------------------
# Meteorological aggregation
# ---------------------------------------------------------------------------

def compute_winter_temp(
    df_met: pd.DataFrame,
    years: list[int],
    temp_col: str = "air temperature [C]",
    winter_start: tuple[int, int] = (9, 1),
    winter_end: tuple[int, int] = (5, 31),
) -> pd.DataFrame:
    """Compute mean winter air temperature for each year.

    A "winter" spans from *winter_start* (month, day) of *year* to
    *winter_end* (month, day) of *year + 1*.

    Parameters
    ----------
    df_met : pandas.DataFrame
        Meteorological forcing data with a DatetimeIndex and at least one
        temperature column.
    years : list of int
        Calendar years for which to compute winter averages.  Each entry *y*
        defines the winter *y* → *y+1*.
    temp_col : str, optional
        Column name for air temperature (default ``'air temperature [C]'``).
    winter_start : (month, day), optional
        Start of the winter period (default ``(9, 1)`` → September 1).
    winter_end : (month, day), optional
        End of the winter period (default ``(5, 31)`` → May 31).

    Returns
    -------
    pandas.DataFrame
        Columns ``['year', 'winter_temp']``, one row per entry in *years*.

    Examples
    --------
    >>> df_winter_T = compute_winter_temp(df_met, range(1981, 2024))
    """
    rows = []
    for yr in years:
        start = pd.Timestamp(yr, *winter_start)
        end   = pd.Timestamp(yr + 1, *winter_end)
        mask  = (df_met.index >= start) & (df_met.index <= end)
        rows.append({"year": yr, "winter_temp": df_met.loc[mask, temp_col].mean()})
    return pd.DataFrame(rows)


def compute_winter_snow(
    df_met: pd.DataFrame,
    years: list[int],
    snow_col: str = "Precipitation Snow [m SWE s^-1]",
    dt_seconds: float = 86400.0,
    to_mm: bool = True,
    winter_start: tuple[int, int] = (9, 1),
    winter_end: tuple[int, int] = (5, 31),
) -> pd.DataFrame:
    """Compute cumulative winter snowfall for each year.

    Parameters
    ----------
    df_met : pandas.DataFrame
        Meteorological forcing data with a DatetimeIndex.
    years : list of int
        Calendar years (each entry *y* defines winter *y* → *y+1*).
    snow_col : str, optional
        Column with snowfall rate (default ``'Precipitation Snow [m SWE s^-1]'``).
    dt_seconds : float, optional
        Time step [s] used to convert rate → depth (default 86 400 for daily).
    to_mm : bool, optional
        Multiply by 1000 to convert m SWE → mm SWE (default ``True``).
    winter_start : (month, day), optional
        Start of the winter period (default ``(9, 1)``).
    winter_end : (month, day), optional
        End of the winter period (default ``(5, 31)``).

    Returns
    -------
    pandas.DataFrame
        Columns ``['year', 'snowfall']``, snowfall in mm SWE (or m SWE if
        *to_mm* is ``False``).

    Examples
    --------
    >>> df_snow = compute_winter_snow(df_met, range(1981, 2024))
    """
    factor = dt_seconds * (1000.0 if to_mm else 1.0)
    rows = []
    for yr in years:
        start = pd.Timestamp(yr, *winter_start)
        end   = pd.Timestamp(yr + 1, *winter_end)
        mask  = (df_met.index >= start) & (df_met.index <= end)
        rows.append({"year": yr, "snowfall": (df_met.loc[mask, snow_col] * factor).sum()})
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Time-series extraction from ATS VisFile outputs
# ---------------------------------------------------------------------------

def compute_ponded_depth_series(
    vis_surf,
    epoch: datetime,
    variable: str = "surface-ponded_depth",
) -> pd.DataFrame:
    """Extract the domain-mean ponded depth at every simulation snapshot.

    Parameters
    ----------
    vis_surf : ats_xdmf.VisFile
        ATS surface VisFile object.
    epoch : datetime.datetime
        Reference date corresponding to simulation time 0 (e.g.
        ``datetime(1980, 1, 1)``).
    variable : str, optional
        ATS variable name for ponded depth
        (default ``'surface-ponded_depth'``).

    Returns
    -------
    pandas.DataFrame
        DatetimeIndex, column ``'ponded_m'`` [m].

    Examples
    --------
    >>> epoch = datetime(1980, 1, 1)
    >>> df_pd = compute_ponded_depth_series(vis_surf_1981_2020, epoch)
    >>> df_pd.plot()
    """
    snap_times = np.array(vis_surf.times, dtype=float)
    snap_dates = np.array([epoch + timedelta(days=int(t)) for t in snap_times])
    ponded = np.empty(len(snap_times), dtype=float)

    for n in tqdm(range(len(snap_times)), desc="Ponded depth", leave=False):
        pd_val = vis_surf.get(variable, vis_surf.cycles[n])
        ponded[n] = float(np.mean(pd_val))

    return (
        pd.DataFrame({"date": snap_dates, "ponded_m": ponded})
        .set_index("date")
        .sort_index()
    )


def compute_unfrozen_area_series(
    vis_obj,
    cell_areas: np.ndarray,
    epoch: datetime,
    date_start: datetime,
    date_end: datetime,
    t_freeze: float = T_FREEZE_DEFAULT,
    t_tol: float = T_TOL_DEFAULT,
) -> pd.DataFrame:
    """Compute the total unfrozen subsurface area over a time window.

    A cell is considered *unfrozen* when its temperature exceeds
    ``t_freeze - t_tol``.

    Parameters
    ----------
    vis_obj : ats_xdmf.VisFile
        ATS subsurface VisFile object.
    cell_areas : numpy.ndarray, shape (n_cells,)
        Area of each mesh cell [m²] (from :func:`mesh_utils.get_cell_areas`).
    epoch : datetime.datetime
        Reference date for simulation time = 0.
    date_start, date_end : datetime.datetime
        Start and end of the analysis window (inclusive).
    t_freeze : float, optional
        Freezing point [K] (default 273.15).
    t_tol : float, optional
        Temperature tolerance for the unfrozen/frozen boundary [K]
        (default 0.5).

    Returns
    -------
    pandas.DataFrame
        DatetimeIndex, columns:

        * ``'snap_time_d'`` — snapshot time [days since epoch]
        * ``'unfrozen_area_m2'`` — total unfrozen area [m²]

    Raises
    ------
    RuntimeError
        If no simulation snapshots fall within [*date_start*, *date_end*].

    Examples
    --------
    >>> df_ua = compute_unfrozen_area_series(
    ...     vis_1981_2020, cell_areas, epoch,
    ...     date_start=datetime(1981, 1, 1),
    ...     date_end=datetime(2020, 12, 31),
    ... )
    """
    snap_times = np.array(vis_obj.times, dtype=float)
    snap_dates = np.array([epoch + timedelta(days=int(t)) for t in snap_times])
    mask = (snap_dates >= date_start) & (snap_dates <= date_end)
    idxs = np.nonzero(mask)[0]

    if len(idxs) == 0:
        raise RuntimeError(f"No snapshots between {date_start} and {date_end}.")

    unfrozen = np.empty(len(idxs), dtype=float)
    desc = f"{date_start.year}–{date_end.year}"
    for k, j in enumerate(tqdm(idxs, desc=f"  {desc}", leave=False)):
        T = vis_obj.get("temperature", vis_obj.cycles[j])
        unfrozen[k] = float(cell_areas[T > (t_freeze - t_tol)].sum())

    return (
        pd.DataFrame({
            "date": snap_dates[idxs],
            "snap_time_d": snap_times[idxs],
            "unfrozen_area_m2": unfrozen,
        })
        .set_index("date")
        .sort_index()
    )
