"""
plotting.py — visualisation helpers for ATS permafrost simulation outputs.

Typical usage
-------------
>>> from utils.plotting import plot_zc_duration
>>> fig, zc_counts = plot_zc_duration(
...     vis=vis_obj,
...     vis_surf=vis_surf_obj,
...     epoch=datetime(1980, 1, 1),
...     target_year=2017,
...     elev=surface_elevation,
... )
"""

from __future__ import annotations

from datetime import datetime

import matplotlib.pyplot as plt
import numpy as np

# Default physical constants
T_FREEZE_DEFAULT: float = 273.15   # [K]
T_TOL_DEFAULT: float = 0.5         # [K]


def plot_zc_duration(
    vis,
    vis_surf,
    epoch: datetime,
    target_year: int,
    elev: np.ndarray,
    *,
    t_freeze: float = T_FREEZE_DEFAULT,
    t_tol: float = T_TOL_DEFAULT,
    winter_start: tuple[int, int] = (9, 1),
    winter_end: tuple[int, int] = (6, 30),
    zlim: tuple[float, float] = (888.2, 892.5),
    clim: tuple[float, float] = (0, 200),
    show_plot: bool = True,
    figsize: tuple[float, float] = (10, 4),
) -> tuple[plt.Figure | None, np.ndarray]:
    """Compute and optionally plot a zero-curtain duration map.

    The zero curtain is defined as the period when subsurface temperature
    stays within *t_tol* of *t_freeze*.  The map shows how many simulation
    snapshots each cell spends in the zero curtain during the winter season
    ``target_year/winter_start … (target_year+1)/winter_end``.

    Parameters
    ----------
    vis : ats_xdmf.VisFile
        ATS subsurface VisFile object for the period that spans the requested
        winter.  The caller is responsible for selecting the correct run.
    vis_surf : ats_xdmf.VisFile
        Matching surface VisFile (used for the elevation line on the plot).
    epoch : datetime.datetime
        Reference date corresponding to simulation time = 0.
    target_year : int
        Calendar year at the start of the winter (e.g. 2017 for the
        2017–2018 winter).
    elev : array-like
        Surface elevation values [m] aligned with ``vis_surf.centroids[:,0]``
        (plotted as a blue line).
    t_freeze : float, optional
        Freezing point [K] (default 273.15).
    t_tol : float, optional
        Half-width of the zero-curtain temperature band [K] (default 0.5).
    winter_start : (month, day), optional
        Start of the winter window (default ``(9, 1)`` → Sep 1).
    winter_end : (month, day), optional
        End of the winter window (default ``(6, 30)`` → Jun 30).
    zlim : (float, float), optional
        Vertical axis limits for the cross-section plot [m] (default
        ``(888.2, 892.5)``).  Adjust to your domain.
    clim : (float, float), optional
        Colorbar limits (min, max) in days (default ``(0, 200)``).
    show_plot : bool, optional
        If ``False`` skip the figure and return ``(None, zc_counts)``
        (default ``True``).
    figsize : (width, height), optional
        Figure size in inches (default ``(10, 4)``).

    Returns
    -------
    fig : matplotlib.figure.Figure or None
        The figure object, or ``None`` when *show_plot* is ``False``.
    zc_counts : numpy.ndarray, shape (n_cells,)
        Number of snapshots each cell spent in the zero curtain.

    Examples
    --------
    Minimal usage with a single simulation run:

    >>> fig, counts = plot_zc_duration(
    ...     vis=vis_1981_2020,
    ...     vis_surf=vis_surf_1981_2020,
    ...     epoch=datetime(1980, 1, 1),
    ...     target_year=2000,
    ...     elev=surface_elev,
    ... )

    Count-only (no plot):

    >>> _, counts = plot_zc_duration(..., show_plot=False)
    """
    zc_start = datetime(target_year,     *winter_start)
    zc_end   = datetime(target_year + 1, *winter_end)

    day_start = (zc_start - epoch).days
    day_end   = (zc_end   - epoch).days

    idx_start = int(np.argmin(np.abs(vis.times - day_start)))
    idx_end   = int(np.argmin(np.abs(vis.times - day_end)))
    idxs      = range(idx_start, idx_end + 1)

    poly0    = vis.getMeshPolygons(cmap="turbo", linewidth=0)
    ncells   = len(poly0.get_paths())
    zc_counts = np.zeros(ncells, dtype=np.int32)

    for n in idxs:
        T    = vis.get("temperature", vis.cycles[n])
        mask = np.abs(T - t_freeze) <= t_tol
        zc_counts += mask.astype(np.int32)

    if not show_plot:
        return None, zc_counts

    plt.rcParams["font.size"] = 15
    fig, ax = plt.subplots(figsize=figsize)

    poly = vis.getMeshPolygons(cmap="turbo", linewidth=0)
    poly.set_array(zc_counts.astype(float))
    poly.set_clim(*clim)

    paths = poly.get_paths()
    xs = np.concatenate([p.vertices[:, 0] for p in paths])
    ax.add_collection(poly)
    ax.set_xlim(xs.min() - 1, xs.max() + 1)
    ax.set_ylim(*zlim)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_title(
        f"Zero-curtain duration — {target_year}–{target_year + 1}", fontsize=14
    )
    plt.colorbar(poly, ax=ax).set_label("ZC duration [days]")
    ax.plot(
        vis_surf.centroids[:, 0],
        elev,
        "b",
        linewidth=1,
        label="Elevation",
    )
    plt.tight_layout()
    plt.show()

    return fig, zc_counts
