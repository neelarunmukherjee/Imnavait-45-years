"""
mesh_utils.py — mesh geometry helpers for ATS polygon outputs.

Typical usage
-------------
>>> from utils.mesh_utils import get_cell_areas, filter_paths, extract_data_for_paths
>>> areas = get_cell_areas(vis_obj)
>>> paths_sub = filter_paths(all_paths, x_min=90, x_max=92)
>>> sat_sub   = extract_data_for_paths(sat_ice, all_paths, x_min=90, x_max=92)
"""

import numpy as np


def get_cell_areas(vis_obj) -> np.ndarray:
    """Compute the area of every mesh cell using the shoelace formula.

    Parameters
    ----------
    vis_obj : ats_xdmf.VisFile
        An ATS VisFile object with a loaded mesh.

    Returns
    -------
    areas : numpy.ndarray, shape (n_cells,)
        Cell areas [m²], one per polygon path returned by
        ``vis_obj.getMeshPolygons()``.

    Raises
    ------
    AssertionError
        If any computed area is non-positive (degenerate polygon).

    Notes
    -----
    The shoelace formula gives the signed area of a simple polygon; taking the
    absolute value makes the result orientation-independent.
    """
    paths = vis_obj.getMeshPolygons(cmap="turbo", linewidth=0).get_paths()
    areas = np.asarray(
        [
            0.5
            * np.abs(
                np.dot(p.vertices[:, 0], np.roll(p.vertices[:, 1], -1))
                - np.dot(p.vertices[:, 1], np.roll(p.vertices[:, 0], -1))
            )
            for p in paths
        ],
        dtype=float,
    )
    assert np.all(areas > 0), "Non-positive cell area detected — check mesh."
    return areas


def filter_paths(paths: list, x_min: float, x_max: float) -> list:
    """Return only those mesh polygon paths whose vertices all lie within [x_min, x_max].

    Parameters
    ----------
    paths : list of matplotlib.path.Path
        Full list of polygon paths from ``PolyCollection.get_paths()``.
    x_min, x_max : float
        Horizontal extent filter [m].

    Returns
    -------
    list of matplotlib.path.Path
        Subset of *paths* where every vertex x-coordinate satisfies
        ``x_min <= x <= x_max``.

    Examples
    --------
    >>> poly   = vis.getMeshPolygons(cmap='jet_r', linewidth=0)
    >>> paths  = poly.get_paths()
    >>> subset = filter_paths(paths, x_min=90, x_max=92)
    """
    return [
        p
        for p in paths
        if np.all((p.vertices[:, 0] >= x_min) & (p.vertices[:, 0] <= x_max))
    ]


def extract_data_for_paths(
    data: np.ndarray, paths: list, x_min: float, x_max: float
) -> list:
    """Extract scalar data values for cells whose polygons fall within [x_min, x_max].

    Parameters
    ----------
    data : array-like, shape (n_cells,)
        One scalar value per mesh cell (e.g. ice saturation, temperature).
    paths : list of matplotlib.path.Path
        Full list of polygon paths (same ordering as *data*).
    x_min, x_max : float
        Horizontal extent filter [m].

    Returns
    -------
    list
        Data values corresponding to the filtered cells, in the same order as
        the filtered paths.

    Examples
    --------
    >>> sat_ice_subset = extract_data_for_paths(sat_ice, paths, 90, 92)
    """
    return [
        data[i]
        for i, p in enumerate(paths)
        if np.all((p.vertices[:, 0] >= x_min) & (p.vertices[:, 0] <= x_max))
    ]
