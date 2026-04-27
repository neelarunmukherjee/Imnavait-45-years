"""
data_io.py — I/O helpers for ATS permafrost simulations.

Typical usage
-------------
>>> from utils.data_io import plot_h5_data, get_temp_profile
>>> stats, data, t = plot_h5_data("forcing.h5")
>>> times, z, T   = get_temp_profile("/path/to/ats/run/")
"""

import numpy as np
import pandas as pd
import h5py
import ats_xdmf  # ATS post-processing library


def plot_h5_data(filename: str):
    """Read an HDF5 forcing file and return basic statistics.

    Parameters
    ----------
    filename : str
        Path to the HDF5 file. The file may contain a ``'time [s]'`` dataset
        and any number of 1-D data vectors.

    Returns
    -------
    stats_df : pandas.DataFrame
        Columns: ``['Dataset Keys', 'Min Value', 'Max Value']``.
    data_vectors : dict[str, numpy.ndarray]
        Mapping of dataset name → 1-D array.
    time_vector : numpy.ndarray or None
        The ``'time [s]'`` array, or ``None`` if not present.

    Examples
    --------
    >>> stats, data, t = plot_h5_data("met_forcing.h5")
    >>> print(stats)
    """
    with h5py.File(filename, "r") as f:
        keys = [k for k in f.keys() if k != "time [s]"]
        time_vector = f["time [s]"][:] if "time [s]" in f else None
        data_vectors = {k: f[k][:] for k in keys if f[k].ndim == 1}

    min_values = [np.min(v) for v in data_vectors.values()]
    max_values = [np.max(v) for v in data_vectors.values()]
    stats_df = pd.DataFrame({
        "Dataset Keys": list(data_vectors.keys()),
        "Min Value": min_values,
        "Max Value": max_values,
    })
    return stats_df, data_vectors, time_vector


def get_temp_profile(directory: str, time_unit: str = "d"):
    """Load the subsurface temperature field from an ATS VisFile.

    Parameters
    ----------
    directory : str
        Path to the ATS run directory containing ``visdump_data.h5`` (or
        equivalent XDMF output).
    time_unit : str, optional
        Time unit to pass to ``ats_xdmf.VisFile`` (default ``'d'`` for days).

    Returns
    -------
    times : numpy.ndarray
        Simulation times in *time_unit*.
    z : numpy.ndarray
        Vertical centroid coordinates [m] of each subsurface cell.
    temperature : numpy.ndarray, shape (n_times, n_cells)
        Temperature field [K].

    Examples
    --------
    >>> times, z, T = get_temp_profile("/data/run_1981_2020/")
    >>> import matplotlib.pyplot as plt
    >>> plt.pcolormesh(times, z, T.T)
    """
    vis = ats_xdmf.VisFile(directory, time_unit=time_unit)
    vis.loadMesh(columnar=True)
    temperature = np.array(vis.getArray("temperature"))
    z = vis.centroids[:, 2]
    times = vis.times
    return times, z, temperature
