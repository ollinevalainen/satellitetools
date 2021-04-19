#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:36:13 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""
import pandas as pd
import numpy as np
from satellitetools.biophys import SNAP_BIO_RMSE


def xr_dataset_to_timeseries(xr_dataset, variables, add_uncertainty=False):
    """Compute timeseries dataframe from xr dataset.

    Parameters
    ----------
    xr_dataset : xarray dataset

    variables : list
        list of varbiale names as string.

    add_uncertainty : bool, default False
        Adds variable {variable}_uncertainty to dataframe. Currently,
        uncertainty is equal to standar error (se) or if variable is biophysical
        variable from biophys_xarray, it sqrt(se^2 + RMSE^2) where RMSE is
        the variable RMSE from the SNAP biophysical processor developers
        (see biophys_xarray.py and linked ATBD).

    Returns
    -------
    df : pandas dataframe
        Pandas dataframe with mean, std, se and percentage of NaNs inside AOI.

    """
    df = pd.DataFrame({"Date": pd.to_datetime(xr_dataset.time.values)})

    for var in variables:
        df[var] = xr_dataset[var].mean(dim=["x", "y"])
        df[var + "_std"] = xr_dataset[var].std(dim=["x", "y"])

        # nans occure due to missging data from 1D to 2D array
        # (pixels outside the polygon),
        # from snap algorihtm nans occure due to input/output ouf of bounds
        # checking.
        # TODO: flaggging with snap biophys algorith or some other solution to
        # check which nan are from snap    algorithm and which from 1d to 2d transformation
        nans = np.isnan(xr_dataset[var]).sum(dim=["x", "y"])
        sample_n = len(xr_dataset[var].x) * len(xr_dataset[var].y) - nans

        df[var + "_se"] = df[var + "_std"] / np.sqrt(sample_n)

        if add_uncertainty:
            df = compute_uncertainty(df, var)

        if hasattr(xr_dataset, "aoi_pixels"):
            # compute how many of the nans are inside aoi (due to snap algorithm)
            out_of_aoi_pixels = (
                len(xr_dataset[var].x) * len(xr_dataset[var].y) - xr_dataset.aoi_pixels
            )
            nans_inside_aoi = nans - out_of_aoi_pixels
            df["aoi_nan_percentage"] = nans_inside_aoi / xr_dataset.aoi_pixels

    return df


def compute_uncertainty(df, var):
    if var in SNAP_BIO_RMSE.keys():
        df[var + "_uncertainty"] = np.sqrt(
            df[var + "_se"] ** 2 + SNAP_BIO_RMSE[var] ** 2
        )
    else:
        df[var + "_uncertainty"] = df[var + "_se"]
    return df
