#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:36:13 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""
import pandas as pd
import numpy as np


def xr_dataset_to_timeseries(xr_dataset, variables):
    """Compute timeseries dataframe from xr dataset.

    Parameters
    ----------
    xr_dataset : xarray dataset

    variables : list
        list of varbiale names as string.

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

        # compute how many of the nans are inside aoi (due to snap algorithm)
        out_of_aoi_pixels = (
            len(xr_dataset[var].x) * len(xr_dataset[var].y) - xr_dataset.aoi_pixels
        )
        nans_inside_aoi = nans - out_of_aoi_pixels
        df["aoi_nan_percentage"] = nans_inside_aoi / xr_dataset.aoi_pixels

        df[var + "_se"] = df[var + "_std"] / np.sqrt(sample_n)

    return df