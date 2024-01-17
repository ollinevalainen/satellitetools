#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:36:13 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""
import sys

import numpy as np
import pandas as pd
import xarray as xr

from satellitetools.biophys.biophys import SNAP_BIO_RMSE


def xr_dataset_to_timeseries(
    xr_dataset,
    variables,
    add_uncertainty=False,
    add_confidence_intervals=False,
    confidence_level="95",
):
    """Compute timeseries dataframe from xr dataset.

    Parameters
    ----------
    xr_dataset : xarray dataset

    variables : list
        list of varbiale names as string.

    add_uncertainty : bool, default False
        Adds variable {variable}_uncertainty and confidence intervals to dataframe.
        Currently, uncertainty is equal to standar error (se) or if variable is
        biophysical variable from biophys_xarray, it sqrt(se^2 + RMSE_mean^2) where
        RMSE_mean is propagated uncertainty for the individual observations/pixels
        uncertainties. Uncertainty for the individual pixels is considered to be the
        variable RMSE from the SNAP biophysical processor developers
        (see biophys_xarray.py and linked ATBD) (i.e. same for all pixels).

    confidence_level : str, default "95"
        Confidence level (%) for calculating the confidence interval bounds.
        Options "90", "95" & "99"

    Returns
    -------
    df : pandas dataframe
        Pandas dataframe with mean, std, se and percentage of NaNs inside AOI.

    """
    df = pd.DataFrame({"Date": pd.to_datetime(xr_dataset.time.values)})

    for var in variables:
        # nans occure due to missging data from 1D to 2D array
        # (pixels outside the polygon),
        # from snap algorihtm nans occure due to input/output ouf of bounds
        # checking.
        # TODO: flaggging with snap biophys algorith or some other solution to
        # check which nan are from snap algorithm and which from 1d to 2d transformation
        nans = np.isnan(xr_dataset[var]).sum(dim=["x", "y"])
        sample_n = len(xr_dataset[var].x) * len(xr_dataset[var].y) - nans

        # Filter dates with all pixels equal to nan (i.e. sample n = 0)
        non_zero_n = np.where(sample_n != 0)[0]
        xr_dataset = xr_dataset.isel(time=non_zero_n)
        sample_n = sample_n.isel(time=non_zero_n)
        nans = nans.isel(time=non_zero_n)
        # Drop also from df
        df = df.iloc[non_zero_n]

        # adjust n if data is resampled
        if var in SNAP_BIO_RMSE.keys():
            # if data is upsampled, take this into account in
            # uncertainty (n "artificially increased")
            # 20 = 20 m which is the SNAP_BIO function standard resolution
            try:
                resampling_ratio = 20 / np.abs(xr_dataset.x[1] - xr_dataset.x[0])
            except IndexError:
                try:
                    resampling_ratio = 20 / np.abs(xr_dataset.y[1] - xr_dataset.y[0])
                except IndexError:
                    # Print assumes that the data is not resampled
                    resampling_ratio = 1

            pixel_multiplier = (
                resampling_ratio**2
            )  # doubling resolution quadruples amount of pixels etc.
            if resampling_ratio > 1:
                sample_n = sample_n / pixel_multiplier

                sample_n[
                    sample_n < 1
                ] = 1  # with small areas the adjusted sample_n might be < 1

        df[var] = xr_dataset[var].mean(dim=["x", "y"])
        df[var + "_F050"] = xr_dataset[var].median(dim=["x", "y"])

        df[var + "_std"] = xr_dataset[var].std(dim=["x", "y"])

        df[var + "_se"] = df[var + "_std"] / np.sqrt(sample_n)

        if add_uncertainty or add_confidence_intervals:
            df = compute_uncertainty(df, xr_dataset, var)

            if add_confidence_intervals:
                df = compute_confidence_intervals(
                    df, var, confidence_level=confidence_level
                )

        if hasattr(xr_dataset, "aoi_pixels"):
            # compute how many of the nans are inside aoi (due to snap algorithm)
            out_of_aoi_pixels = (
                len(xr_dataset[var].x) * len(xr_dataset[var].y) - xr_dataset.aoi_pixels
            )
            nans_inside_aoi = nans - out_of_aoi_pixels
            df["aoi_nan_percentage"] = nans_inside_aoi / xr_dataset.aoi_pixels

    return df


def propagate_rmse(n, rmse):
    propagated_rmse = np.sqrt(np.sum([rmse**2] * int(n))) / n
    return propagated_rmse


def compute_uncertainty(df, xr_dataset, var):
    if var in SNAP_BIO_RMSE.keys():
        nans = np.isnan(xr_dataset[var]).sum(dim=["x", "y"])
        sample_n = len(xr_dataset[var].x) * len(xr_dataset[var].y) - nans

        # adjust n if data is resampled
        if var in SNAP_BIO_RMSE.keys():
            # if data is upsampled, take this into account in
            # uncertainty (n "artificially increased")
            # 20 = 20 m which is the SNAP_BIO function standard resolution
            try:
                resampling_ratio = 20 / np.abs(xr_dataset.x[1] - xr_dataset.x[0])
            except IndexError:
                try:
                    resampling_ratio = 20 / np.abs(xr_dataset.y[1] - xr_dataset.y[0])
                except IndexError:
                    # Print assumes that the data is not resampled
                    resampling_ratio = 1

            pixel_multiplier = (
                resampling_ratio**2
            )  # doubling resolution quadruples amount of pixels etc.
            if resampling_ratio > 1:
                sample_n = sample_n / pixel_multiplier

                sample_n[
                    sample_n < 1
                ] = 1  # with small areas the adjusted sample_n might be < 1

        # propagated RMSE for the
        # mean value (RMSE as uncertainty for individual observations/pixels)
        rmse_mean = xr.apply_ufunc(
            propagate_rmse,
            sample_n,
            kwargs={"rmse": SNAP_BIO_RMSE[var]},
            vectorize=True,
        )
        # rmse_means = np.sqrt(np.sum([SNAP_BIO_RMSE[var] ** 2 ] * n)) / n

        # if data is upsampled, take this into account in
        # uncertainty (n "artificially increased")
        # 20 = 20 m which is the SNAP_BIO function standard resolution
        # resampling_ratio = 20 / np.abs(xr_dataset.x[1] - xr_dataset.x[0])
        # pixel_multiplier = (
        #     resampling_ratio ** 2
        # )  # doubling resolution quadruples amount of pixels etc.
        # if resampling_ratio > 1:
        #     rmse_mean = rmse_mean * pixel_multiplier / np.sqrt(pixel_multiplier)

        df[var + "_uncertainty"] = np.sqrt(df[var + "_std"] ** 2 + rmse_mean**2)

    else:
        df[var + "_uncertainty"] = df[var + "_std"]

    return df


def compute_confidence_intervals(df, var, confidence_level="95"):
    if confidence_level == "90":
        z_score = 1.645
        ci_min = "_F005"
        ci_max = "_F095"
    elif confidence_level == "95":
        z_score = 1.96
        ci_min = "_F0025"
        ci_max = "_F0975"
    elif confidence_level == "99":
        z_score = 2.576
        ci_min = "_F0005"
        ci_max = "_F0995"
    else:
        sys.exit("Unknown confidence level")

    # uncertainty to confidence intervals
    df[var + ci_min] = df[var] - z_score * df[var + "_uncertainty"]
    df[var + ci_max] = df[var] + z_score * df[var + "_uncertainty"]

    if var in SNAP_BIO_RMSE.keys():
        # Cap unrealistic lower bound negative values to 0
        df.loc[df[var + ci_min] < 0, var + ci_min] = 0

    return df
