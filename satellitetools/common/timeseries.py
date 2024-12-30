#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for timeseries data handling.

@author: Olli Nevalainen (Finnish Meteorological Institute)

"""
import logging
from typing import List, Union

try:
    # breaking change introduced in python 3.11
    from enum import StrEnum
except ImportError:
    from enum import Enum

    class StrEnum(str, Enum):
        pass


import numpy as np
import pandas as pd
import xarray as xr

from satellitetools.biophys.biophys import (
    SNAP_BIO_RMSE,
    BiophysVariable,
    VegetationIndex,
)

logger = logging.getLogger(__name__)


class ConfidenceLevel(StrEnum):
    """Confidence level for calculating the confidence interval bounds."""

    C90 = "90"
    C95 = "95"
    C99 = "99"


def xr_dataset_to_timeseries(
    xr_dataset: xr.Dataset,
    variables: List[Union[BiophysVariable, VegetationIndex]],
    add_uncertainty: bool = False,
    add_confidence_intervals: bool = False,
    confidence_level: ConfidenceLevel = ConfidenceLevel.C95,
) -> pd.DataFrame:
    """Compute timeseries dataframe from xarray dataset.

    Parameters
    ----------
    xr_dataset : xr.Dataset
        xarray dataset with Sentinel-2 data.

    variables : List[BiophysVariable]
        List of variables to compute timeseries for.

    add_uncertainty : bool, default False
        Adds variable {variable}_uncertainty and confidence intervals to dataframe.
        Currently, uncertainty is equal to standar error (se) or if variable is
        biophysical variable from biophys_xarray, it sqrt(se^2 + RMSE_mean^2) where
        RMSE_mean is propagated uncertainty for the individual observations/pixels
        uncertainties. Uncertainty for the individual pixels is considered to be the
        variable RMSE from the SNAP biophysical processor developers
        (see biophys_xarray.py and linked ATBD) (i.e. same for all pixels).

    confidence_level : ConfidenceLevel, default ConfidenceLevel.C95
        Confidence level (%) for calculating the confidence interval bounds.
        Options "90", "95" & "99"

    Returns
    -------
    df : pandas dataframe
        Pandas dataframe with mean, std, se and percentage of NaNs inside AOI.

    """
    df = pd.DataFrame(index=pd.to_datetime(xr_dataset.time.values))

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
        if var in SNAP_BIO_RMSE:
            # if data is upsampled, take this into account in
            # uncertainty (n "artificially increased")
            sample_n = _adjust_sample_size(xr_dataset, sample_n)

        df[var] = xr_dataset[var].mean(dim=["x", "y"])
        df[var + "_F050"] = xr_dataset[var].median(dim=["x", "y"])

        df[var + "_std"] = xr_dataset[var].std(dim=["x", "y"])

        df[var + "_se"] = df[var + "_std"] / np.sqrt(sample_n)

        if add_uncertainty or add_confidence_intervals:
            df = compute_uncertainty(df, xr_dataset, var, sample_n)

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
            df[var + "_aoi_nan_percentage"] = nans_inside_aoi / xr_dataset.aoi_pixels

    return df


def _adjust_sample_size(xr_dataset: xr.Dataset, sample_n: float) -> float:
    """Get sample size for uncertainty calculation.

    Parameters
    ----------
    xr_dataset : xr.Dataset
        xarray dataset with Sentinel-2 data.
    sample_n : float
        Original sample size.

    Returns
    -------
    sample_n : float
        Adjusted sample size.

    """
    # 20 = 20 m which is the SNAP_BIO function standard resolution
    try:
        resampling_ratio = 20 / np.abs(xr_dataset.x[1] - xr_dataset.x[0])
    except IndexError:
        try:
            resampling_ratio = 20 / np.abs(xr_dataset.y[1] - xr_dataset.y[0])
        except IndexError:
            # Assume that the data is not resampled
            logger.warning(
                "Could not determine resampling ratio, assuming no resampling"
            )
            resampling_ratio = 1

    pixel_multiplier = (
        resampling_ratio**2
    )  # doubling resolution quadruples amount of pixels etc.
    if resampling_ratio > 1:
        sample_n = sample_n / pixel_multiplier

        sample_n[sample_n < 1] = (
            1  # with small areas the adjusted sample_n might be < 1
        )
    return sample_n


def propagate_rmse(n: int, rmse: float) -> float:
    """Propagate RMSE for the mean value.

    Parameters
    ----------
    n : int
        Sample size.
    rmse : float
        Root mean square error.

    Returns
    -------
    propagated_rmse : float
        Propagated RMSE for the mean value.

    """

    propagated_rmse = np.sqrt(np.sum([rmse**2] * int(n))) / n
    return propagated_rmse


def compute_uncertainty(
    df: pd.DataFrame, xr_dataset: xr.Dataset, var: BiophysVariable, sample_n: float
) -> pd.DataFrame:
    """Compute uncertainty for the variable.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with the variable.
    xr_dataset : xr.Dataset
        xarray dataset with Sentinel-2 data.
    var : BiophysVariable
        Biophysical variable.
    sample_n : float
        Sample size for uncertainty calculation.
    Returns
    -------
    df : pd.DataFrame
        Dataframe with the variable and uncertainty.

    """
    if var in SNAP_BIO_RMSE:
        sample_n = _adjust_sample_size(xr_dataset, sample_n)

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


def compute_confidence_intervals(
    df: pd.DataFrame, var: BiophysVariable, confidence_level=ConfidenceLevel.C95
) -> pd.DataFrame:
    """Compute confidence intervals for the variable.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with the variable.
    var : BiophysVariable
        Biophysical variable.
    confidence_level : ConfidenceLevel, default ConfidenceLevel.C95
        Confidence level (%) for calculating the confidence interval bounds.
        Options "90", "95" & "99"

    Returns
    -------
    df : pd.DataFrame
        Dataframe with the variable and confidence intervals.

    """
    if confidence_level == ConfidenceLevel.C90:
        z_score = 1.645
        ci_min = "_F005"
        ci_max = "_F095"
    elif confidence_level == ConfidenceLevel.C95:
        z_score = 1.96
        ci_min = "_F0025"
        ci_max = "_F0975"
    elif confidence_level == ConfidenceLevel.C99:
        z_score = 2.576
        ci_min = "_F0005"
        ci_max = "_F0995"
    else:
        logger.error("Unknown confidence level")
        raise ValueError("Unknown confidence level")

    # uncertainty to confidence intervals
    df[var + ci_min] = df[var] - z_score * df[var + "_uncertainty"]
    df[var + ci_max] = df[var] + z_score * df[var + "_uncertainty"]

    if var in SNAP_BIO_RMSE:
        # Cap unrealistic lower bound negative values to 0
        df.loc[df[var + ci_min] < 0, var + ci_min] = 0

    return df
