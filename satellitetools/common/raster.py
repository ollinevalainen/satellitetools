#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for raster data handling.

@author: Olli Nevalainen (Finnish Meteorological Institute)

"""
import os
from typing import Optional, Tuple, Union

import numpy as np
import rasterio
from rasterio import mask
from rasterio.enums import Resampling
from rasterio.io import MemoryFile
from rasterio.windows import get_data_window
from shapely.geometry import Polygon


def mask_raster(
    raster: Union[MemoryFile, os.PathLike],
    aoi_geometry: Polygon,
    no_data: Union[float, int],
) -> Tuple[np.ndarray, dict]:
    """Mask raster data with area of interest geometry.

    Parameters:
    ----------------
    raster: Union[MemoryFile, os.PathLike]
        Raster data.
    aoi_geometry: Polygon
        Area of interest geometry.
    no_data: Union[float, int]
        No data value.

    Returns:
    ----------------
    Tuple[np.ndarray, dict]
        Tuple of masked data and metadata.
    """

    with raster.open() as dataset:
        kwds = dataset.profile
        masked_data, masked_transform = mask.mask(
            dataset,
            [aoi_geometry],
            all_touched=False,
            invert=False,
            nodata=no_data,
            filled=True,
            crop=True,  # default False
            pad=False,
        )

        masked_kwds = kwds
        masked_kwds.update(
            transform=masked_transform,
            height=masked_data.shape[-2],
            width=masked_data.shape[-1],
            nodata=no_data,
        )
        if np.isnan(no_data):
            # get_data_window does not work with np.nan, change temporary to -99999
            tmp_no_data = -99999
            masked_data[np.isnan(masked_data)] = tmp_no_data
        else:
            tmp_no_data = no_data

        # crop nodata row and columns
        crop_window = get_data_window(masked_data, nodata=tmp_no_data)

        cropped_data = masked_data[
            0,
            crop_window.row_off : crop_window.row_off + crop_window.height,
            crop_window.col_off : crop_window.col_off + crop_window.width,
        ]

        # change nodata to user-defined
        cropped_data[cropped_data == tmp_no_data] = no_data

        cropped_kwds = masked_kwds.copy()
        cropped_kwds.update(
            height=crop_window.height,
            width=crop_window.width,
            nodata=no_data,
            transform=rasterio.windows.transform(crop_window, masked_kwds["transform"]),
        )
    return cropped_data, cropped_kwds


def resample_raster(
    raster: Union[MemoryFile, os.PathLike],
    target_gsd: float,
    resampling_method: Resampling,
    target_height: Optional[int] = None,
    target_width: Optional[int] = None,
) -> Tuple[np.ndarray, dict]:
    """Resample raster data to target ground sample distance.

    Parameters:
    ----------------
    raster: Union[MemoryFile, os.PathLike]
        Raster data.
    target_gsd: float
        Target ground sample distance.
    resampling_method: Resampling
        Resampling method.
    target_height: Optional[int]
        Target height.
    target_width: Optional[int]
        Target width.

    Returns:
    ----------------
    Tuple[np.ndarray, dict]
        Tuple of resampled data and metadata.
    """
    with raster.open() as dataset:
        # a first parameter in Affine i.e. the pixel x size
        scale_factor = dataset.profile["transform"].a / target_gsd
        # resample data to target shape
        if target_height:
            new_height = target_height
        else:
            new_height = round(dataset.height * scale_factor)

        if target_width:
            new_width = target_width
        else:
            new_width = round(dataset.width * scale_factor)

        data = dataset.read(
            out_shape=(dataset.count, new_height, new_width),
            resampling=resampling_method,
        )

        # scale image transform
        # new_transform = dataset.transform * dataset.transform.scale(
        #     (dataset.width / data.shape[-1]), (dataset.height / data.shape[-2])
        # )
        new_transform = dataset.transform * dataset.transform.scale(1 / scale_factor)

        new_kwds = dataset.profile
        new_kwds.update(
            transform=new_transform,
            driver="GTiff",
            height=data.shape[-2],
            width=data.shape[-1],
        )
    return data, new_kwds
