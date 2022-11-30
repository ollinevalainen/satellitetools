#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for raster data handling.
TODO: UPDATE docstrings!!

@author: Olli Nevalainen (Finnish Meteorological Institute)
Created on Tue Mar 16 10:45:05 2021
"""
import rasterio
from rasterio import mask
from rasterio.enums import Resampling
from rasterio.windows import get_data_window
import numpy as np


def mask_raster(raster, aoi_geometry, no_data=0):
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
            # get_data_window does not work with np.nan, change temporally
            NODATA = -99999
            masked_data[np.isnan(masked_data)] = NODATA
        else:
            NODATA = no_data

        # crop nodata row and columns
        crop_window = get_data_window(masked_data, nodata=NODATA)

        cropped_data = masked_data[
            0,
            crop_window.row_off : crop_window.row_off + crop_window.height,
            crop_window.col_off : crop_window.col_off + crop_window.width,
        ]

        # change nodata to user-defined
        cropped_data[cropped_data == NODATA] = no_data

        cropped_kwds = masked_kwds.copy()
        cropped_kwds.update(
            height=crop_window.height,
            width=crop_window.width,
            nodata=no_data,
            transform=rasterio.windows.transform(crop_window, masked_kwds["transform"]),
        )
    return cropped_data, cropped_kwds


def resample_raster(raster, target_gsd, target_height=None, target_width=None):
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
            resampling=Resampling.bilinear,
        )

        # scale image transform
        new_transform = dataset.transform * dataset.transform.scale(
            (dataset.width / data.shape[-1]), (dataset.height / data.shape[-2])
        )

        new_kwds = dataset.profile
        new_kwds.update(
            transform=new_transform, driver="GTiff", height=new_height, width=new_width
        )
    return data, new_kwds
