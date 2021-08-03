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
        )

        masked_kwds = kwds
        masked_kwds.update(
            transform=masked_transform,
            height=masked_data.shape[-2],
            width=masked_data.shape[-1],
            nodata=no_data,
        )
    return masked_data, masked_kwds


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
        new_transform = rasterio.transform.Affine(
            dataset.profile["transform"].a / scale_factor,
            dataset.profile["transform"].b,
            dataset.profile["transform"].c,
            dataset.profile["transform"].d,
            dataset.profile["transform"].e / scale_factor,
            dataset.profile["transform"].f,
        )

        # new_transform = dataset.transform * dataset.transform.scale(
        #     (dataset.width / data.shape[-1]), (dataset.height / data.shape[-2])
        # )
        new_kwds = dataset.profile
        new_kwds.update(
            transform=new_transform, driver="GTiff", height=new_height, width=new_width
        )
    return data, new_kwds
