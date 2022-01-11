#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module to get Sentinel-2 data from AWS Open data registry,
where Sentinel-2 (level 2A) data is available as cloud-optimized
geotiffs (https://registry.opendata.aws/sentinel-2-l2a-cogs/).


@author: Olli Nevalainen (Finnish Meteorological Institute)
Created on Fri Mar 12 15:47:31 2021
"""
import sys
import satsearch
import rasterio
import numpy as np
import pandas as pd
import xarray as xr
import urllib
import datetime
from collections import Counter

import xmltodict
from rasterio import MemoryFile
from satellitetools.common.vector import (
    transform_crs,
    expand_bounds,
    create_coordinate_arrays,
)
from satellitetools.common.raster import mask_raster, resample_raster
from satellitetools.common.sentinel2 import (
    S2_SCL_CLASSES,
    S2_REFL_TRANS,
    S2_FILTER1,
    filter_s2_qi_dataframe,
)


def search_s2_cogs(aoi, req_params):
    print(
        "Searching S2 data from {} to {} for area {}".format(
            req_params.datestart, req_params.dateend, aoi.name
        )
    )
    bbox = list(aoi.geometry.bounds)
    dates = "{}/{}".format(req_params.datestart, req_params.dateend)
    URL = "https://earth-search.aws.element84.com/v0"
    search = satsearch.Search(
        url=URL, collections=["sentinel-s2-l2a-cogs"], datetime=dates, bbox=bbox
    )
    if search.found() == 0:
        print("No available data for specified time!")
        items = None
    else:
        items = search.items()
        print("Found {} available S2 acquisition dates.".format(len(items)))
    return items


def get_xml_metadata(item):
    with urllib.request.urlopen(item.assets["metadata"]["href"]) as url:
        metadata = xmltodict.parse(url.read().decode())
    metadata = metadata.popitem(last=False)[1]
    return metadata


def get_angles(item):
    metadata = get_xml_metadata(item)
    sunzen = np.float64(
        metadata["n1:Geometric_Info"]["Tile_Angles"]["Mean_Sun_Angle"]["ZENITH_ANGLE"][
            "#text"
        ]
    )
    sunaz = np.float64(
        metadata["n1:Geometric_Info"]["Tile_Angles"]["Mean_Sun_Angle"]["AZIMUTH_ANGLE"][
            "#text"
        ]
    )

    viewangles = metadata["n1:Geometric_Info"]["Tile_Angles"][
        "Mean_Viewing_Incidence_Angle_List"
    ]["Mean_Viewing_Incidence_Angle"]

    viewazs = np.array([np.float64(d["AZIMUTH_ANGLE"]["#text"]) for d in viewangles])
    viewaz = np.mean(viewazs)

    viewzens = np.array([np.float64(d["ZENITH_ANGLE"]["#text"]) for d in viewangles])
    viewzen = np.mean(viewzens)

    angles_dict = {
        "sun_azimuth": sunaz,
        "sun_zenith": sunzen,
        "view_azimuth": viewaz,
        "view_zenith": viewzen,
    }
    return angles_dict


def cog_get_s2_scl_data(aoi, item):
    SCL_NODATA = 99
    cog_crs = "EPSG:{}".format(item.properties["proj:epsg"])
    aoi_geometry_cog_crs = transform_crs(aoi.geometry, aoi.geometry_crs, cog_crs)
    bbox_cog_crs = list(aoi_geometry_cog_crs.bounds)

    # Scene Classification Band
    band = "SCL"
    # Transform aoi to pixel coordinates/window
    cog_transform = rasterio.transform.Affine(*item.assets[band]["proj:transform"][:-3])
    window = rasterio.windows.from_bounds(*bbox_cog_crs, cog_transform).round_offsets()

    # Get windowed data
    file_url = item.assets[band]["href"]
    # loop trough bands (file_url) here
    with rasterio.open(file_url) as src:
        kwds = src.profile
        raster_data = src.read(1, window=window, boundless=True)

        # Form a new clipped rasterio dataset
        transform = rasterio.windows.transform(window, kwds["transform"])
        height = raster_data.shape[-2]
        width = raster_data.shape[-1]

        new_kwds = kwds.copy()
        new_kwds.update(
            transform=transform,
            driver="GTiff",
            height=height,
            width=width,
            dtype=str(raster_data.dtype),
        )

    # Mask
    with MemoryFile() as memfile:
        with memfile.open(**new_kwds) as dataset:
            dataset.write(raster_data, indexes=1)

        masked_data, masked_kwds = mask_raster(
            memfile, aoi_geometry_cog_crs, no_data=SCL_NODATA
        )

    scl_dict = {"data": masked_data, "profile": masked_kwds}
    return scl_dict


def cog_generate_qi_dict(aoi, item, scl_data):

    date = pd.to_datetime(
        datetime.datetime.strptime(
            item.properties["sentinel:product_id"].split("_")[2], "%Y%m%dT%H%M%S"
        )
    )
    projection = {
        "type": "Projection",
        "crs": "EPSG:{}".format(item.properties["proj:epsg"]),
        "transform": item.assets["SCL"]["proj:transform"][:-3],
    }

    qi_dict = {
        "Date": date,
        "name": aoi.name,
        "tileid": "{}{}{}".format(
            item.properties["sentinel:utm_zone"],
            item.properties["sentinel:latitude_band"],
            item.properties["sentinel:grid_square"],
        ),
        "assetid": item.id,
        "productid": item.properties["sentinel:product_id"],
        "projection": projection,
        "datasource": "aws_cog",
    }

    # number of pixels inside aoi, excludes out-of-aoi pixels
    num_of_aoi_pixels = np.sum(scl_data != 99)
    for i, scl_class in enumerate(S2_SCL_CLASSES):

        class_percentage = np.sum(scl_data == i) / num_of_aoi_pixels

        qi_dict.update({scl_class: class_percentage})

    return qi_dict


def cog_get_s2_quality_info(aoi, req_params, items):

    qi_df = pd.DataFrame()
    for item in items:
        # print("Retrieving QI for item {}...".format(item.id))
        scl_dict = cog_get_s2_scl_data(aoi, item)
        qi_dict = cog_generate_qi_dict(aoi, item, scl_dict["data"])
        qi_df = qi_df.append(qi_dict, ignore_index=True)

    qi_df = qi_df.sort_values("Date").reset_index(drop=True)
    # Move "Date" to first columns
    cols = list(qi_df)
    cols.insert(0, cols.pop(cols.index("Date")))
    qi_df = qi_df.loc[:, cols]
    return qi_df


def cog_create_data_dict(aoi, item):
    date = pd.to_datetime(
        datetime.datetime.strptime(
            item.properties["sentinel:product_id"].split("_")[2], "%Y%m%dT%H%M%S"
        )
    )
    data_dict = {
        "Date": date,
        "name": aoi.name,
        "tileid": "{}{}{}".format(
            item.properties["sentinel:utm_zone"],
            item.properties["sentinel:latitude_band"],
            item.properties["sentinel:grid_square"],
        ),
        "assetid": item.id,
        "productid": item.properties["sentinel:product_id"],
        "projection": "EPSG:{}".format(item.properties["proj:epsg"]),
    }
    return data_dict


def cog_get_s2_band_data(
    aoi,
    req_params,
    items,
    qi_dataframe,
    qi_threshold=0.02,
    qi_filter=S2_FILTER1,
):

    filtered_qi = filter_s2_qi_dataframe(qi_dataframe, qi_threshold, qi_filter)
    if len(filtered_qi) == 0:
        print("No data to be retrieved for area %s !" % aoi.name)
        return None

    if aoi.tile is None:
        min_tile = min(filtered_qi["tileid"].values)
        filtered_qi = filtered_qi[filtered_qi["tileid"] == min_tile]
        aoi.tile = min_tile
    else:
        filtered_qi = filtered_qi[filtered_qi["tileid"] == aoi.tile]

    if len(filtered_qi) == 0:
        print("No qualified observations with used tile")
        return None
    print("{} good observations available. Retrieving data...".format(len(filtered_qi)))

    # Dataframe to collect all data
    data_df = pd.DataFrame()

    for item in items:
        if item.id not in filtered_qi["assetid"].values.tolist():
            continue
        print("Retrieving band data for item {}".format(item.id))

        # Data dict to hold single date data
        data_dict = cog_create_data_dict(aoi, item)
        # store all profiles to check that all profiles are equal as they should be
        data_dict["profile"] = []

        # Tranform aoi_geometry to raster crs
        cog_crs = "EPSG:{}".format(item.properties["proj:epsg"])
        aoi_geometry_cog_crs = transform_crs(aoi.geometry, aoi.geometry_crs, cog_crs)
        bbox_cog_crs = list(aoi_geometry_cog_crs.bounds)
        # Get larger window for resampling
        bbox_cog_crs = expand_bounds(bbox_cog_crs, 80)

        # currently always includes "SCL" data
        for band in req_params.bands + ["SCL"]:
            cog_transform = rasterio.transform.Affine(
                *item.assets[band]["proj:transform"][:-3]
            )
            window = rasterio.windows.from_bounds(
                *bbox_cog_crs, cog_transform
            ).round_offsets()

            file_url = item.assets[band]["href"]
            # loop trough bands (file_url) here
            with rasterio.open(file_url) as src:
                kwds = src.profile
                if band == "SCL":
                    raster_data = src.read(1, window=window, boundless=True)
                else:
                    raster_data = (
                        src.read(1, window=window, boundless=True) / S2_REFL_TRANS
                    )
                # Form a new clipped rasterio dataset
                transform = rasterio.windows.transform(window, kwds["transform"])
                height = raster_data.shape[-2]
                width = raster_data.shape[-1]

                new_kwds = kwds.copy()
                new_kwds.update(
                    transform=transform,
                    driver="GTiff",
                    height=height,
                    width=width,
                    dtype=str(raster_data.dtype),
                )
            # Resample data if needed
            with MemoryFile() as memfile:
                with memfile.open(**new_kwds) as dataset:  # Open as DatasetWriter
                    dataset.write(raster_data, 1)

                if req_params.target_gsd != new_kwds["transform"].a:
                    resampled_data, resampled_kwds = resample_raster(
                        memfile, req_params.target_gsd
                    )

                    data = resampled_data
                    new_kwds = resampled_kwds
                else:
                    data = raster_data[np.newaxis, ...]

            # Mask data
            with MemoryFile() as memfile:
                with memfile.open(**new_kwds) as dataset:
                    dataset.write(data)

                if band == "SCL":
                    no_data = 99
                else:
                    no_data = np.nan
                masked_data, masked_kwds = mask_raster(
                    memfile, aoi_geometry_cog_crs, no_data=no_data
                )

            data_dict.update({band: masked_data})
            data_dict["profile"].append(masked_kwds)

        angles_dict = get_angles(item)
        data_dict.update(angles_dict)
        data_df = data_df.append(data_dict, ignore_index=True)
    data_df = data_df.sort_values("Date").reset_index(drop=True)

    # Transform to xarray dataset
    data_ds = cog_s2_data_to_xarray(aoi, req_params, data_df)
    return data_ds


def check_shapes(dataframe, bands):

    # band to band comparison
    dataframe_cp = dataframe.copy()
    drop_these = []
    for image in dataframe_cp.itertuples(index=True):
        for band in bands[1:]:
            if getattr(image, bands[0]).shape != getattr(image, band).shape:

                print(
                    """Shape of bands doesn't match for image {} """.format(
                        image.productid
                    )
                )
                drop_these.append(image.Index)
                break
    # Drop images that have inconsistend shapes of bands
    dataframe.drop(index=drop_these, inplace=True)

    # solve most common shape
    shapes = list(dataframe[bands[0]].apply(lambda x: np.shape(x)))
    counts = Counter(shapes)
    most_common_shape = counts.most_common(1)[0][0]

    dataframe_cp = dataframe.copy()
    drop_these = []
    # image to image comparison
    for image in dataframe_cp.itertuples(index=True):

        if getattr(image, bands[0]).shape != most_common_shape:
            print("""Image {} uncommon shape. Dropping it.""".format(image.productid))
            print("Most common shape = {}".format(most_common_shape))
            print("Shape counts {}".format(counts))
            drop_these.append(image.Index)

    cleaned_dataframe = dataframe.drop(index=drop_these).reset_index(drop=True)

    return cleaned_dataframe


def cog_s2_data_to_xarray(aoi, req_params, dataframe):
    """"""

    #  2D data
    bands = req_params.bands

    #  1D data
    list_vars = [
        "assetid",
        "productid",
        "sun_azimuth",
        "sun_zenith",
        "view_azimuth",
        "view_zenith",
    ]

    # check that shapes of images and bands match
    dataframe = check_shapes(dataframe, bands)

    # crs from projection
    crs = dataframe["projection"][0]
    tileid = dataframe["tileid"][0]

    profile = dataframe.loc[0]["profile"][0]

    x_coords, y_coords = create_coordinate_arrays(profile)

    # translate to pixel center coordinates for netcdf/xarray dataset
    dx = profile["transform"].a
    dy = profile["transform"].e
    x_coords_center = x_coords + dx / 2
    y_coords_center = y_coords + dy / 2

    array = dataframe[bands].values
    # this will stack the array to ndarray with
    # dimension order (time, band, x,y)
    narray = np.stack(
        [np.stack(array[:, b], axis=2) for b in range(len(bands))], axis=2
    ).transpose()

    aoi_pixels = np.size(narray[0, 0, :, :]) - np.sum(np.isnan(narray[0, 0, :, :]))

    scl_array = np.stack(dataframe["SCL"].values, axis=2).transpose().astype(np.int16)

    coords = {
        "time": dataframe["Date"].values.astype(np.datetime64),
        "band": [
            b.replace("B0", "B") for b in bands
        ],  # switch to be consistent with gee
        "y": y_coords_center,
        "x": x_coords_center,
    }

    dataset_dict = {
        "band_data": (["time", "band", "x", "y"], narray),
        "SCL": (["time", "x", "y"], scl_array),
    }
    var_dict = {var: (["time"], dataframe[var]) for var in list_vars}
    source_series = pd.Series([req_params.datasource] * len(coords["time"]))
    var_dict.update(datasource=(["time"], source_series))
    dataset_dict.update(var_dict)

    ds = xr.Dataset(
        dataset_dict,
        coords=coords,
        attrs={
            "name": aoi.name,
            "crs": crs,
            "tile_id": tileid,
            "aoi_geometry": aoi.geometry.wkt,
            "aoi_pixels": aoi_pixels,
        },
    )
    ds = ds.transpose("time", "band", "y", "x")
    return ds
