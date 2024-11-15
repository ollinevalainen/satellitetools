#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module to retrieve Sentinel-2 data from Google Earth Engine (GEE).

@author: Olli Nevalainen (olli.nevalainen@fmi.fi),
 Finnish Meteorological Institute)

Created on Thu Feb  6 15:24:12 2020

"""
import datetime
from typing import Dict, List, Tuple, Union
from warnings import warn

import ee
import numpy as np
import pandas as pd

from satellitetools.common.classes import AOI, DataSource
from satellitetools.common.sentinel2 import (
    S2_FILTER1,
    S2_REFL_TRANS,
    SCL_NODATA,
    Coordinates,
    S2Band,
    SCLClass,
    Sentinel2DataCollection,
    Sentinel2Item,
    Sentinel2Metadata,
    Sentinel2ObservationGeometry,
    Sentinel2RequestParams,
)

NO_DATA = -99999
GEE_DATASET = "COPERNICUS/S2_SR_HARMONIZED"
GEE_SCL_BAND = S2Band.SCL.to_gee()


class GEESentinel2DataCollection(Sentinel2DataCollection):
    """Class to retrieve Sentinel-2 data from Google Earth Engine (GEE).

    Attributes
    ----------
    aoi : AOI
        Area of interest.
    req_params : Sentinel2RequestParams
        Request parameters.

    """

    def __init__(self, aoi: AOI, req_params: Sentinel2RequestParams):
        """Initialize GEESentinel2DataCollection object.

        Parameters
        ----------
        aoi : AOI
            Area of interest.
        req_params : Sentinel2RequestParams
            Request parameters.
        """
        super().__init__(aoi, req_params)

    def create_s2_items(self, feature: Dict):
        """Create Sentinel-2 items from GEE feature.

        Parameters
        ----------
        feature : dict
            Feature dictionary from GEE.

        """

        s2_items = []
        properties = feature["properties"]
        productids = properties["productid"]
        for i, productid in enumerate(productids):

            acquisition_time = datetime.datetime.strptime(
                productid.split("_")[2], "%Y%m%dT%H%M%S"
            )

            metadata_dict = {
                "acquisition_time": acquisition_time,
                "tileid": properties["tileid"][i],
                "assetid": properties["assetid"][i],
                "productid": productid,
                "projection": properties["projection"][i]["crs"],
                "datasource": DataSource.GEE,
            }
            # If observation geometry is available, set it
            try:
                observation_geometry = Sentinel2ObservationGeometry(
                    sun_azimuth=properties["sun_azimuth"][i],
                    sun_zenith=properties["sun_zenith"][i],
                    view_azimuth=properties["view_azimuth"][i],
                    view_zenith=properties["view_zenith"][i],
                )
                metadata_dict["observation_geometry"] = observation_geometry
            except KeyError:
                # There is no observation geometry information available, so let's not set it
                pass

            class_percentages = {
                scl_class.name: properties[scl_class.name][i]
                for scl_class in SCLClass
                if scl_class.name in properties
            }
            if class_percentages:
                metadata_dict["class_percentages"] = class_percentages
            metadata = Sentinel2Metadata(**metadata_dict)
            s2_item = Sentinel2Item(metadata)
            s2_items.append(s2_item)

        self.s2_items = s2_items
        self.sort_s2_items()

    def create_s2_items_from_quality_information(self):
        """Create Sentinel-2 items from quality information."""

        feature_dict = {}
        feature_dict["properties"] = self.quality_information.to_dict("list")
        self.create_s2_items(feature_dict)
        self.sort_s2_items()

    def set_s2_item_data(self, feature: Dict):
        """Set Sentinel-2 item data from GEE feature.

        Parameters
        ----------
        feature : dict
            Feature dictionary from GEE.

        """
        self.sort_s2_items()
        properties = feature["properties"]
        bands = [b.value for b in self.req_params.bands]
        for i, productid in enumerate(properties["productid"]):
            coordinates = {}
            band_data = {}

            try:
                observation_geometry = Sentinel2ObservationGeometry(
                    sun_azimuth=properties["sun_azimuth"][i],
                    sun_zenith=properties["sun_zenith"][i],
                    view_azimuth=properties["view_azimuth"][i],
                    view_zenith=properties["view_zenith"][i],
                )
            except KeyError:  # No observation geometry available
                observation_geometry = None

            # Check data consistency
            band_lens = [len(properties[band][i]) for band in bands]
            coord_lens = [
                len(properties["x_coords"][i]),
                len(properties["y_coords"][i]),
            ]
            unique_lens = set(band_lens + coord_lens)
            if len(unique_lens) > 1:
                print(
                    f"Data lengths are not consistent for {productid}. "
                    "Dropping product."
                )
                self.drop_s2_item_with_productid(productid)
                continue
            for band in bands:
                ys = properties["y_coords"][i]
                xs = properties["x_coords"][i]
                band_arr = properties[band][i]
                if band == S2Band.SCL.value:
                    xs, ys, band_arr = gee_lists_to_2D_array(
                        xs, ys, band_arr, no_data=SCL_NODATA
                    )
                    band_arr = band_arr.astype(np.int32)
                else:
                    xs, ys, band_arr = gee_lists_to_2D_array(
                        xs, ys, band_arr, no_data=np.nan
                    )
                    band_arr = band_arr.astype(np.float64) / S2_REFL_TRANS

                coordinates.update({band: Coordinates(xs, ys)})
                band_data.update({band: band_arr})
            s2_item = self.get_s2_item_with_productid(productid)

            if s2_item:
                s2_item.metadata.observation_geometry = observation_geometry
                s2_item.data = band_data
                s2_item.coordinates = coordinates

    def get_s2_item_with_productid(self, productid: str):
        """Get Sentinel-2 item with productid.

        Parameters
        ----------
        productid : str
            Sentinel-2 product ID.

        Returns
        -------
        s2_item : Sentinel2Item
            Sentinel-2 item.

        """
        for s2_item in self.s2_items:
            if s2_item.metadata.productid == productid:
                return s2_item
        return None

    def drop_s2_item_with_productid(self, productid: str):
        """Drop Sentinel-2 item with productid.

        Parameters
        ----------
        productid : str
            Sentinel-2 product ID.

        """
        self.s2_items = [
            s2_item
            for s2_item in self.s2_items
            if s2_item.metadata.productid != productid
        ]

    def ee_feature_from_s2_items(self) -> ee.Feature:
        """Create GEE feature from Sentinel-2 items.

        Returns
        -------
        feature : ee.Feature
            GEE feature.

        """

        if self.s2_items is None:
            print(f"No S2 items available for {self.aoi.name}. Get quality info first.")
            return None

        gee_assetids = []
        for s2_item in self.s2_items:
            gee_assetids.append(GEE_DATASET + "/" + s2_item.metadata.assetid)

        image_list = [ee.Image(asset_id) for asset_id in gee_assetids]
        crs = self.s2_items[0].metadata.projection
        feature = ee.Feature(
            ee.Geometry.Polygon(list(self.aoi.geometry.exterior.coords)),
            {"name": self.aoi.name, "image_list": image_list, "crs": crs},
        )

        return feature

    def ee_feature_from_aoi(self) -> ee.Feature:
        """Create GEE feature from AOI.

        Returns
        -------
        feature : ee.Feature
            GEE feature.

        """

        feature = ee.Feature(
            ee.Geometry.Polygon(list(self.aoi.geometry.exterior.coords)),
            {"name": self.aoi.name},
        )
        return feature

    def get_quality_info(self):
        """Get quality information from GEE."""

        multi_data_collection = MultiGEESentinel2DataCollection(
            [self.aoi], self.req_params
        )
        print(
            "Searching and computing quality information for S2 data "
            "from {} to {} for {}".format(
                self.req_params.datestart, self.req_params.dateend, self.aoi.name
            )
        )
        multi_data_collection.ee_get_s2_quality_info()

        self.s2_items = multi_data_collection.data_collections[self.aoi.name].s2_items
        self.quality_information = multi_data_collection.data_collections[
            self.aoi.name
        ].quality_information
        print(f"Found {len(self.s2_items)} items.")

    def get_s2_data(self):
        """Get S2 data (level L2A, bottom of atmosphere data) from GEE."""
        multi_data_collection = MultiGEESentinel2DataCollection(
            [self.aoi], self.req_params
        )
        multi_data_collection.data_collections[self.aoi.name].quality_information = (
            self.quality_information
        )
        multi_data_collection.data_collections[self.aoi.name].s2_items = self.s2_items

        if self.s2_items:
            print("Retrieving data...")
        multi_data_collection.ee_get_s2_data()
        self.s2_items = multi_data_collection.data_collections[self.aoi.name].s2_items

    def search_s2_items(self):
        multi_data_collection = MultiGEESentinel2DataCollection(
            [self.aoi], self.req_params
        )
        print(
            "Searching S2 data from {} to {} for area {}".format(
                self.req_params.datestart, self.req_params.dateend, self.aoi.name
            )
        )
        multi_data_collection.ee_search_s2_products()
        self.s2_items = multi_data_collection.data_collections[self.aoi.name].s2_items
        print(f"Found {len(self.s2_items)} items.")


class MultiGEESentinel2DataCollection:
    """Class to retrieve Sentinel-2 data from Google Earth Engine (GEE) for
    multiple AOIs.

    Attributes
    ----------
    aois : List[AOI]
        List of AOIs.
    req_params : Sentinel2RequestParams
        Request parameters.

    """

    def __init__(self, aois: List[AOI], req_params: Sentinel2RequestParams):
        """Initialize MultiGEESentinel2DataCollection object.

        Parameters
        ----------
        aois : List[AOI]
            List of AOIs.
        req_params : Sentinel2RequestParams
            Request parameters.
        """
        self.data_collections = {
            aoi.name: GEESentinel2DataCollection(aoi, req_params) for aoi in aois
        }
        self.req_params = req_params

    def ee_search_s2_products(self):

        req_params = self.req_params
        features = [
            data_collection.ee_feature_from_aoi()
            for data_collection in self.data_collections.values()
        ]

        feature_collection = ee.FeatureCollection(features)

        def ee_get_s2_image_collection_metadata(feature):
            """Get S2 metadata for images.

            Parameters
            ----------
            feature : ee.Feature
                GEE feature.

            """

            area = feature.geometry()
            image_collection = (
                ee.ImageCollection(GEE_DATASET)
                .filterBounds(area)
                .filterDate(req_params.datestart, req_params.dateend)
                .select([GEE_SCL_BAND])
            )

            def ee_get_s2_image_metadata(img):
                productid = img.get("PRODUCT_ID")
                assetid = img.id()
                tileid = img.get("MGRS_TILE")
                system_index = img.get("system:index")
                proj = img.select(GEE_SCL_BAND).projection()

                tmpfeature = (
                    ee.Feature(ee.Geometry.Point([0, 0]))
                    .set("productid", productid)
                    .set("system_index", system_index)
                    .set("assetid", assetid)
                    .set("tileid", tileid)
                    .set("projection", proj)
                )
                return tmpfeature

            s2_qi_image_collection = image_collection.map(ee_get_s2_image_metadata)

            return (
                feature.set(
                    "productid", s2_qi_image_collection.aggregate_array("productid")
                )
                .set(
                    "system_index",
                    s2_qi_image_collection.aggregate_array("system_index"),
                )
                .set("assetid", s2_qi_image_collection.aggregate_array("assetid"))
                .set("tileid", s2_qi_image_collection.aggregate_array("tileid"))
                .set("projection", s2_qi_image_collection.aggregate_array("projection"))
            )

        s2_meta_feature_collection = feature_collection.map(
            ee_get_s2_image_collection_metadata
        ).getInfo()

        for feature in s2_meta_feature_collection["features"]:
            name = feature["properties"]["name"]
            self.data_collections[name].create_s2_items(feature)

    def ee_get_s2_quality_info(self):
        """Get S2 quality information from GEE."""

        req_params = self.req_params
        features = [
            data_collection.ee_feature_from_aoi()
            for data_collection in self.data_collections.values()
        ]

        feature_collection = ee.FeatureCollection(features)

        def ee_get_s2_quality_info_feature(feature):
            """Get S2 quality information from GEE feature.

            Parameters
            ----------
            feature : ee.Feature
                GEE feature.

            """

            area = feature.geometry()
            image_collection = (
                ee.ImageCollection(GEE_DATASET)
                .filterBounds(area)
                .filterDate(req_params.datestart, req_params.dateend)
                .select([GEE_SCL_BAND])
            )

            def ee_get_s2_quality_info_image(img):
                productid = img.get("PRODUCT_ID")
                assetid = img.id()
                tileid = img.get("MGRS_TILE")
                system_index = img.get("system:index")
                proj = img.select(GEE_SCL_BAND).projection()

                # apply reducer to list
                img = img.reduceRegion(
                    reducer=ee.Reducer.toList(),
                    geometry=area,
                    maxPixels=1e8,
                    scale=req_params.qi_evaluation_scale,
                )

                # get data into arrays
                classdata = ee.Array(
                    ee.Algorithms.If(
                        img.get(GEE_SCL_BAND),
                        ee.Array(img.get(GEE_SCL_BAND)),
                        ee.Array([0]),
                    )
                )

                totalcount = classdata.length()
                classpercentages = {
                    scl.name: classdata.eq(scl.value)
                    .reduce(ee.Reducer.sum(), [0])
                    .divide(totalcount)
                    .get([0])
                    for scl in SCLClass
                }

                tmpfeature = (
                    ee.Feature(ee.Geometry.Point([0, 0]))
                    .set("productid", productid)
                    .set("system_index", system_index)
                    .set("assetid", assetid)
                    .set("tileid", tileid)
                    .set("projection", proj)
                    .set(classpercentages)
                )
                return tmpfeature

            s2_qi_image_collection = image_collection.map(ee_get_s2_quality_info_image)

            return (
                feature.set(
                    "productid", s2_qi_image_collection.aggregate_array("productid")
                )
                .set(
                    "system_index",
                    s2_qi_image_collection.aggregate_array("system_index"),
                )
                .set("assetid", s2_qi_image_collection.aggregate_array("assetid"))
                .set("tileid", s2_qi_image_collection.aggregate_array("tileid"))
                .set("projection", s2_qi_image_collection.aggregate_array("projection"))
                .set(
                    {
                        scl_class.name: s2_qi_image_collection.aggregate_array(
                            scl_class.name
                        )
                        for scl_class in SCLClass
                    }
                )
            )

        qi_feature_collection = feature_collection.map(
            ee_get_s2_quality_info_feature
        ).getInfo()

        for feature in qi_feature_collection["features"]:
            name = feature["properties"]["name"]
            self.data_collections[name].create_s2_items(feature)
            self.data_collections[name].create_quality_information()

    def ee_get_s2_data(
        self,
    ):
        """Get S2 data (level L2A, bottom of atmosphere data) from GEE.

        Warning! All bands are resampled to resolution specified by
        req_params.target_gsd.

        """
        datestart = self.req_params.datestart
        dateend = self.req_params.dateend
        bands = [b.value for b in self.req_params.bands]  # Convert to strings for GEE
        spectral_bands = [b for b in bands if b != S2Band.SCL.value]
        resolution = self.req_params.target_gsd

        features = []
        for data_collection in self.data_collections.values():
            feature = data_collection.ee_feature_from_s2_items()
            features.append(feature)

        if len(features) == 0:
            print("No data to be retrieved.")
            return None

        feature_collection = ee.FeatureCollection(features)

        def ee_get_s2_feature_data(feature):
            crs = feature.get("crs")
            geom = feature.geometry(0.01, crs)
            image_collection = (
                ee.ImageCollection.fromImages(feature.get("image_list"))
                .filterBounds(geom)
                .filterDate(datestart, dateend)
                .select(bands)
            )

            def ee_get_s2_image_data(img):
                # img = img.clip(geom)
                productid = img.get("PRODUCT_ID")
                assetid = img.id()
                tileid = img.get("MGRS_TILE")
                system_index = img.get("system:index")
                proj = img.select(bands[0]).projection()
                sun_azimuth = img.get("MEAN_SOLAR_AZIMUTH_ANGLE")
                sun_zenith = img.get("MEAN_SOLAR_ZENITH_ANGLE")
                if spectral_bands:
                    view_azimuth = (
                        ee.Array(
                            [
                                img.get("MEAN_INCIDENCE_AZIMUTH_ANGLE_%s" % b)
                                for b in spectral_bands
                            ]
                        )
                        .reduce(ee.Reducer.mean(), [0])
                        .get([0])
                    )
                    view_zenith = (
                        ee.Array(
                            [
                                img.get("MEAN_INCIDENCE_ZENITH_ANGLE_%s" % b)
                                for b in spectral_bands
                            ]
                        )
                        .reduce(ee.Reducer.mean(), [0])
                        .get([0])
                    )

                # img = img.resample("bilinear").reproject(crs=crs, scale=resolution)

                # get the lat lon
                image_grid = ee.Image.pixelCoordinates(ee.Projection(crs)).addBands(
                    [img.select(b) for b in bands]
                )

                # apply reducer to list
                image_grid = image_grid.reduceRegion(
                    reducer=ee.Reducer.toList(),
                    geometry=geom,
                    maxPixels=1e8,
                    scale=resolution,
                )

                # get data into 1D arrays
                x_coords = ee.Array(image_grid.get("x"))
                y_coords = ee.Array(image_grid.get("y"))
                band_data = {b: ee.Array(image_grid.get("%s" % b)) for b in bands}

                # perform LAI et al. computation possibly here!

                tmpfeature = (
                    ee.Feature(ee.Geometry.Point([0, 0]))
                    .set("productid", productid)
                    .set("system_index", system_index)
                    .set("assetid", assetid)
                    .set("tileid", tileid)
                    .set("projection", proj)
                    .set("sun_zenith", sun_zenith)
                    .set("sun_azimuth", sun_azimuth)
                    .set("x_coords", x_coords)
                    .set("y_coords", y_coords)
                    .set(band_data)
                )
                if spectral_bands:
                    tmpfeature = tmpfeature.set("view_zenith", view_zenith)
                    tmpfeature = tmpfeature.set("view_azimuth", view_azimuth)

                return tmpfeature

            s2_data_feature = image_collection.map(ee_get_s2_image_data)

            feature = (
                feature.set("productid", s2_data_feature.aggregate_array("productid"))
                .set("system_index", s2_data_feature.aggregate_array("system_index"))
                .set("assetid", s2_data_feature.aggregate_array("assetid"))
                .set("tileid", s2_data_feature.aggregate_array("tileid"))
                .set("projection", s2_data_feature.aggregate_array("projection"))
                .set("sun_zenith", s2_data_feature.aggregate_array("sun_zenith"))
                .set("sun_azimuth", s2_data_feature.aggregate_array("sun_azimuth"))
                .set("x_coords", s2_data_feature.aggregate_array("x_coords"))
                .set("y_coords", s2_data_feature.aggregate_array("y_coords"))
                .set({b: s2_data_feature.aggregate_array(b) for b in bands})
            )
            if spectral_bands:
                feature = feature.set(
                    "view_zenith", s2_data_feature.aggregate_array("view_zenith")
                )
                feature = feature.set(
                    "view_azimuth", s2_data_feature.aggregate_array("view_azimuth")
                )
            return feature

        s2_data_feature_collection = feature_collection.map(
            ee_get_s2_feature_data
        ).getInfo()

        for feature in s2_data_feature_collection["features"]:
            name = feature["properties"]["name"]
            self.data_collections[name].set_s2_item_data(feature)


def gee_lists_to_2D_array(
    x_coords: Union[np.ndarray, list],
    y_coords: Union[np.ndarray, list],
    data: Union[np.ndarray, list],
    no_data: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert 1D lists of coordinates and corresponding values to 2D array.

    Parameters
    ----------
    x_coords : Union[np.ndarray, list]
        X coordinates.
    y_coords : Union[np.ndarray, list]
        Y coordinates.
    data : Union[np.ndarray, list]
        Data values.
    no_data : float
        No data value.

    Returns
    -------
    x_vals : np.ndarray
        X values.
    y_vals : np.ndarray
        Y values.
    arr : np.ndarray
        2D array of data values.

    """

    y_vals, y_idx = np.unique(y_coords, return_inverse=True)
    y_vals = np.flip(y_vals)
    x_vals, x_idx = np.unique(x_coords, return_inverse=True)
    arr = np.empty(y_vals.shape + x_vals.shape)
    arr.fill(no_data)
    arr[y_idx, x_idx] = np.array(data)
    arr = np.flipud(arr)
    return x_vals, y_vals, arr


def get_copernicus_dem_elevation(lat: float, lon: float):
    """Get elevation from Copernicus DEM.

    Parameters
    ----------
    lat : float
        Latitude.
    lon : float
        Longitude.

    Returns
    -------
    elevation : float
        Elevation.
    """

    point = ee.Geometry.Point(lon, lat)
    dataset = (
        ee.ImageCollection("COPERNICUS/DEM/GLO30").select("DEM").filterBounds(point)
    )

    # returns list of two lists where first list is headers and the second are data
    point_data: list = dataset.getRegion(point, 30).getInfo()
    for key, value in zip(point_data[0], point_data[1], strict=True):
        if key == "DEM":
            return value


# To be deprecated functions


def ee_get_s2_quality_info(
    aois: Union[str, List[AOI]], req_params: Sentinel2RequestParams
) -> Dict[str, pd.DataFrame]:
    """Get S2 quality information from GEE.

    Parameters
    ----------
    aois : Union[str, List[AOI]]
        Area of interest.
    req_params : Sentinel2RequestParams
        Request parameters.

    Returns
    -------
    quality_dict : Dict[str, pd.DataFrame]
        Quality information.

    """
    DEPRECATION_WARNING_TEXT = (
        "Function ee_get_s2_quality_info will be deprecated and removed at some point."
        "Use GEESentinel2DataCollection, MultiGEESentinel2DataCollection or function "
        "get_s2_qi_and_data() in satellitetools.common.wrappers instead."
    )
    warn(DEPRECATION_WARNING_TEXT, DeprecationWarning, stacklevel=2)

    if isinstance(aois, str):
        aois = [aois]

    multi_data_collection = MultiGEESentinel2DataCollection(aois, req_params)
    multi_data_collection.ee_get_s2_quality_info()

    quality_dict = {}
    for aoi_name, data_collection in multi_data_collection.data_collections.items():
        quality_dict[aoi_name] = data_collection.quality_information

    return quality_dict


def ee_get_s2_data(
    aois: Union[str, List[AOI]],
    req_params: Sentinel2RequestParams,
    quality_dict: Dict[str, pd.DataFrame],
    qi_threshold: float = 0,
    qi_filter: List[SCLClass] = S2_FILTER1,
) -> Dict[str, pd.DataFrame]:
    """Get S2 data (level L2A, bottom of atmosphere data) from GEE.

    Parameters
    ----------
    aois : Union[str, List[AOI]]
        Area of interest.
    req_params : Sentinel2RequestParams
        Request parameters.
    quality_dict : Dict[str, pd.DataFrame]
        Quality information.
    qi_threshold : float, optional
        Quality index threshold, by default 0.
    qi_filter : List[SCLClass], optional
        Quality index filter, by default S2_FILTER1.

    Returns
    -------
    data_dict : Dict[str, pd.DataFrame]
        Data.

    """
    DEPRECATION_WARNING_TEXT = (
        "Function ee_get_s2_data will be deprecated and removed at some point."
        "Use GEESentinel2DataCollection, MultiGEESentinel2DataCollection or function"
        " get_s2_qi_and_data() in satellitetools.common.wrappers instead."
    )
    warn(DEPRECATION_WARNING_TEXT, DeprecationWarning, stacklevel=2)

    if isinstance(aois, str):
        aois = [aois]

    multi_data_collection = MultiGEESentinel2DataCollection(aois, req_params)
    for data_collection in multi_data_collection.data_collections.values():
        data_collection.quality_information = quality_dict[data_collection.aoi.name]
        data_collection.create_s2_items_from_quality_information()
        data_collection.filter_s2_items(qi_threshold, qi_filter)

    multi_data_collection.ee_get_s2_data()

    data_dict = {}
    for aoi_name, data_collection in multi_data_collection.data_collections.items():
        data_collection.data_to_xarray()
        data_dict[aoi_name] = data_collection.xr_dataset
    return data_dict
