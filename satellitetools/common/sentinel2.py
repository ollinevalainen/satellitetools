#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 10:53:15 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional

import numpy as np
import pandas as pd
import xarray as xr

from satellitetools.common.classes import AOI, DataSource

S2_REFL_TRANS = 10000
SCL_NODATA = 99


# Use GEE band naming
class S2Band(str, Enum):
    """Sentinel-2 bands."""

    B1 = "B1"
    B2 = "B2"
    B3 = "B3"
    B4 = "B4"
    B5 = "B5"
    B6 = "B6"
    B7 = "B7"
    B8 = "B8"
    B8A = "B8A"
    B9 = "B9"
    B11 = "B11"
    B12 = "B12"
    AOT = "AOT"
    WVP = "WVP"
    SCL = "SCL"

    def to_aws(self) -> str:
        """Convert band name to AWS band name.
        Returns:
        ----------------
        str
            Band name in AWS.
        """
        band_name_in_aws = S2_BANDS_GEE_TO_AWS[self.value]
        return band_name_in_aws

    def to_gee(self) -> str:
        """Convert band name to GEE band name.

        Returns:
        ----------------
        str
            Band name in GEE.
        """
        return self.value

    @classmethod
    def get_10m_to_20m_bands(cls) -> List["S2Band"]:
        """Get 10-20 m bands for the band.

        Returns:
        ----------------
        List[S2Band]
            List of 10-20 m bands.
        """
        return [S2Band(b) for b in S2_BANDS_10_20_GEE]

    @classmethod
    def get_all_bands(cls) -> List["S2Band"]:
        """Get all bands for the band.

        Returns:
        ----------------
        List[S2Band]
            List of all bands.
        """
        return [b for b in S2Band]


class SCLClass(Enum):
    """Sentinel-2 Scene Classification Layer (SCL) classes."""

    NODATA = 0
    SATURATED_DEFECTIVE = 1
    DARK_FEATURE_SHADOW = 2
    CLOUD_SHADOW = 3
    VEGETATION = 4
    NOT_VEGETATED = 5
    WATER = 6
    UNCLASSIFIED = 7
    CLOUD_MEDIUM_PROBA = 8
    CLOUD_HIGH_PROBA = 9
    THIN_CIRRUS = 10
    SNOW_ICE = 11


S2_BANDS_GEE = [
    "B1",
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "B8",
    "B8A",
    "B9",
    "B11",
    "B12",
    "AOT",
    "WVP",
    "SCL",
]
S2_BANDS_10_20_GEE = [
    "B2",
    "B3",
    "B4",
    "B5",
    "B6",
    "B7",
    "B8A",
    "B11",
    "B12",
]
# Old v0 earth search


# S2_BANDS_COG = [
#     "B01",
#     "B02",
#     "B03",
#     "B04",
#     "B05",
#     "B06",
#     "B07",
#     "B08",
#     "B8A",
#     "B09",
#     "B11",
#     "B12",
# ]

# # 10-20 m bands
# S2_BANDS_10_20_COG = ["B02", "B03", "B04", "B05", "B06", "B07", "B8A", "B11", "B12"]

S2_BANDS_COG = [
    "coastal",  # B1
    "blue",  # B2
    "green",  # B3
    "red",  # B4
    "rededge1",  # B5
    "rededge2",  # B6
    "rededge3",  # B7
    "nir",  # B8
    "nir08",  # B8A
    "nir09",  # B9
    "swir16",  # B11
    "swir22",  # B12
    "aot",
    "wvp",
    "scl",
]
S2_BANDS_10_20_COG = [
    "blue",
    "green",
    "red",
    "rededge1",
    "rededge2",
    "rededge3",
    "nir08",
    "swir16",
    "swir22",
]

S2_BANDS_AWS_TO_GEE = {
    aws: gee for aws, gee in zip(S2_BANDS_COG, S2_BANDS_GEE, strict=True)
}
# S2_BANDS_10_20_AWS_TO_GEE = {
#     aws: gee for aws, gee in zip(S2_BANDS_10_20_COG, S2_BANDS_10_20_GEE)
# }

S2_BANDS_GEE_TO_AWS = {
    gee: aws for gee, aws in zip(S2_BANDS_GEE, S2_BANDS_COG, strict=True)
}
# S2_BANDS_10_20_GEE_TO_AWS = {
#     gee: aws for gee, aws in zip(S2_BANDS_10_20_GEE, S2_BANDS_10_20_COG)
# }

S2_SCL_CLASSES = [c.name for c in SCLClass]

S2_FILTER1 = [
    SCLClass.NODATA,
    SCLClass.SATURATED_DEFECTIVE,
    SCLClass.CLOUD_SHADOW,
    SCLClass.UNCLASSIFIED,
    SCLClass.CLOUD_MEDIUM_PROBA,
    SCLClass.CLOUD_HIGH_PROBA,
    SCLClass.THIN_CIRRUS,
    SCLClass.SNOW_ICE,
]

S2_FILTER2 = [
    SCLClass.NODATA,
    SCLClass.SATURATED_DEFECTIVE,
    SCLClass.CLOUD_SHADOW,
    SCLClass.CLOUD_MEDIUM_PROBA,
    SCLClass.CLOUD_HIGH_PROBA,
    SCLClass.THIN_CIRRUS,
    SCLClass.SNOW_ICE,
]


class Sentinel2RequestParams:
    """S2 data request paramaters.

    Attributes:
    ----------------
    datestart : str
        Start date of the data request.
    dateend : str
        End date of the data request.
    datasource : DataSource
        Data source for the request.
    bands : List[S2Band]
        List of Sentinel-2 bands to request.
    target_gsd : float
        Target ground sampling distance (GSD) in meters.
    qi_evaluation_scale : float
        Quality indicator evaluation scale.

    """

    def __init__(
        self,
        datestart: str,
        dateend: str,
        datasource: DataSource,
        bands: List[S2Band] = None,
        target_gsd: float = 20,
        qi_evaluation_scale: float = 20,
    ):
        """Initialize the Sentinel2RequestParams class.

        Parameters:
        ----------------
        datestart : str
            Start date of the data request.
        dateend : str
            End date of the data request.
        datasource : DataSource
            Data source for the request.
        bands : List[S2Band]
            List of Sentinel-2 bands to request.
        target_gsd : float
            Target ground sampling distance (GSD) in meters.
        qi_evaluation_scale : float
            Quality indicator evaluation scale.
        """

        self.datestart = datestart
        self.dateend = dateend
        self.datasource = datasource
        if bands:
            self.bands = bands
        else:
            self.bands = [band for band in S2Band]
        self.target_gsd = target_gsd
        self.qi_evaluation_scale = qi_evaluation_scale

    def __repr__(self) -> str:
        return (
            f"RequestParams(datestart={self.datestart}, dateend={self.dateend}, "
            f"datasource={self.datasource}, bands={self.bands}, "
            f"target_gsd={self.target_gsd}, qi_evaluation_scale={self.qi_evaluation_scale})"
        )


class Sentinel2Metadata:
    """Sentinel-2 metadata class.

    Attributes:
    ----------------
    acquisition_time : pd.Timestamp
        Acquisition time of the image.
    tileid : str
        Sentinel-2 tile ID.
    assetid : str
        Asset ID of the image.
    productid : str
        Product ID of the image.
    projection : str
        Projection of the image.
    datasource : DataSource
        Data source of the image.
    observation_geometry : Optional[Sentinel2ObservationGeometry]
        Observation geometry of the image.
    class_percentages : Optional[Dict[SCLClass, float]]
        Class percentages of the image.

    """

    def __init__(
        self,
        acquisition_time: pd.Timestamp,
        tileid: str,
        assetid: str,
        productid: str,
        projection: str,
        datasource: DataSource,
        observation_geometry: Optional["Sentinel2ObservationGeometry"] = None,
        class_percentages: Optional[Dict[SCLClass, float]] = None,
    ):
        """Initialize the Sentinel2Metadata class.

        Parameters:
        ----------------
        acquisition_time : pd.Timestamp
            Acquisition time of the image.
        tileid : str
            Sentinel-2 tile ID.
        assetid : str
            Asset ID of the image.
        productid : str
            Product ID of the image.
        projection : str
            Projection of the image.
        datasource : DataSource
            Data source of the image.

        """

        self.acquisition_time = acquisition_time
        self.tileid = tileid
        self.assetid = assetid
        self.productid = productid
        self.projection = projection
        self.datasource = datasource

        # Placeholder for observation geometry
        self.observation_geometry = (
            observation_geometry if observation_geometry else None
        )
        # Placeholder for class percentages
        self.class_percentages = class_percentages if class_percentages else None

    def __repr__(self) -> str:
        return (
            f"Sentinel2Metadata(assetid={self.assetid}, "
            f"acquisition_time={self.acquisition_time}),"
            f"datasource={self.datasource}"
        )


@dataclass
class Sentinel2ObservationGeometry:
    """Sentinel-2 observation geometry class.

    Attributes:
    ----------------
    sun_azimuth : float
        Sun azimuth angle.
    sun_zenith : float
        Sun zenith angle.
    view_azimuth : float
        View azimuth angle.
    view_zenith : float
        View zenith angle.

    """

    sun_azimuth: float
    sun_zenith: float
    view_azimuth: float
    view_zenith: float


@dataclass
class Coordinates:
    """Coordinates class.

    Attributes:
    ----------------
    x : list
        List of x coordinates.
    y : list
        List of y coordinates.

    """

    x: list
    y: list


class Sentinel2Item:
    """Sentinel-2 data item class.

    Attributes:
    ----------------
    metadata : Sentinel2Metadata
        Metadata of the data item.
    data : Optional[Dict[S2Band, np.ndarray]]
        Data of the data item.
    coordinates : Dict[S2Band, Coordinates]
        Coordinates of the data item.

    """

    def __init__(
        self,
        metadata: Sentinel2Metadata,
        data: Optional[Dict[S2Band, np.ndarray]] = None,
    ):
        """Initialize the Sentinel2Item class.

        Parameters:
        ----------------
        metadata : Sentinel2Metadata
            Metadata of the data item.
        data : Optional[Dict[S2Band, np.ndarray]]
            Data of the data item.
        """
        self.metadata = metadata
        self.data = data if data else {}
        self.coordinates: Dict[S2Band, Coordinates] = {}

    def __repr__(self):
        return (
            f"Sentinel2Item(productid={self.metadata.productid}, "
            f"time={self.metadata.acquisition_time},"
            f"data={self.data.keys()})"
        )

    def set_coordinates(self, band: S2Band, xs: list, ys: list):
        """Set coordinates for the data item.

        Parameters:
        ----------------
        band : S2Band
            Sentinel-2 band.
        xs : list
            List of x coordinates.
        ys : list
            List of y coordinates.

        """
        self.coordinates[band] = Coordinates(x=xs, y=ys)


class Sentinel2DataCollection:
    """Sentinel-2 data collection class.

    Attributes:
    ----------------
    aoi : AOI
        Area of interest.
    req_params : Sentinel2RequestParams
        Request parameters.
    quality_information : pd.DataFrame
        Quality information.
    xr_dataset : xr.Dataset
        Xarray dataset.
    s2_items : List[Sentinel2Item]
        List of Sentinel-2 items.

    """

    def __init__(
        self,
        aoi: AOI,
        req_params: Sentinel2RequestParams,
    ):
        """Initialize the Sentinel2DataCollection class.

        Parameters:
        ----------------
        aoi : AOI
            Area of interest.
        req_params : Sentinel2RequestParams
            Request parameters.

        """
        self.aoi = aoi
        self.req_params = req_params
        self.quality_information: pd.DataFrame = None
        self.xr_dataset: xr.Dataset = None

        self.s2_items: List[Sentinel2Item] = None

    # Print information about the data collection
    def __repr__(self):
        return (
            f"Sentinel2DataCollection(AOI.name={self.aoi.name}, "
            f"Sentinel2RequestParams={self.req_params})"
        )

    def create_quality_information(self):
        """Create quality information dataframe."""
        qi_dicts = []
        for s2_item in self.s2_items:
            qi_dict = s2_item.metadata.__dict__
            qi_dict.update(s2_item.metadata.class_percentages)
            # Remove class percentages key as class percentages are added to the dictionary
            qi_dict.pop("class_percentages")
            # Remove observation geometry
            qi_dict.pop("observation_geometry")

            qi_dicts.append(qi_dict)

        # Make qi dataframe
        df_qi = pd.DataFrame(qi_dicts)
        df_qi.sort_values("acquisition_time", inplace=True)
        # Set acquisition time as index
        df_qi.set_index("acquisition_time", inplace=True)
        self.quality_information = df_qi

    def filter_s2_items(
        self, qi_threshold: float = 0, qi_filter: List[SCLClass] = S2_FILTER1
    ):
        """Filter S2 items based on quality information.

        Parameters:
        ----------------
        qi_threshold : float
            Quality indicator threshold.
        qi_filter : List[SCLClass]
            Quality indicator filter.

        """
        filtered_qi = filter_s2_qi_dataframe(
            self.quality_information, qi_threshold, qi_filter
        )

        # Filter to specified tile or use the first tile
        if self.aoi.tile is None:
            min_tile = min(filtered_qi["tileid"].values)
            filtered_qi = filtered_qi[filtered_qi["tileid"] == min_tile]
            self.aoi.tile = min_tile
        else:
            filtered_qi = filtered_qi[filtered_qi["tileid"] == self.aoi.tile]

        # IDs for images passing the quality filter
        assetids = filtered_qi["assetid"].values.tolist()

        # Remove items that do not pass the quality filter
        self.s2_items = [
            s2_item for s2_item in self.s2_items if s2_item.metadata.assetid in assetids
        ]

    def sort_s2_items(self):
        """Sort S2 items by acquisition time."""
        self.s2_items = sorted(self.s2_items, key=lambda x: x.metadata.acquisition_time)

    def data_to_xarray(self):
        """Convert data to xarray dataset."""

        # Sort s2_items by acquisition time
        self.sort_s2_items()

        #  2D data
        bands = self.req_params.bands
        bands_str = [b.value for b in bands]
        spectral_bands = [b for b in bands if b != S2Band.SCL]

        dataset_dict = {}
        multiband_arrays = []
        scl_arrays = []
        acquisition_times = []
        aoi_pixels = None

        all_metadata = {
            "assetid": [],
            "productid": [],
            "datasource": [],
            "sun_azimuth": [],
            "sun_zenith": [],
            "view_azimuth": [],
            "view_zenith": [],
        }

        for s2_item in self.s2_items:

            # Handle spectral data
            if len(spectral_bands) > 0:
                multiband_array = np.stack([s2_item.data[band] for band in bands])
                multiband_arrays.append(multiband_array)

            if S2Band.SCL in bands:
                scl_array = s2_item.data[S2Band.SCL]
                scl_arrays.append(scl_array)

            acquisition_times.append(np.datetime64(s2_item.metadata.acquisition_time))

            all_metadata["assetid"].append(s2_item.metadata.assetid)
            all_metadata["productid"].append(s2_item.metadata.productid)
            all_metadata["datasource"].append(s2_item.metadata.datasource.value)

            all_metadata["sun_azimuth"].append(
                s2_item.metadata.observation_geometry.sun_azimuth
            )
            all_metadata["sun_zenith"].append(
                s2_item.metadata.observation_geometry.sun_zenith
            )
            all_metadata["view_azimuth"].append(
                s2_item.metadata.observation_geometry.view_azimuth
            )
            all_metadata["view_zenith"].append(
                s2_item.metadata.observation_geometry.view_zenith
            )

        if multiband_arrays:
            multitemporal_array = np.stack(multiband_arrays)
            dataset_dict["band_data"] = (
                ["time", "band", "y", "x"],
                multitemporal_array,
            )

            aoi_pixels = np.size(multitemporal_array[0, 0, :, :]) - np.sum(
                np.isnan(multitemporal_array[0, 0, :, :])
            )

        if scl_arrays:
            multitemporal_scl_array = np.stack(scl_arrays).astype(np.int16)
            dataset_dict[S2Band.SCL.value] = (
                ["time", "y", "x"],
                multitemporal_scl_array,
            )
            # count aoi pixels if not already counted with spectral data
            if not aoi_pixels:
                aoi_pixels = np.size(multitemporal_scl_array[0, :, :]) - np.sum(
                    multitemporal_scl_array[0, :, :] == SCL_NODATA
                )

        # Add metadata to dataset
        dataset_dict.update(
            {
                var: (["time"], metadata_var)
                for var, metadata_var in all_metadata.items()
            }
        )

        # crs from projection
        crs = self.s2_items[0].metadata.projection
        tileid = self.s2_items[0].metadata.tileid

        coords = {
            "time": acquisition_times,
            "band": bands_str,
            "y": np.float64(self.s2_items[0].coordinates[bands[0]].y),
            "x": np.float64(self.s2_items[0].coordinates[bands[0]].x),
        }

        ds = xr.Dataset(
            dataset_dict,
            coords=coords,
            attrs={
                "name": self.aoi.name,
                "crs": crs,
                "tile_id": tileid,
                "aoi_geometry": self.aoi.geometry.wkt,
                "aoi_pixels": aoi_pixels,
            },
        )
        # ds = ds.transpose("time", "band", "y", "x") # TODO: Check if needed
        self.xr_dataset = ds


def filter_s2_qi_dataframe(
    s2_qi_dataframe: pd.DataFrame,
    qi_thresh: float,
    s2_filter: List[SCLClass] = S2_FILTER1,
) -> pd.DataFrame:
    """Filter qi dataframe.

    Parameters:
    ----------------
    s2_qi_dataframe : pd.DataFrame
        Quality information dataframe.
    qi_thresh : float
        Quality indicator threshold.
    s2_filter : List[SCLClass]
        Quality indicator filter.

    Returns:
    ----------------
    pd.DataFrame
        Filtered quality information dataframe.

    """

    # Drop if SCL data is nan
    s2_qi_dataframe = s2_qi_dataframe.dropna(axis=0, subset=S2_SCL_CLASSES)
    s2_filter_str = [scl.name for scl in s2_filter]
    filtered_s2_qi_df = s2_qi_dataframe.loc[
        s2_qi_dataframe[s2_filter_str].sum(axis=1) <= qi_thresh
    ]

    return filtered_s2_qi_df
