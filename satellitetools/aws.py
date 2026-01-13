#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module to get Sentinel-2 data from AWS Open data registry,
where Sentinel-2 (level 2A) data is available as cloud-optimized
geotiffs (https://registry.opendata.aws/sentinel-2-l2a-cogs/).


@author: Olli Nevalainen (Finnish Meteorological Institute)
"""

import datetime
import logging
import urllib
from multiprocessing import Pool
from typing import Dict, List, Optional, Tuple, Union
from warnings import warn

try:
    # breaking change introduced in python 3.11
    from enum import StrEnum
except ImportError:
    from enum import Enum

    class StrEnum(str, Enum):
        pass


import numpy as np
import pandas as pd
import rasterio
import xmltodict
from pystac.item import Item
from pystac_client import Client
from rasterio import MemoryFile

from satellitetools.common.classes import AOI, DataSource
from satellitetools.common.raster import mask_raster, reproject_data_to_profile
from satellitetools.common.sentinel2 import (
    SCL_NODATA,
    SPECTRAL_BAND_NO_DATA,
    S2Band,
    SCLClass,
    Sentinel2DataCollection,
    Sentinel2Item,
    Sentinel2Metadata,
    Sentinel2ObservationGeometry,
    Sentinel2RequestParams,
)
from satellitetools.common.vector import (
    coordinate_arrays_from_profile,
    expand_bounds,
    transform_crs,
)

logger = logging.getLogger(__name__)

BUFFER_MULTIPLIER = 8


class EarthSearchCollection(StrEnum):
    SENTINEL2_C1_L2A = "sentinel-2-c1-l2a"
    SENTINEL2_L2A = "sentinel-2-l2a"


class EarthSearch:
    """Class to handle search for items in EarthSearch.

    Attributes:
    -----------
    datestart: Union[str, pd.Timestamp, datetime.datetime]
        Start date for search
    dateend: Union[str, pd.Timestamp, datetime.datetime]
        End date for search
    bbox: List[float]
        Bounding box coordinates [minx, miny, maxx, maxy]
    collection: EarthSearchCollection
        Collection to search
    limit: int
        Limit for search results

    Note:
    -----
    EarthSearch API documentation at:
    https://earth-search.aws.element84.com/v1/api.html#tag/Item-Search/operation/getItemSearch

    """

    EARTH_SEARCH_ENDPOINT = "https://earth-search.aws.element84.com/v1"
    DEFAULT_REQUEST_LIMIT = 10000

    def __init__(
        self,
        datestart: Union[str, pd.Timestamp, datetime.datetime],
        dateend: Union[str, pd.Timestamp, datetime.datetime],
        bbox: List[float],
        collection: EarthSearchCollection,
    ):
        """Initialize EarthSearch object.

        Parameters:
        -----------
        datestart: Union[str, pd.Timestamp, datetime.datetime]
            Start date for search
        dateend: Union[str, pd.Timestamp, datetime.datetime]
            End date for search
        bbox: List[float]
            Bounding box coordinates [minx, miny, maxx, maxy]
        collection: EarthSearchCollection
            Collection to search

        """

        self.datestart = pd.to_datetime(datestart)
        self.dateend = pd.to_datetime(dateend)
        if self.datestart > self.dateend:
            logger.error("datestart must be before dateend.")
            raise ValueError("datestart must be before dateend.")
        self.bbox = bbox
        self.collection = collection
        self.limit = self.DEFAULT_REQUEST_LIMIT

    def search_collection(
        self,
        datestart: pd.Timestamp,
        dateend: pd.Timestamp,
        collection: EarthSearchCollection,
    ) -> List[Item]:
        """Search for items in EarthSearch collection.

        Parameters:
        -----------
        datestart: pd.Timestamp
            Start date for search
        dateend: pd.Timestamp
            End date for search
        collection: EarthSearchCollection
            Collection to search

        Returns:
        --------
        all_items: List[Item]
            List of items

        """

        # Split queries to half year time ranges
        time_ranges = split_time_range(datestart, dateend)
        all_items = []
        for time_range in time_ranges:
            dates = "{}/{}".format(
                time_range[0].isoformat() + "Z",
                time_range[1].isoformat() + "Z",
            )
            # Search
            client = Client.open(self.EARTH_SEARCH_ENDPOINT)
            search = client.search(
                collections=[collection],
                datetime=dates,
                bbox=self.bbox,
                limit=self.limit,
            )

            if search.matched() > 0:
                all_items.extend(search.item_collection())

        logger.info("Found {} items.".format(len(all_items)))
        return all_items

    def get_items(self) -> List[Item]:
        """Get items from EarthSearch.

        Returns:
        --------
        all_items: List[Item]
            List of items

        """

        all_items = []
        # Search for items
        items = self.search_collection(self.datestart, self.dateend, self.collection)
        if items:
            all_items.extend(items)

        # Search for items from 2022 in EarthSearchCollection.SENTINEL2_L2A since at the
        # moment (2024-10-23) SENTINEL2_C1_L2A collection is still missing that year
        if (
            self.collection == EarthSearchCollection.SENTINEL2_C1_L2A
            and self.datestart.year <= 2022
            and self.dateend.year >= 2022
        ):
            if self.datestart.year < 2022:
                datestart_2022 = pd.Timestamp("2022-01-01")
            else:
                datestart_2022 = self.datestart
            if self.dateend.year > 2022:
                dateend_2022 = pd.Timestamp("2022-12-31")
            else:
                dateend_2022 = self.dateend
            items_2022 = self.search_collection(
                datestart_2022, dateend_2022, EarthSearchCollection.SENTINEL2_L2A
            )
            if items_2022:
                all_items.extend(items_2022)
                # Check if there's duplicate items from SENTINEL2_C1_L2A and
                # SENTINEL2_L2A. Keep the ones from SENTINEL2_C1_L2A.
                all_items = remove_duplicate_items(all_items)
        return all_items


def remove_duplicate_items(items) -> List[Item]:
    """Remove duplicate items from list of items.

    Parameters:
    -----------
    items: List[Item]
        List of items

    Returns:
    --------
    filtered_items: List[Item]
        Filtered list of items

    """

    # Find duplicate items (same "s2:product_uri" )
    duplicate_product_ids = []
    all_product_ids = []
    for item in items:
        product_id = item.properties["s2:product_uri"]
        all_product_ids.append(product_id)
        if product_id in all_product_ids:
            duplicate_product_ids.append(product_id)

    # For duplicate items keep the one with
    # properties["processing:software"] == "sentinel-2-c1-l2a-to-stac"
    # and remove the one with properties["processing:software"] == "sentinel2-to-stac"
    filtered_items = []
    for item in items:
        product_id = item.properties["s2:product_uri"]
        if product_id in duplicate_product_ids:
            if "sentinel-2-c1-l2a-to-stac" in item.properties["processing:software"]:
                filtered_items.append(item)
        else:
            filtered_items.append(item)
    return items


class AWSSentinel2DataCollection(Sentinel2DataCollection):
    """Class to handle Sentinel-2 data from AWS Open data registry.

    Attributes:
    -----------
    aoi: AOI
        Area of interest
    req_params: Sentinel2RequestParams
        Request parameters
    s2_items: List[AWSSentinel2Item]
        List of Sentinel-2 items
    multiprocessing: Optional[int]
        Number of processes to use in multiprocessing

    """

    def __init__(
        self,
        aoi: AOI,
        req_params: Sentinel2RequestParams,
        multiprocessing: Optional[int] = None,
    ):
        """Initialize AWSSentinel2DataCollection object.

        Parameters:
        -----------
        aoi: AOI
            Area of interest
        req_params: Sentinel2RequestParams
            Request parameters
        multiprocessing: Optional[int]
            Number of processes to use in multiprocessing
        """
        super().__init__(aoi, req_params)
        self.s2_items: List[AWSSentinel2Item] = None
        self.multiprocessing = multiprocessing

    def search_s2_items(self):
        """Search for Sentinel-2 items from AWS Open data registry."""

        logger.info(
            "Searching S2 data from {} to {} for {}".format(
                self.req_params.datestart, self.req_params.dateend, self.aoi.name
            )
        )
        bbox = list(self.aoi.geometry.bounds)
        items = EarthSearch(
            datestart=self.req_params.datestart,
            dateend=self.req_params.dateend,
            bbox=bbox,
            collection=EarthSearchCollection.SENTINEL2_C1_L2A,
        ).get_items()

        self.s2_items = [AWSSentinel2Item(item) for item in items]
        self.sort_s2_items()

    def get_quality_info(self):
        """Get quality information for Sentinel-2 items."""

        # Check that s2_items are available
        if not self.s2_items:
            logger.info("No Sentinel-2 items available.")
            return None

        logger.info("Computing S2 quality information...")
        if self.multiprocessing is not None:
            self.s2_items = _multiprocess_get_scl_data(
                self.s2_items,
                self.aoi,
                self.req_params.qi_evaluation_scale,
                self.multiprocessing,
            )
        else:
            for s2_item in self.s2_items:
                s2_item.get_item_data(
                    self.aoi, [S2Band.SCL], self.req_params.qi_evaluation_scale
                )
                s2_item.add_class_percentages()

        self.create_quality_information()

    def get_s2_data(self):
        """Get Sentinel-2 data."""
        # Check that s2_items are available
        if not self.check_s2_items_exist():
            return None

        self.sort_s2_items()
        logger.info(f"Retrieving S2 data from {len(self.s2_items)} products...")
        if self.multiprocessing is not None:
            self.s2_items = _multiprocess_get_item_s2_data(
                self.s2_items,
                self.aoi,
                self.req_params.bands,
                self.req_params.target_gsd,
                self.multiprocessing,
            )

        else:
            for s2_item in self.s2_items:
                logger.info("Get data for item {}".format(s2_item.metadata.assetid))

                # Get band data
                s2_item.get_item_data(
                    self.aoi,
                    self.req_params.bands,
                    self.req_params.target_gsd,
                )

                # Get observation geometry
                s2_item.get_observation_geometry()


class AWSSentinel2Metadata(Sentinel2Metadata):
    """Class to handle Sentinel-2 metadata from AWS Open data registry.

    Attributes:
    -----------
    profiles: Dict[S2Band, dict]
        Dictionary to store all profiles for bands
    """

    def __init__(self, item: Item):
        """Initialize AWSSentinel2Metadata object.

        Parameters:
        -----------
        item: Item
            Item object from pystac_client
        """

        s2_product_id = item.properties["s2:product_uri"].replace(".SAFE", "")
        date = pd.to_datetime(
            datetime.datetime.strptime(s2_product_id.split("_")[2], "%Y%m%dT%H%M%S")
        )
        time = date
        tileid = "{}{}{}".format(
            item.properties["mgrs:utm_zone"],
            item.properties["mgrs:latitude_band"],
            item.properties["mgrs:grid_square"],
        )
        assetid = item.id
        productid = s2_product_id
        try:
            projection = "EPSG:{}".format(item.properties["proj:epsg"])  # pystac<1.12.0
        except KeyError:
            try:
                projection = item.properties["proj:code"]  # pystac >= 1.12.0
            except KeyError as e:
                logger.error("Projection information not found. Check pySTAC version.")
                raise e

        datasource = DataSource.AWS

        self.profiles: Dict[S2Band, dict] = {}  # to store all profiles for bands
        super().__init__(time, tileid, assetid, productid, projection, datasource)

    def get_observation_geometry(self, item: Item):
        """Get observation geometry for Sentinel-2 item.

        Parameters:
        -----------
        item: Item
            Item object from pystac_client
        """

        self.observation_geometry = get_observation_geometry(item)

    def get_reference_band(self, target_resolution: float) -> S2Band:
        spatial_resolutions = {
            band: abs(profile["transform"][0])
            for band, profile in self.profiles.items()
        }
        # Get first band with spatial resolution equal to target resolution
        reference_band = [
            band
            for band, resolution in spatial_resolutions.items()
            if resolution == target_resolution
        ]

        if not reference_band:
            # Get band with smallest spatial resolution
            reference_band = min(spatial_resolutions, key=spatial_resolutions.get)
        else:
            reference_band = reference_band[0]

        return reference_band


class AWSSentinel2Item(Sentinel2Item):
    """Class to handle Sentinel-2 item from AWS Open data registry.

    Attributes:
    -----------
    source_item: Item
        Item object from pystac_client
    """

    def __init__(self, item: Item):
        """Initialize AWSSentinel2Item object.

        Parameters:
        -----------
        item: Item
            Item object from pystac_client
        """

        self.source_item = item
        super().__init__(AWSSentinel2Metadata(item))

    def get_observation_geometry(self):
        """Get observation geometry for Sentinel-2 item."""
        self.metadata.get_observation_geometry(self.source_item)

    def get_band_data(
        self,
        aoi: AOI,
        band: S2Band,
    ):
        """Get band data for Sentinel-2 item.

        Parameters:
        -----------
        aoi: AOI
            Area of interest
        band: S2Band
            Sentinel-2 band

        """
        DEFAULT_BUFFER = 100  # meters
        band_aws = band.to_aws()
        aoi_geometry_data_crs = transform_crs(
            aoi.geometry, aoi.geometry_crs, self.metadata.projection
        )
        bbox_data_crs = list(aoi_geometry_data_crs.bounds)

        band_metadata = self.source_item.assets[band_aws].extra_fields["raster:bands"][
            0
        ]

        # Currently buffer used for all bands, even though certain bands might not
        # resampling. Otherwise, the data dimensions might not match with
        # GEE data source. Occationally there was one pixeld difference in x or y dim
        # with certain polygons.

        # spatial_resolution = band_metadata["spatial_resolution"]
        # if spatial_resolution != target_resolution:  # resampling needed, use buffer
        buffer = DEFAULT_BUFFER
        bbox_data_crs = expand_bounds(bbox_data_crs, buffer)

        # # Transform aoi to pixel coordinates/window
        data_transform = rasterio.transform.Affine(
            *self.source_item.assets[band_aws].extra_fields["proj:transform"]
        )

        window = (
            rasterio.windows.from_bounds(*bbox_data_crs, data_transform)
            .round_offsets()
            .round_lengths()
        )
        # Get windowed data
        file_url = self.source_item.assets[band_aws].href
        with rasterio.open(file_url) as src:
            profile = src.profile

            if band == S2Band.SCL:
                band_data = src.read(1, window=window, boundless=True)
            else:
                # Reflectance transformation and apply offset if not applied
                # Not applied necessarily in the old collection
                if "earthsearch:boa_offset_applied" in self.source_item.properties:
                    offset_applied = self.source_item.properties[
                        "earthsearch:boa_offset_applied"
                    ]
                else:
                    offset_applied = False
                offset = band_metadata["offset"] if not offset_applied else 0
                scale = band_metadata["scale"]
                band_data = src.read(1, window=window, boundless=True) * scale + offset

            # Form a new rasterio dataset
            transform = rasterio.windows.transform(window, profile["transform"])
            height = band_data.shape[-2]
            width = band_data.shape[-1]

            new_profile = profile.copy()
            new_profile.update(
                transform=transform,
                driver="GTiff",
                height=height,
                width=width,
                dtype=str(band_data.dtype),
                nodata=SCL_NODATA if band == S2Band.SCL else SPECTRAL_BAND_NO_DATA,
            )

        self.metadata.profiles[band] = new_profile
        self.data[band] = band_data

    def get_item_data(
        self,
        aoi: AOI,
        bands: List[S2Band],
        target_resolution: float,
    ):
        """Get data for all bands for Sentinel-2 item.

        Parameters:
        -----------
        aoi: AOI
            Area of interest
        bands: List[S2Band]
            List of Sentinel-2 bands
        target_resolution: float
            Target resolution
        """

        for band in bands:
            self.get_band_data(aoi, band)

        aoi_geometry_item_crs = transform_crs(
            aoi.geometry, aoi.geometry_crs, self.metadata.projection
        )
        # Resample all bands to the same resolution and reproject to same shape
        reference_band = self.metadata.get_reference_band(target_resolution)
        reference_profile = self.metadata.profiles[reference_band]

        for band in bands:
            if (
                self.metadata.profiles[band]["transform"]
                != reference_profile["transform"]
            ):
                src_profile = self.metadata.profiles[band]
                src_data = self.data[band]
                # Don't use directly the reference profile, since it might have
                # different data type than the reprojected band
                new_profile = src_profile.copy()
                new_profile.update(
                    transform=reference_profile["transform"],
                    driver="GTiff",
                    height=reference_profile["height"],
                    width=reference_profile["width"],
                )
                resampling = (
                    rasterio.enums.Resampling.nearest
                    if band == S2Band.SCL
                    else rasterio.enums.Resampling.bilinear
                )
                reproj_data = reproject_data_to_profile(
                    src_data, src_profile, new_profile, resampling
                )

                self.data[band] = reproj_data
                self.metadata.profiles[band] = new_profile

            # Clip data to AOI
            no_data = SCL_NODATA if band == S2Band.SCL else SPECTRAL_BAND_NO_DATA
            band_data = self.data[band]
            profile = self.metadata.profiles[band]
            with MemoryFile() as memfile:
                with memfile.open(**profile) as dataset:
                    dataset.write(band_data, 1)

                band_data, new_profile = mask_raster(
                    memfile, aoi_geometry_item_crs, no_data=no_data
                )

            self.data[band] = band_data
            self.metadata.profiles[band] = new_profile
            self.create_coordinates(band)

    def add_class_percentages(self):
        """Add class percentages for Sentinel-2 item.

        Class percentages are calculated based on the SCL band data.
        """

        # Check that SCL data is available
        if S2Band.SCL not in self.data:
            raise ValueError("No SCL data available for class percentage calculation")
        else:
            scl_data = self.data[S2Band.SCL]
            # Sometimes SCL image is faulty and doesn't contain data at the
            # area of interest. Set class percentages to nan in this case.
            if scl_data.size == 0:
                class_percentages = {scl_class.name: np.nan for scl_class in SCLClass}
            else:
                num_of_aoi_pixels = np.sum(scl_data != SCL_NODATA)
                class_percentages = {}
                for scl_class in SCLClass:
                    class_percentage = (
                        np.sum(scl_data == scl_class.value) / num_of_aoi_pixels
                    )
                    class_percentages[scl_class.name] = class_percentage
            self.metadata.class_percentages = class_percentages

    def create_coordinates(self, band: S2Band):
        """Create coordinates for Sentinel-2 item.

        Parameters:
        -----------
        band: S2Band
            Sentinel-2 band
        """

        profile = self.metadata.profiles[band]
        x_coords, y_coords = coordinate_arrays_from_profile(profile)

        # translate to pixel center coordinates for netcdf/xarray dataset
        dx = profile["transform"].a
        dy = profile["transform"].e
        x_coords_center = x_coords + dx / 2
        y_coords_center = y_coords + dy / 2
        self.set_coordinates(band, x_coords_center, y_coords_center)


def get_xml_metadata(item: Item) -> Dict:
    """Get XML metadata for Sentinel-2 item.

    Parameters:
    -----------
    item: Item
        Item object from pystac_client

    Returns:
    --------
    metadata: Dict
        XML metadata
    """

    with urllib.request.urlopen(item.assets["granule_metadata"].href) as url:
        metadata = xmltodict.parse(url.read().decode())
    metadata = metadata.popitem()[1]
    return metadata


def get_observation_geometry(item: Item) -> Sentinel2ObservationGeometry:
    """Get observation geometry for Sentinel-2 item.

    Parameters:
    -----------
    item: Item
        Item object from pystac_client

    Returns:
    --------
    observation_geometry: Sentinel2ObservationGeometry
        Observation geometry
    """

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

    observation_geometry = Sentinel2ObservationGeometry(
        sunaz,
        sunzen,
        viewaz,
        viewzen,
    )
    return observation_geometry


# functions for multiprocessing
def _get_s2_data_single(
    s2_item: AWSSentinel2Item, aoi: AOI, bands: List[S2Band], target_resolution: float
):
    logger.info("Get data for item {}".format(s2_item.metadata.assetid))

    # Get band data
    s2_item.get_item_data(
        aoi,
        bands,
        target_resolution,
    )

    # Get observation geometry
    s2_item.get_observation_geometry()

    return s2_item


def _get_scl_data_single(s2_item: AWSSentinel2Item, aoi: AOI, target_resolution: float):
    s2_item.get_item_data(
        aoi,
        [S2Band.SCL],
        target_resolution,
    )
    s2_item.add_class_percentages()
    return s2_item


def _multiprocess_get_item_s2_data(
    s2_items: List[AWSSentinel2Item],
    aoi: AOI,
    bands: List[S2Band],
    target_resolution: float,
    processes: int,
) -> List[AWSSentinel2Item]:
    multiprocess_s2_items = [
        (s2_item, aoi, bands, target_resolution) for s2_item in s2_items
    ]

    with Pool(processes) as p:
        results = p.starmap(_get_s2_data_single, multiprocess_s2_items)
    return results


def _multiprocess_get_scl_data(
    s2_items: List[AWSSentinel2Item], aoi: AOI, target_resolution: float, processes: int
) -> List[AWSSentinel2Item]:
    multiprocess_s2_items = [(s2_item, aoi, target_resolution) for s2_item in s2_items]

    with Pool(processes) as p:
        results = p.starmap(_get_scl_data_single, multiprocess_s2_items)
    return results


def split_time_range(
    datestart: pd.Timestamp, dateend: pd.Timestamp
) -> List[Tuple[pd.Timestamp, pd.Timestamp]]:
    """Split time range to half year time ranges.

    Parameters:
    -----------
    datestart: pd.Timestamp
        Start date
    dateend: pd.Timestamp
        End date

    Returns:
    --------
    time_ranges: List[Tuple[pd.Timestamp, pd.Timestamp]]
        List of time ranges

    """
    time_ranges = []
    current_start = datestart

    while current_start < dateend:
        next_end = current_start + pd.Timedelta(days=91)  # Approx. 3 months
        if next_end > dateend:
            next_end = dateend
        time_ranges.append((current_start, next_end))
        current_start = next_end

    return time_ranges


# To be deprecated functions


def search_s2_cogs(aoi: AOI, req_params: Sentinel2RequestParams) -> List[Item]:
    """Search for Sentinel-2 items from AWS Open data registry.

    Parameters:
    -----------
    aoi: AOI
        Area of interest
    req_params: Sentinel2RequestParams
        Request parameters

    Returns:
    --------
    items: List[Item]
        List of Sentinel-2 items
    """

    DEPRECATION_WARNING_TEXT = (
        "Function search_s2_cogs() will be deprecated and removed at some point."
        "Use AWSSentinel2DataCollection, and"
        "AWSSentinel2DataCollection.search_s2_items() instead."
    )
    logger.warning(DEPRECATION_WARNING_TEXT)
    warn(DEPRECATION_WARNING_TEXT, DeprecationWarning, stacklevel=2)

    data_collection = AWSSentinel2DataCollection(aoi, req_params)
    data_collection.search_s2_items()
    items = [s2_item.source_item for s2_item in data_collection.s2_items]
    return items


def cog_get_s2_quality_info(
    aoi: AOI, req_params: Sentinel2RequestParams, items: List[Item]
) -> pd.DataFrame:
    """Get quality information for Sentinel-2 items.

    Parameters:
    -----------
    aoi: AOI
        Area of interest
    req_params: Sentinel2RequestParams
        Request parameters
    items: List[Item]
        List of Sentinel-2 items

    Returns:
    --------
    qi_df: pd.DataFrame
        Quality information dataframe
    """

    DEPRECATION_WARNING_TEXT = (
        "Function cog_get_s2_quality_info() will be deprecated and removed at some "
        "point. Use AWSSentinel2DataCollection and "
        "AWSSentinel2DataCollection.get_quality_info() or function get_s2_qi_and_data()"
        "in satellitetools.common.wrappers instead."
    )
    logger.warning(DEPRECATION_WARNING_TEXT)
    warn(DEPRECATION_WARNING_TEXT, DeprecationWarning, stacklevel=2)

    data_collection = AWSSentinel2DataCollection(aoi, req_params)
    data_collection.s2_items = [AWSSentinel2Item(item) for item in items]
    data_collection.get_quality_info()
    return data_collection.quality_information


def cog_get_s2_band_data(
    aoi: AOI,
    req_params: Sentinel2RequestParams,
    items: List[Item],
    qi_df: pd.DataFrame,
) -> pd.DataFrame:
    """Get Sentinel-2 data.

    Parameters:
    -----------
    aoi: AOI
        Area of interest
    req_params: Sentinel2RequestParams
        Request parameters
    items: List[Item]
        List of Sentinel-2 items
    qi_df: pd.DataFrame
        Quality information dataframe

    Returns:
    --------
    xr_dataset: pd.DataFrame
        Sentinel-2 data as xarray dataset
    """
    DEPRECATION_WARNING_TEXT = (
        "Function cog_get_s2_band_data() will be deprecated and removed at some point."
        "Use AWSSentinel2DataCollection, and AWSSentinel2DataCollection.get_s2_data() "
        " or function get_s2_qi_and_data() in satellitetools.common.wrappers instead."
    )
    logger.warning(DEPRECATION_WARNING_TEXT)
    warn(DEPRECATION_WARNING_TEXT, DeprecationWarning, stacklevel=2)

    data_collection = AWSSentinel2DataCollection(aoi, req_params)
    data_collection.s2_items = [AWSSentinel2Item(item) for item in items]
    data_collection.quality_information = qi_df
    data_collection.get_s2_data()
    data_collection.data_to_xarray()
    return data_collection.xr_dataset
