#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:05:36 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""
from typing import List, Optional, Tuple

import pandas as pd
import xarray as xr

from satellitetools.aws import AWSSentinel2DataCollection
from satellitetools.common.classes import AOI
from satellitetools.common.sentinel2 import (
    S2_FILTER1,
    DataSource,
    SCLClass,
    Sentinel2RequestParams,
)
from satellitetools.gee import GEESentinel2DataCollection


def get_s2_qi_and_data(
    aoi: AOI,
    req_params: Sentinel2RequestParams,
    qi_threshold: Optional[float] = 0.02,
    qi_filter: Optional[List[SCLClass]] = S2_FILTER1,
) -> Tuple[pd.DataFrame, xr.Dataset]:
    """Get Sentinel-2 quality information and data.

    Parameters:
    ----------------
    aoi: AOI
        Area of interest.
    req_params: Sentinel2RequestParams
        Request parameters.
    qi_threshold: Optional[float]
        Quality information threshold.
    qi_filter: Optional[List[SCLClass]]
        Sentinel-2 scene classification filter.

    Returns:
    ----------------
    Tuple[pd.DataFrame, xr.Dataset]
        Tuple of quality information DataFrame and xarray Dataset.
    """
    if req_params.datasource == DataSource.GEE:
        data_collection = GEESentinel2DataCollection(aoi, req_params)

    elif req_params.datasource == DataSource.AWS_COG:
        data_collection = AWSSentinel2DataCollection(aoi, req_params)
        data_collection.search_s2_items()
    else:
        raise ValueError(f"Unknown datasource {req_params.datasource}")

    print("Computing S2 quality information...")
    data_collection.get_quality_info()
    data_collection.filter_s2_items(qi_threshold, qi_filter)

    if len(data_collection.s2_items) == 0:
        print("No new observations for area %s" % aoi.name)
    else:
        data_collection.get_s2_data()
        data_collection.data_to_xarray()

    return data_collection.quality_information, data_collection.xr_dataset
