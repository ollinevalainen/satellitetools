#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 16:20:30 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""

# %%
import os
import time

import geopandas as gpd

import satellitetools.aws_cog.aws_cog as aws_cog
import satellitetools.biophys.biophys as biophys
from satellitetools.common.classes import AOI, RequestParams
from satellitetools.common.sentinel2 import S2_BANDS_10_20_COG
from satellitetools.common.wrappers import get_s2_qi_and_data

start_time = time.time()
package_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))
# should also work with geojson and kml/kmz
geomfile = os.path.join(package_dir, "geometry-files/ruukki_blocks_wgs84.shp")
gdf = gpd.read_file(geomfile)

ruukki_block1 = gdf.loc[[0]]
aoi_name = ruukki_block1.iloc[0].Name.replace(" ", "").lower()
geom = ruukki_block1.iloc[0].geometry
crs = ruukki_block1.crs.to_string()

aoi = AOI(aoi_name, geom, crs)
req_params = RequestParams(
    "2023-06-01", "2023-06-10", "aws_cog", bands=S2_BANDS_10_20_COG, target_gsd=10.0
)
# %%
items = aws_cog.search_s2_cogs(aoi, req_params)

# %%
qi_df = aws_cog.cog_get_s2_quality_info(aoi, req_params, items)
dataset = aws_cog.cog_get_s2_band_data(aoi, req_params, items, qi_df)

aoi.qi = qi_df
aoi.data = dataset

end_time = time.time()

print(end_time - start_time)
# %% Alternative using wrapper

datasource = "aws_cog"
qi_threshold = 0.02

request = RequestParams(
    "2020-06-01",
    "2020-06-10",
    datasource=datasource,
    bands=S2_BANDS_10_20_COG,
    target_gsd=10.0,
)
aoi.qi, aoi.data = get_s2_qi_and_data(aoi, request, qi_threshold=qi_threshold)

# %% Test biohys

aoi.data = biophys.run_snap_biophys(aoi.data, "LAI")

# %%
