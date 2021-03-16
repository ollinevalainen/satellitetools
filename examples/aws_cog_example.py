#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 16:20:30 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""


#%%
import satellitetools.aws_cog as aws_cog
from satellitetools.common.classes import AOI, RequestParams
from satellitetools.common.sentinel2 import S2_BANDS_10_20_COG
import geopandas as gpd
import time

start_time = time.time()
ruukki_df = gpd.read_file(
    "/Users/olli/python/satellitetools/geometry-files/ruukki_blocks_wgs84.shp"
)
ruukki_block1 = ruukki_df.loc[[0]]
aoi_name = ruukki_block1.iloc[0].Name.replace(" ", "").lower()
geom = ruukki_block1.iloc[0].geometry
crs = ruukki_block1.crs.to_string()

aoi = AOI(aoi_name, geom, crs)
req_params = RequestParams(
    "2020-06-01", "2020-06-10", "aws_cog", bands=S2_BANDS_10_20_COG, target_gsd=10.0
)
#%%
items = aws_cog.search_s2_cogs(aoi, req_params)
qi_df = aws_cog.cog_get_s2_quality_info(aoi, req_params, items)
dataset = aws_cog.cog_get_s2_band_data(aoi, req_params, items, qi_df)

aoi.qi = qi_df
aoi.data = dataset

end_time = time.time()

print(end_time - start_time)
#%% Alternative using wrapper
from satellitetools.common.wrappers import get_s2_qi_and_data

datasource = "aws_cog"
qi_threshold = 0.02

request = RequestParams(
    "2020-06-01",
    "2020-06-10",
    datasource=datasource,
    bands=S2_BANDS_10_20_COG,
    target_gd=10.0,
)
aoi.qi, aoi.data = get_s2_qi_and_data(aoi, request, qi_threshold=qi_threshold)

#%% Test biohys
import satellitetools.biophys as biophys

aoi.data = biophys.run_snap_biophys(aoi.data, "LAI")
