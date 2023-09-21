#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:30:39 2020

@author: Olli Nevalainen (olli.nevalainen@fmi.fi),
 Finnish Meteorological Institute)

Example of GEE usage
"""
# %%
import os

import geopandas as gpd

import satellitetools.biophys.biophys as biophys
import satellitetools.gee.gee as gee
from satellitetools.common.classes import AOI, RequestParams
from satellitetools.common.sentinel2 import S2_BANDS_10_20_GEE
from satellitetools.common.wrappers import get_s2_qi_and_data
from satellitetools.common.xrtools import xr_dataset_to_timeseries

# %% Get data
package_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# should also work with geojson and kml/kmz
geomfile = os.path.join(package_dir, "geometry-files/ruukki_blocks_wgs84.shp")

datestart = "2019-06-01"
dateend = "2019-06-30"
qi_threshold = 0.02

# the example geometry file has eight different polygons
gdf = gpd.read_file(geomfile)

request = RequestParams(datestart, dateend, datasource="gee", bands=S2_BANDS_10_20_GEE)


areas = []
for block in gdf.itertuples(index=True, name="Block"):
    aoi = AOI(block.Name, block.geometry, gdf.crs.to_string())
    areas.append(aoi)

print("Computing S2 quality information...")
qi_df_dict = gee.ee_get_s2_quality_info(areas, request)
print("Retrieving S2 data...")
data_dict = gee.ee_get_s2_data(areas, request, qi_df_dict, qi_threshold=qi_threshold)
for aoi in areas:
    aoi.qi = qi_df_dict[aoi.name]
    aoi.data = data_dict[aoi.name]

# %% If you get errors such as "Computation timed out." This means that you are
# trying to get too much data (too many AOIs or too long time range) from GEE at
# the same time. Solution to this is to get data for AOIs sequentially, like this:
geomfile = "/Users/olli/python/satellitetools/geometry-files/ruukki_blocks_wgs84.shp"

datestart = "2019-06-01"
dateend = "2019-06-30"
qi_threshold = 0.02

# the example geometry file has eight different polygons
gdf = gpd.read_file(geomfile)

request = RequestParams(datestart, dateend, datasource="gee", bands=S2_BANDS_10_20_GEE)

areas = []
for block in gdf.itertuples(index=True, name="Block"):
    aoi = AOI(block.Name, block.geometry, gdf.crs.to_string())
    print("Computing S2 quality information...")
    qi_df_dict = gee.ee_get_s2_quality_info(aoi, request)
    print("Retrieving S2 data...")
    data_dict = gee.ee_get_s2_data(aoi, request, qi_df_dict, qi_threshold=qi_threshold)
    aoi.qi = qi_df_dict[aoi.name]
    aoi.data = data_dict[aoi.name]
    areas.append(aoi)

# %% Compute NDVI, LAI, FAPAR
for aoi in areas:
    # compute NDVI, fapar and LAI
    aoi.data = biophys.compute_ndvi(aoi.data)
    aoi.data = biophys.run_snap_biophys(aoi.data, "LAI")
    aoi.data = biophys.run_snap_biophys(aoi.data, "FAPAR")

# %% Compute timeseries from xarrays
timeseries = {}
timeseries_variables = ["lai", "fapar", "ndvi"]
for aoi in areas:
    timeseries[aoi.name] = xr_dataset_to_timeseries(aoi.data, timeseries_variables)

# %% save xarrays as netcdf
out_dir = ""
for aoi in areas:
    aoi.data.to_netcdf(os.path.join(out_dir, aoi.name + "example_netcdf.nc"))
    # timeseries

# %% Alternative way to retrieve data using the wrapper, works only for single AOI instance
package_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), os.pardir))

# should also work with geojson and kml/kmz
geomfile = os.path.join(package_dir, "geometry-files/ruukki_blocks_wgs84.shp")

datestart = "2019-06-01"
dateend = "2019-06-30"
qi_threshold = 0.02

# the example geometry file has eight different polygons
gdf = gpd.read_file(geomfile)

request = RequestParams(datestart, dateend, datasource="gee", bands=S2_BANDS_10_20_GEE)


areas = []
for block in gdf.itertuples(index=True, name="Block"):
    aoi = AOI(block.Name, block.geometry, gdf.crs.to_string())
    areas.append(aoi)

aoi = areas[0]
datasource = "gee"
request = RequestParams(
    datestart, dateend, datasource=datasource, bands=S2_BANDS_10_20_GEE
)
aoi.qi, aoi.data = get_s2_qi_and_data(aoi, request, qi_threshold=qi_threshold)

# %%
