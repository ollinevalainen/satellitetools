#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 16:30:39 2020

@author: Olli Nevalainen (olli.nevalainen@fmi.fi),
 Finnish Meteorological Institute)

Example of GEE usage
"""
import satellitetools.gee as gee
import satellitetools.biophys as biophys
from satellitetools.common.classes import AOI, RequestParams
from satellitetools.common.sentinel2 import S2_BANDS_10_20_GEE
import geopandas as gpd
import os

#%%

geomfile = "geometry-files/ruukki_blocks_wgs84.shp"

datestart = "2019-01-01"
dateend = "2019-12-31"
qi_threshold = 0.02

# the example geometry file has eight different polygons
df = gpd.read_file(geomfile)

request = RequestParams(datestart, dateend, S2_BANDS_10_20_GEE)


areas = []
for block in df.itertuples(index=True, name="Block"):
    aoi = AOI(block.Name, block.geometry, df.crs.to_string())
    areas.append(aoi)

print("Computing S2 quality information...")
gee.ee_get_s2_quality_info(areas, request)
print("Retrieving S2 data...")
gee.ee_get_s2_data(areas, request, qi_threshold=qi_threshold)

#%%
for aoi in areas:
    # compute NDVI, fapar and LAI
    aoi.data = biophys.compute_ndvi(aoi.data)
    aoi.data = biophys.run_snap_biophys(aoi.data, "LAI")
    aoi.data = biophys.run_snap_biophys(aoi.data, "FAPAR")

#%%
timeseries = {}
timeseries_variables = ["lai", "fapar", "ndvi"]
for aoi in areas:
    aoi.data.to_netcdf(os.path.join(out_dir, aoi.name + "example_netcdf.nc"))
    # timeseries
    timeseries[aoi.name] = gee.xr_dataset_to_timeseries(aoi.data, timeseries_variables)
