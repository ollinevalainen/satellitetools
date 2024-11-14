#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 09:37:54 2020

@author:
    Olli Nevalainen (olli.nevalainen@fmi.fi), Finnish Meteorological Institute)
"""
from zipfile import ZipFile

import geopandas as gpd
import shapely


def read_kmx_file(kml_or_kmz_file) -> gpd.GeoDataFrame:
    """Read KML or KMZ file.

    Parameters
    ----------
    kml_or_kmz_file : str
        Path to KML or KMZ file.

    Returns
    -------
    gdf : gpd.GeoDataFrame
        GeoDataFrame.

    """
    gpd.io.file.fiona.drvsupport.supported_drivers["KML"] = "rw"

    if kml_or_kmz_file.endswith(".kmz"):
        kmz = ZipFile(kml_or_kmz_file, "r")
        kml = kmz.open("doc.kml", "r")
        gdf = gpd.read_file(kml, driver="KML")
    else:
        gdf = gpd.read_file(kml_or_kmz_file, driver="KML")
    gdf.geometry = gdf.geometry.map(
        lambda polygon: shapely.ops.transform(lambda x, y, z: (x, y), polygon)
    )
    return gdf
