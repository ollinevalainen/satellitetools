#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common classes for gee and aws_cog process.
TODO: UPDATE docstrings!!

@author: Olli Nevalainen (Finnish Meteorological Institute)
Created on Tue Mar 16 10:45:05 2021
"""
import sys
from satellitetools.common.sentinel2 import (
    S2_BANDS_COG,
    S2_BANDS_GEE,
)


class RequestParams:
    """S2 data request paramaters.

    Attributes
    ----------
    datestart : str
        Starting date for data request in form "2019-01-01".
    dateend : str
        Starting date for data request in form "2019-12-31".
    bands : list, optional
        List of strings with band name.
        the default is ['B3', 'B4', 'B5',
        'B6', 'B7', 'B8A', 'B11', 'B12'].
    """

    def __init__(self, datestart, dateend, datasource, bands, target_gsd=20):
        """.

        Parameters
        ----------
        datestart : str
            Starting date for data request in form "2019-01-01".
        dateend : str
            Starting date for data request in form "2019-12-31".
        datasource : str
            Source of data. Current options "gee" og "aws_cog".
        bands : list, optional
            List of strings with band name.
            The default is ['B3', 'B4', 'B5',
            'B6', 'B7', 'B8A', 'B11', 'B12'].
        target_gsd : float
            Requested Ground Sampling Distance (GSD). All bands will be resampled to this resolution.
            Default 20m.
        Returns
        -------
        None.

        """

        self.datestart = datestart
        self.dateend = dateend
        self.datasource = datasource  # "gee" or  "aws_cog"
        if bands:
            self.bands = bands
        else:
            if datasource == "gee":
                self.bands = S2_BANDS_GEE
            elif datasource == "aws_cog":
                self.bands = S2_BANDS_COG
            else:
                sys.exit("""Bands not given and unknown datasource.""")
        self.target_gsd = target_gsd


class AOI:
    """Area of interest for area info and data.

    TODO: UPDATE docstring

    Attributes
    ----------
    name : str
        Name of the area.
    geometry : str
        Geometry of the area of interest e.g. from geopandas.
        Currently only polygons tested. The default is None.
    coordinate_list : list, optional
        List of coordinates of a polygon
        (loop should be closed). Computed from geometry if not
        provided. The default is None.
    tile : str, optional
        Tile id as string for the data. Used to keep the data in
        same crs because an area can be in multiple tiles with
        different crs. The default is None.
    qi : pandas dataframe
        Dataframe with quality information about available imagery for the AOI.
        qi is empty at init and can be computed with
        ee_get_s2_quality_info function.
    data : pandas dataframe or xarray
        Dataframe holding data retrieved from GEE. Data can be computed using
        function
        qi is empty at init and can be computed with ee_get_s2_data and
        converted to xarray using s2_data_to_xarray function.

    Methods
    -------
    __init__
    """

    def __init__(self, name, geometry, geometry_crs, tile=None):
        """.
        TODO: UPDATE docstring
        Parameters
        ----------
        name : str
            Name of the area.
        geometry : geometry in wkt, optional
            Geometry of the area of interest e.g. from geopandas.
            Currently only polygons tested. The default is None.
        coordinate_list : list, optional
            List of coordinates of a polygon
            (loop should be closed). Computed from geometry if not
            provided. The default is None.
        tile : str, optional
            Tile id as string for the data. Used to keep the data in
            same crs because an area can be in multiple tiles with
            different crs. The default is None.

        Returns
        -------
        None.

        """

        self.name = name
        self.geometry = geometry
        self.geometry_crs = geometry_crs
        self.tile = tile
