#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common classes for gee and aws_cog process.
TODO: UPDATE docstrings!!

@author: Olli Nevalainen (Finnish Meteorological Institute)
Created on Tue Mar 16 10:45:05 2021
"""
from enum import Enum
from typing import Union

from shapely.geometry import Polygon


class DataSource(str, Enum):
    """Data source for the data request."""

    GEE = "gee"
    AWS_COG = "aws_cog"


class AOI:
    """Area of interest (AOI) class.

    Attributes:
    name: str
        Name of the AOI.
    geometry: Union[str, Polygon]
        Geometry of the AOI. Can be either a WKT string or a shapely Polygon.
    geometry_crs: str
        Coordinate reference system of the geometry.
    tile: str
        Sentinel-2 tile ID.

    """

    def __init__(
        self,
        name: str,
        geometry: Union[str, Polygon],
        geometry_crs: str,
        tile: str = None,
    ):
        """Initialize the AOI class.

        Parameters:
        name: str
            Name of the AOI.
        geometry: Union[str, Polygon]
            Geometry of the AOI. Can be either a WKT string or a shapely Polygon.
        geometry_crs: str
            Coordinate reference system of the geometry.
        tile: str
            Sentinel-2 tile ID.
        """

        self.name = name
        self.geometry = geometry
        self.geometry_crs = geometry_crs
        self.tile = tile

    def __repr__(self) -> str:
        return (
            f"AOI(name={self.name}, geometry={self.geometry}, "
            f"geometry_crs={self.geometry_crs}, tile={self.tile})"
        )
