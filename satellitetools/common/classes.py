#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Common classes for gee and aws_cog process.

@author: Olli Nevalainen (Finnish Meteorological Institute)

"""
from typing import Union

try:
    # breaking change introduced in python 3.11
    from enum import StrEnum
except ImportError:
    from enum import Enum

    class StrEnum(str, Enum):
        pass


from shapely.geometry import Polygon


class DataSource(StrEnum):
    """Data source for the data request."""

    GEE = "gee"
    AWS = "aws_cog"


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


# Deprecated
class RequestParams:
    def __init__(self, *args, **kwargs):
        raise DeprecationWarning(
            "RequestParams is deprecated. Replace with Sentinel2RequestParams from"
            "common.sentinel2 with same parameters."
        )
