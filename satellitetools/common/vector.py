#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:23:20 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""
from typing import Tuple

import numpy as np
import pyproj
from shapely.geometry import Polygon
from shapely.ops import transform


def transform_crs(aoi_geometry: Polygon, src_crs: str, dst_crs: str) -> Polygon:
    """Transform area of interest geometry to another coordinate reference system.

    Parameters
    ----------
    aoi_geometry : Polygon
        Geometry of the area of interest.
    src_crs : str
        Source coordinate reference system.

    dst_crs : str
        Destination coordinate reference system.

    Returns
    -------
    aoi_geometry_transformed : Polygon

    """
    # Tranform aoi_geometry to raster crs

    project = pyproj.Transformer.from_proj(
        pyproj.Proj(src_crs),  # source coordinate system
        pyproj.Proj(dst_crs),
        always_xy=True,
    )
    aoi_geometry_transformed = transform(project.transform, aoi_geometry)
    return aoi_geometry_transformed


def expand_bounds(bbox_coordinates: list, amount: float) -> list:
    """Expand bounding box coordinates by amount.

    Parameters
    ----------
    bbox_coordinates : list
        Bounding box coordinates [xmin, ymin, xmax, ymax].
    amount : float
        Amount to expand the bounding box.

    Returns
    -------
    bbox_coordinates : list
        Expanded bounding box coordinates [xmin, ymin, xmax, ymax].

    """
    bbox_coordinates[0] = bbox_coordinates[0] - amount
    bbox_coordinates[1] = bbox_coordinates[1] - amount
    bbox_coordinates[2] = bbox_coordinates[2] + amount
    bbox_coordinates[3] = bbox_coordinates[3] + amount
    return bbox_coordinates


def coordinate_arrays_from_profile(profile: dict) -> Tuple[np.ndarray, np.ndarray]:
    """Get coordinate arrays from rasterio profile.

    Parameters
    ----------
    profile : dict
        Rasterio profile.

    Returns
    -------
    xs : np.ndarray
        Array of x coordinates.
    ys : np.ndarray
        Array of y coordinates.

    """

    dx = profile["transform"].a
    dy = profile["transform"].e

    # upperleft corner coordinates
    x_ul = profile["transform"].c
    y_ul = profile["transform"].f

    width = profile["width"]
    height = profile["height"]
    xs = np.array([x_ul + i * dx for i in range(width)])
    ys = np.array([y_ul + i * dy for i in range(height)])

    return xs, ys
