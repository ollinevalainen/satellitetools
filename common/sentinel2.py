#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 10:53:15 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""

S2_REFL_TRANS = 10000

S2_BANDS_COG = [
    "B01",
    "B02",
    "B03",
    "B04",
    "B05",
    "B06",
    "B07",
    "B08",
    "B8A",
    "B09",
    "B11",
    "B12",
]
# 10-20 m bands
S2_BANDS_10_20_COG = ["B02", "B03", "B04", "B05", "B06", "B07", "B8A", "B11", "B12"]

S2_BANDS_GEE = [b.replace("B0", "B") for b in S2_BANDS_COG]
S2_BANDS_10_20_GEE = [b.replace("B0", "B") for b in S2_BANDS_10_20_COG]

s2_qi_labels = [
    "NODATA",
    "SATURATED_DEFECTIVE",
    "DARK_FEATURE_SHADOW",
    "CLOUD_SHADOW",
    "VEGETATION",
    "NOT_VEGETATED",
    "WATER",
    "UNCLASSIFIED",
    "CLOUD_MEDIUM_PROBA",
    "CLOUD_HIGH_PROBA",
    "THIN_CIRRUS",
    "SNOW_ICE",
]

S2_SCL_CLASSES = [
    "NODATA",
    "SATURATED_DEFECTIVE",
    "DARK_FEATURE_SHADOW",
    "CLOUD_SHADOW",
    "VEGETATION",
    "NOT_VEGETATED",
    "WATER",
    "UNCLASSIFIED",
    "CLOUD_MEDIUM_PROBA",
    "CLOUD_HIGH_PROBA",
    "THIN_CIRRUS",
    "SNOW_ICE",
]

S2_FILTER1 = [
    "NODATA",
    "SATURATED_DEFECTIVE",
    "CLOUD_SHADOW",
    "UNCLASSIFIED",
    "CLOUD_MEDIUM_PROBA",
    "CLOUD_HIGH_PROBA",
    "THIN_CIRRUS",
    "SNOW_ICE",
]


def filter_s2_qi_dataframe(s2_qi_dataframe, qi_thresh, s2_filter=S2_FILTER1):
    """Filter qi dataframe.

    Parameters
    ----------
    s2_qi_dataframe : pandas dataframe
        S2 quality information dataframe (AOI instance qi attribute).
    qi_thresh : float
        Threshold value to filter images based on used qi filter.
        qi filter holds labels of classes whose percentages within the AOI
        is summed. If the sum is larger then the qi_threhold, data will not be
        retrieved for that date/image. The default is 1, meaning all data is
        retrieved.
    s2_filter : list
        List of strings with class labels (of unwanted classes) used to compute qi value,
        see qi_threhold. The default is s2_filter1 = ['NODATA',
              'SATURATED_DEFECTIVE',
              'CLOUD_SHADOW',
              'UNCLASSIFIED',
              'CLOUD_MEDIUM_PROBA',
              'CLOUD_HIGH_PROBA',
              'THIN_CIRRUS',
              'SNOW_ICE'].

    Returns
    -------
    filtered_s2_qi_df : pandas dataframe
        Filtered dataframe.

    """

    # Drop if SCL data is nan
    s2_qi_dataframe = s2_qi_dataframe.dropna(axis=0, subset=S2_SCL_CLASSES)
    filtered_s2_qi_df = s2_qi_dataframe.loc[
        s2_qi_dataframe[s2_filter].sum(axis=1) < qi_thresh
    ]

    return filtered_s2_qi_df
