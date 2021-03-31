#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 16 11:05:36 2021

@author: Olli Nevalainen (Finnish Meteorological Institute)
"""
import sys
from satellitetools.common.sentinel2 import S2_FILTER1


def get_s2_qi_and_data(aoi, req_params, qi_threshold=None, qi_filter=S2_FILTER1):
    if req_params.datasource == "gee":
        import satellitetools.gee as gee

        print("Computing S2 quality information...")
        qi_df_dict = gee.ee_get_s2_quality_info(aoi, req_params)
        qi_df = qi_df_dict[aoi.name]
        if qi_df is None or qi_df.empty:
            print("No new observations for area %s" % aoi.name)
            dataset = None
        else:
            print("Retrieving S2 data...")
            dataset_dict = gee.ee_get_s2_data(
                aoi,
                req_params,
                qi_df_dict,
                qi_threshold=qi_threshold,
                qi_filter=qi_filter,
            )
            if dataset_dict is None:
                dataset = None
            else:
                dataset = dataset_dict[aoi.name]

    elif req_params.datasource == "aws_cog":
        import satellitetools.aws_cog as aws_cog

        items = aws_cog.search_s2_cogs(aoi, req_params)
        if items is None:
            qi_df = None
        else:
            qi_df = aws_cog.cog_get_s2_quality_info(aoi, req_params, items)

        if qi_df is None or qi_df.empty:
            print("No new observations for area %s" % aoi.name)
            dataset = None
        else:
            print("Retrieving S2 data...")
            dataset = aws_cog.cog_get_s2_band_data(
                aoi,
                req_params,
                items,
                qi_df,
                qi_threshold=qi_threshold,
                qi_filter=qi_filter,
            )
    else:
        sys.exit("Unknown data source %s" % req_params.datasource)

    return qi_df, dataset
