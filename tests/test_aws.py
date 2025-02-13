import logging
import os
from pathlib import Path

from shapely.geometry import MultiPolygon, Polygon

import satellitetools as sattools

logger = logging.getLogger(__name__)
logger.setLevel(logging.DEBUG)

qvidja_ec_field_geom = [
    [22.3913931, 60.295311],
    [22.3917056, 60.2951721],
    [22.3922131, 60.2949717],
    [22.3927016, 60.2948124],
    [22.3932251, 60.2946874],
    [22.3931117, 60.2946416],
    [22.3926039, 60.2944037],
    [22.3920127, 60.2941585],
    [22.3918447, 60.2940601],
    [22.391413, 60.2937852],
    [22.3908102, 60.2935286],
    [22.390173, 60.2933897],
    [22.389483, 60.2933106],
    [22.3890777, 60.293541],
    [22.3891442, 60.2936358],
    [22.3889863, 60.2940313],
    [22.3892131, 60.2941537],
    [22.3895462, 60.2942468],
    [22.3899066, 60.2944289],
    [22.3903881, 60.2946329],
    [22.3904738, 60.2948121],
    [22.3913931, 60.295311],
]

qvidja_polygon = Polygon(qvidja_ec_field_geom)
test_AOI = sattools.AOI("qvidja_ec", qvidja_polygon, "EPSG:4326")
test_datestart = "2023-06-01"
test_dateend = "2023-06-15"
test_bands = [sattools.S2Band.B4, sattools.S2Band.B8A, sattools.S2Band.SCL]
test_bands_10_20 = sattools.S2Band.get_10m_to_20m_bands()
test_multiprocessing = 4


class TestAWS:
    def test_get_s2_data(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.AWS, test_bands
        )
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(
            test_AOI, req_params
        )
        s2_data_collection.search_s2_items()
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands])

        s2_data_collection.data_to_xarray()

    def test_get_s2_data_with_multiprocessing(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.AWS, test_bands
        )
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(
            test_AOI, req_params, test_multiprocessing
        )
        s2_data_collection.search_s2_items()
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands])

        s2_data_collection.data_to_xarray()

    def test_get_s2_data_2022_with_multiprocessing(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.AWS, test_bands
        )
        req_params.datestart = "2022-06-01"
        req_params.dateend = "2022-06-15"
        req_params.bands = [sattools.S2Band.B4, sattools.S2Band.B8A]
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(
            test_AOI, req_params, test_multiprocessing
        )
        s2_data_collection.search_s2_items()
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        s2_data_collection.data_to_xarray()
        s2_data_collection.xr_dataset.to_netcdf(
            Path(__file__).parent / "test_data" / "test_aws_2022.nc"
        )

    def test_search_s2_data_2022_2023(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.AWS, test_bands
        )
        req_params.datestart = "2021-11-01"
        req_params.dateend = "2022-06-15"
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(
            test_AOI, req_params
        )
        s2_data_collection.search_s2_items()
        assert s2_data_collection.s2_items is not None

    def test_get_s2_data_10_20_bands_20m(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.AWS, test_bands
        )
        req_params.bands = test_bands_10_20
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(
            test_AOI, req_params, os.cpu_count()
        )
        s2_data_collection.search_s2_items()
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands])

        s2_data_collection.data_to_xarray()
        s2_data_collection.xr_dataset.to_netcdf(
            Path(__file__).parent / "test_data" / "test_aws_20m.nc"
        )

    def test_get_s2_data_10_20_bands_10m(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.AWS, test_bands
        )
        req_params.bands = test_bands_10_20
        req_params.target_gsd = 10
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(
            test_AOI, req_params, os.cpu_count()
        )
        s2_data_collection.search_s2_items()
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands])

        s2_data_collection.data_to_xarray()
        s2_data_collection.xr_dataset.to_netcdf(
            Path(__file__).parent / "test_data" / "test_aws_10m.nc"
        )

    def test_get_s2_data_with_multiprocessing_2020(self):
        req_params = sattools.Sentinel2RequestParams(
            "2020-06-01", "2020-06-15", sattools.DataSource.AWS, test_bands
        )
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(
            test_AOI, req_params, test_multiprocessing
        )
        s2_data_collection.search_s2_items()
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands])

    def test_search_for_long_timeperiod_in_tile_buffer_zone(self):
        difficult_polygon = Polygon(
            [
                [23.82573048401551, 60.75855515465485],
                [23.82544917989838, 60.75680426756835],
                [23.82657367536738, 60.75691684707225],
                [23.82691587445285, 60.75850501120518],
                [23.82573048401551, 60.75855515465485],
            ]
        )
        aoi = sattools.AOI("test", difficult_polygon, "EPSG:4326")
        req_params = sattools.Sentinel2RequestParams(
            "2019-01-01", "2021-12-31", sattools.DataSource.AWS, test_bands
        )
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(
            aoi, req_params, test_multiprocessing
        )
        s2_data_collection.search_s2_items()
        assert s2_data_collection.s2_items is not None

    def test_get_scl_data_without_quality_information(self):
        year = 2023
        date_start = f"{year}-06-01"
        date_end = f"{year}-07-01"
        data_source = sattools.DataSource.AWS
        bands = [sattools.S2Band.SCL]
        request = sattools.Sentinel2RequestParams(
            date_start,
            date_end,
            data_source,
            bands=bands,
            target_gsd=20,
        )
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(test_AOI, request)
        s2_data_collection.search_s2_items()
        s2_data_collection.filter_s2_items_by_tile()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in bands])
        s2_data_collection.data_to_xarray()

    def test_product_with_missing_data(self):
        aoi_coords = [
            [25.874837428345387, 64.789844062914838],
            [25.87563257650196, 64.790851981054146],
            [25.875984787908642, 64.790814842307171],
            [25.875167562154864, 64.789786508355235],
            [25.874837428345387, 64.789844062914838],
        ]
        aoi_polygon = Polygon(aoi_coords)
        aoi = sattools.AOI("problematic", aoi_polygon, "EPSG:4326")
        date_start = "2020-08-20"
        date_end = "2020-09-10"
        data_source = sattools.DataSource.AWS
        bands = sattools.S2Band.get_10m_to_20m_bands()
        request = sattools.Sentinel2RequestParams(
            date_start,
            date_end,
            data_source,
            bands=bands,
            target_gsd=20,
        )
        # Removed multiprocessing as there's potentially an issue in the AWS data
        # that occur with this specific AOI and time period. Processes seem to be
        # reading a same image occasionally causing  error
        # rasterio._err.CPLE_AppDefinedError: SCL.tif, band 1: IReadBlock failed at...
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(aoi, request)
        s2_data_collection.search_s2_items()
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        s2_data_collection.data_to_xarray()

    def test_multipolygon(self):
        qvidja_ca15 = [
            [22.3873751, 60.2854051],
            [22.389999, 60.2870115],
            [22.3904474, 60.2868637],
            [22.3879298, 60.28518369999998],
            [22.3873751, 60.2854051],
        ]
        qvidja_ca6 = [
            [22.3889573, 60.28739560000001],
            [22.3896006, 60.2871823],
            [22.3870542, 60.2855162],
            [22.3865752, 60.2857636],
            [22.3889573, 60.28739560000001],
        ]

        aoi_polygon = MultiPolygon([Polygon(qvidja_ca15), Polygon(qvidja_ca6)])
        aoi = sattools.AOI("multipolygon", aoi_polygon, "EPSG:4326")
        data_source = sattools.DataSource.AWS
        bands = sattools.S2Band.get_10m_to_20m_bands()
        request = sattools.Sentinel2RequestParams(
            test_datestart,
            test_dateend,
            data_source,
            bands=bands,
            target_gsd=10,
        )
        s2_data_collection = sattools.aws.AWSSentinel2DataCollection(aoi, request)
        s2_data_collection.search_s2_items()
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        s2_data_collection.data_to_xarray()
