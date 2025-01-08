import os
import pickle  # noqa
from pathlib import Path

import ee
from shapely.geometry import Polygon

import satellitetools as sattools


def initialize_gee():
    ee_project = os.environ.get("EE_PROJECT_PYTEST")
    ee_service_account = os.environ.get("EE_SERVICE_ACCOUNT_PYTEST")
    ee_service_account_credentials = os.environ.get(
        "EE_SERVICE_ACCOUNT_CREDENTIALS_PYTEST"
    )
    ee_service_account_credentials_file = os.environ.get(
        "EE_SERVICE_ACCOUNT_CREDENTIALS_FILE_PYTEST"
    )

    if ee_project is not None:
        ee.Initialize(project=ee_project)
    elif ee_service_account is not None and ee_service_account_credentials is not None:
        credentials = ee.ServiceAccountCredentials(
            ee_service_account, key_data=ee_service_account_credentials
        )
        ee.Initialize(credentials=credentials)
    elif (
        ee_service_account_credentials_file is not None
        and ee_service_account_credentials_file is not None
    ):
        credentials = ee.ServiceAccountCredentials(
            ee_service_account, key_file=ee_service_account_credentials_file
        )
        ee.Initialize(credentials=credentials)
    else:
        # Try to initialize with default credentials
        ee.Initialize()


initialize_gee()

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

test_polygon_lake = Polygon(
    [
        [23.76072665998865, 60.78877854773374],
        [23.760445355871518, 60.78702931064117],
        [23.76156985134052, 60.7871417840537],
        [23.76191205042599, 60.788728451538454],
        [23.76072665998865, 60.78877854773374],
    ]
)
test_lake_AOI = sattools.AOI("lake", test_polygon_lake, "EPSG:4326")


class TestGEE:
    def test_get_s2_data(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.GEE, test_bands
        )
        s2_data_collection = sattools.gee.GEESentinel2DataCollection(
            test_AOI, req_params
        )
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands])
        s2_data_collection.data_to_xarray()

    def test_get_s2_data_10_20_bands_20m(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.GEE, test_bands
        )
        req_params.bands = test_bands_10_20
        s2_data_collection = sattools.gee.GEESentinel2DataCollection(
            test_AOI, req_params
        )
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands_10_20])
        s2_data_collection.data_to_xarray()
        s2_data_collection.xr_dataset.to_netcdf(
            Path(__file__).parent / "test_data" / "test_gee_20m.nc"
        )

    def test_get_s2_data_2022(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.GEE, test_bands
        )
        req_params.datestart = "2022-06-01"
        req_params.dateend = "2022-06-15"
        req_params.bands = [sattools.S2Band.B4, sattools.S2Band.B8A]
        s2_data_collection = sattools.gee.GEESentinel2DataCollection(
            test_AOI, req_params
        )
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        s2_data_collection.data_to_xarray()
        s2_data_collection.xr_dataset.to_netcdf(
            Path(__file__).parent / "test_data" / "test_gee_2022.nc"
        )

    def test_get_s2_data_10_20_bands_10m(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.GEE, test_bands
        )
        req_params.bands = test_bands_10_20
        req_params.target_gsd = 10
        s2_data_collection = sattools.gee.GEESentinel2DataCollection(
            test_AOI, req_params
        )
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands_10_20])
        s2_data_collection.data_to_xarray()
        s2_data_collection.xr_dataset.to_netcdf(
            Path(__file__).parent / "test_data" / "test_gee_10m.nc"
        )

    def test_get_s2_data_lake(self):
        req_params = sattools.Sentinel2RequestParams(
            "2021-01-01", "2021-12-31", sattools.DataSource.GEE, test_bands
        )
        s2_data_collection = sattools.gee.GEESentinel2DataCollection(
            test_lake_AOI, req_params
        )
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in test_bands])
        s2_data_collection.data_to_xarray()

    def test_search_s2_items(self):
        req_params = sattools.Sentinel2RequestParams(
            test_datestart, test_dateend, sattools.DataSource.GEE, test_bands
        )
        s2_data_collection = sattools.gee.GEESentinel2DataCollection(
            test_AOI, req_params
        )
        s2_data_collection.search_s2_items()
        assert s2_data_collection.s2_items is not None
        # Save s2_data_collection to file
        with open(
            Path(__file__).parent / "test_data" / "s2_data_collection_gee.pkl", "wb"
        ) as f:
            pickle.dump(s2_data_collection, f)

    def test_get_scl_data_without_quality_information(self):
        year = 2022
        date_start = f"{year}-01-01"
        date_end = f"{year}-12-31"
        data_source = sattools.DataSource.GEE
        bands = [sattools.S2Band.SCL]
        request = sattools.Sentinel2RequestParams(
            date_start,
            date_end,
            data_source,
            bands=bands,
            target_gsd=20,
        )
        s2_data_collection = sattools.gee.GEESentinel2DataCollection(test_AOI, request)
        s2_data_collection.search_s2_items()
        s2_data_collection.filter_s2_items_by_tile()
        s2_data_collection.get_s2_data()
        assert s2_data_collection.s2_items is not None
        assert all([b in s2_data_collection.s2_items[0].data for b in bands])
        s2_data_collection.data_to_xarray()

    def test_aoi_with_y_len_one_pixel(self):
        aoi_coords = [
            (0.7928222017683817, 42.63188925867777),
            (0.7932474148375968, 42.63183255500702),
            (0.7935161231648058, 42.63207447646401),
            (0.7933808033796685, 42.63216285078755),
            (0.7931923944390995, 42.63208507482328),
            (0.7929224148153224, 42.63210705264266),
            (0.7927627009200733, 42.6320269199325),
            (0.7924680675306872, 42.63195239434005),
            (0.7928222017683817, 42.63188925867777),
        ]
        aoi_polygon = Polygon(aoi_coords)
        aoi = sattools.AOI("problematic", aoi_polygon, "EPSG:4326")
        date_start = "2018-01-01"
        date_end = "2018-12-31"
        data_source = sattools.DataSource.GEE
        bands = sattools.S2Band.get_10m_to_20m_bands()
        request = sattools.Sentinel2RequestParams(
            date_start,
            date_end,
            data_source,
            bands=bands,
            target_gsd=20,
        )
        s2_data_collection = sattools.gee.GEESentinel2DataCollection(aoi, request)
        s2_data_collection.get_quality_info()
        s2_data_collection.filter_s2_items()
        s2_data_collection.get_s2_data()
        s2_data_collection.data_to_xarray()
