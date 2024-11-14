import cProfile
import pstats

from shapely.geometry import Polygon

import satellitetools as sattools

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


def main():
    qvidja_polygon = Polygon(qvidja_ec_field_geom)
    test_AOI = sattools.AOI("qvidja_ec", qvidja_polygon, "EPSG:4326")
    test_datestart = "2023-06-01"
    test_dateend = "2023-06-15"
    test_bands = [sattools.S2Band.B4, sattools.S2Band.B8A, sattools.S2Band.SCL]
    req_params = sattools.Sentinel2RequestParams(
        test_datestart, test_dateend, sattools.DataSource.AWS_COG, test_bands
    )
    s2_data_collection = sattools.aws.AWSSentinel2DataCollection(test_AOI, req_params)
    s2_data_collection.search_s2_items()
    s2_data_collection.get_quality_info()
    s2_data_collection.filter_s2_items()
    s2_data_collection.get_s2_data()

    s2_data_collection.data_to_xarray()


if __name__ == "__main__":
    profiler = cProfile.Profile()
    profiler.enable()
    main()
    profiler.disable()
    stats = pstats.Stats(profiler).sort_stats("tottime")
    stats.print_stats()
