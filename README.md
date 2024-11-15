[![DOI](https://zenodo.org/badge/270676132.svg)](https://zenodo.org/badge/latestdoi/270676132)

# satellitetools 🛰️ 
This package provide methods to get Sentinel-2 L2A data from Google Earth Engine 
(https://developers.google.com/earth-engine/) or from cloud-optimized geotiffs 
maintained by Element84 from AWS Open data registry 
(https://registry.opendata.aws/sentinel-2-l2a-cogs/, https://github.com/Element84/earth-search).

It has been mainly used in agricultural projects at the Climate System Research unit at 
the Finnish Meteorological Institute.

The package is at the moment under constant development. 
Apologies for the lack of documentation at the moment. Please feel 
free to contact me in case you want more information about this package.

 Package also includes python implementation of ESA's SNAP Biophysical processor v1 
 which can be used to compute biophysical parameters, such as LAI (Original Java code in
  SNAP here: https://github.com/senbox-org/s2tbx/tree/master/s2tbx-biophysical/src/main/java/org/esa/s2tbx/biophysical).
 ATBD: http://step.esa.int/docs/extra/ATBD_S2ToolBox_L2B_V1.1.pdf.

 ## Major changes from merging develop-2024-update to develop (develop is now the master)! ##
 **Breaking changes:** ⚠️

* Renamed xrtools.py to timeseries.py
* RequestParams class is now Sentinel2RequestParams in sentinel2 submodule
* Restructured files and removed nesting. For example gee imported as satellitetools.gee instead of satellitetools.gee.gee
* The quality information dataframe doesn't have anymore "Date" column, but instead index named "acquisition_time" as UTC aware timestamp (from pd.to_datetime).

**New features:** 🔧

* Refactored and made codes more object-oriented and modular:
* There's now parent classes Sentinel2DataCollection, Sentinel2Item and Sentinel2Metadata which have datasource specific child classes:
    - GEESentinel2DataCollection
    - AWSSentinel2DataCollection, AWSSentinel2Item, AWSSentinel2Metadata
* The parent classes have methods that are common for both data sources and the child classes have methods that are specific to the data source.
* Improved handling of Sentinel-2 bands and scene classification classes with S2Band and SCLClass classes
* Biophysical processor is now a class SNAPBiophysProcessor, also the biophysical variables and vegetation indices are now Enum classes.
* Enabled easier importing and access of classes and submodules. For example, you can define the data source with satellitetools.DataSource.GEE and the bands with satellitetools.S2Band.B4.
* Added tests.
* Improved docstrings, added examples and updated README.


## Installation ##

```console
pip install git+https://github.com/ollinevalainen/satellitetools.git
```

## About Google Earth Engine ##
Earth Engine must bee authenticated and initialized by the user.
Guides to authentication: https://developers.google.com/earth-engine/guides/auth

Remember GEE's licence terms.

## About AWS cloud-optimized geotiffs ##
Currently (2024-10-23) 2022 data is missing still from the sentinel-2-c1-l2a collection
because ESA has not yet processed them. 2022 is instead fetched from the sentinel-2-l2a
collection. Follow discussion and issues related to this at: https://github.com/Element84/earth-search/issues 

## Warning ##
The biophysical processor implementation in this package does not currently use the 
convex hull check (see the Java code and ATBD) and does not have as extensive 
input/output validity flagging as the original version in SNAP.


## Usage examples ##


**Define area of interest and request parameters:**
```python
from shapely.geometry import Polygon
import satellitetools as sattools
from satellitetools.common.sentinel2 import S2_BANDS_10_20_GEE, S2_FILTER1


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
aoi = sattools.AOI("qvidja_ec", qvidja_polygon, "EPSG:4326")
datestart = "2023-06-01"
dateend = "2023-06-30"
bands = sattools.S2Band.get_10m_to_20m_bands()

```

**GEE:**

You need to initialize (and possibly authenticate) Google Earth Engine and define GEE 
data source.
```python
import ee
ee.Initialize(project="your-ee-project-name")
datasource = sattools.DataSource.GEE
```


**AWS:**

Define AWS data source.
```python
datasource = sattools.DataSource.AWS
```

**Evaluate quality and get the data using a wrapper function:**
```python
req_params = sattools.Sentinel2RequestParams(
    datestart, dateend, datasource, bands
)

df_quality_information, ds_s2_data = sattools.wrappers.get_s2_qi_and_data(
    aoi=aoi,
    req_params=req_params,
    qi_threshold=0.02,
    qi_filter=S2_FILTER1,
)
# df_quality_information will include quality information for all available Sentinel-2 images.

# ds_s2_data will include Sentinel-2 data for the requested bands and metadata for dates passing quality filter (see below more on the quality evaluation).
```
**Alternatively without the wrapper:**

This way you can access the data collection class with some extra metadata and with AWS
the source data item. With AWS data source you can speed-up process by concurrently fetching 
multiple images with multiprocessing parameter.
```python
if req_params.datasource == DataSource.GEE:
    s2_data_collection = GEESentinel2DataCollection(aoi, req_params)
if req_params.datasource == DataSource.AWS:
    import os
    number_of_processes = os.cpu_count() - 1 # or something else
    s2_data_collection = AWSSentinel2DataCollection(aoi, req_params, multiprocessing=number_of_processes)
    s2_data_collection.search_s2_items()

s2_data_collection.get_quality_info()
qi_threshold=0.02 # without specifying default 0 will be used
qi_filter=S2_FILTER1
s2_data_collection.filter_s2_items(qi_threshold, qi_filter)
s2_data_collection.get_s2_data()
s2_data_collection.data_to_xarray()

df_quality_information = s2_data_collection.quality_information
ds_s2_data = s2_data_collection.xr_dataset
```

**Calculate vegetation variables:**
```python
# Example to calculate LAI and NDVI.
variable = sattools.biophys.BiophysVariable.LAI
ds_s2_data = sattools.biophys.run_snap_biophys(ds_s2_data, variable)
vi = sattools.biophys.VegetationIndex.NDVI
ds_s2_data = sattools.biophys.compute_vegetation_index(ds_s2_data, vi)
# Both LAI and NDVI are now included in the dataset with variables names lai and ndvi.

# Other biophysical variables and vegetation indices also available:
class VegetationIndex(str, Enum):
    NDVI = "ndvi"
    CI_RED_EDGE = "ci_red_edge"
    GCC = "gcc"


class BiophysVariable(str, Enum):
    FAPAR = "fapar"
    FCOVER = "fcover"
    LAI = "lai"
    LAI_Cab = "lai_cab"
    LAI_Cw = "lai_cw"
```

## Testing
Tests have been written for pytest.

If you need to specify Google Earth Engine project for ee.Initialize(), for test you 
need to do it by having "EE_PROJECT_PYTEST" environment variable set. For example:
```console
export EE_PROJECT_PYTEST=your-ee-project
```

## About the satellite data filtering process ##

The satellite data is filtered based on the information available in the scene 
classification band (SCL) of the Sentinel-2 Level 2A products. 
In the SCL band every pixel is classified into one of the following classes: 

```python
class SCLClass(Enum):
    """Sentinel-2 Scene Classification Layer (SCL) classes."""
    NODATA = 0
    SATURATED_DEFECTIVE = 1
    DARK_FEATURE_SHADOW = 2
    CLOUD_SHADOW = 3
    VEGETATION = 4
    NOT_VEGETATED = 5
    WATER = 6
    UNCLASSIFIED = 7
    CLOUD_MEDIUM_PROBA = 8
    CLOUD_HIGH_PROBA = 9
    THIN_CIRRUS = 10
    SNOW_ICE = 11
```

The percentage of each of those classes within the polygon/area-of-interest(AOI) is 
calculated from SCL data. This information is/can be used to to 
filter the dates for which the actual data is requested from GEE or AWS. 
There's two parameters, qi_threhold and qi_filter, which are used to filter the data. 
The qi_filter is a list of the (unwanted) classes whose percentages within the AOI is 
summed. Currently, qi_filter is by default this:

```python
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
```
These classes are considered unwanted and user can change the classes included in that
list by the qi_filter parameter. If the sum of the percentages of the classes defined 
in the qi_filter is less than the qi_threshold, the data is considered good/acceptable. 
I have been using the 0.02 which I have empirically tested and found to be ok for our  
purposes.

The SCL band provided in the S2 level-2A products is automatically generated and it has 
errors/misclassifications in it. This affects also the accuracy of the filtering 
processs and some bad data acquisition dates might pass the filter.


## TODO ##
* Documentation with sphinx and push to read the docs
* Examples to docstrings
* Point-based AOI (currently AOI must be a polygon)
* ...
