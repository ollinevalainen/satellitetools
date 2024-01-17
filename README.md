[![DOI](https://zenodo.org/badge/270676132.svg)](https://zenodo.org/badge/latestdoi/270676132)

# satellitetools
This package is at the moment under constant development and use in our research projects at FMI (Finnish Meteorological Institute). Currently developed and used in agricultural study areas. Sorry for lack of documentation at the moment. Please feel free to contact me in case you want more information about this package.

This package has currently tools for getting Sentinel-2 data from Google Earth Engine (GEE) or alternatively from AWS Open Data registry where Setinel-2 data (Level-2A) data is available as cloud-optimized geotiffs (https://registry.opendata.aws/sentinel-2-l2a-cogs/). There's also my python implementation of ESA's SNAP Biophysical processor which can be used to compute biophysical parameters, such as LAI (Original Java code in SNAP here: https://github.com/senbox-org/s2tbx/tree/master/s2tbx-biophysical/src/main/java/org/esa/s2tbx/biophysical). ATBD: http://step.esa.int/docs/extra/ATBD_S2ToolBox_L2B_V1.1.pdf.

See examples/gee_example.py and examples/aws_cog examples for simple usage.

**WARNING**:
* GEE data is currently retrieved with 10m resolution (scale=10), so the 20m resolution bands are resampled.
* The biophysical processor implementation in this package does not currently use the convex hull check (see the Java code and ATBD) and does not have as extensive input/output validity flagging as the original version in SNAP.

**TODO**:
* Proper testing
* Point-based AOI (currently AOI must be a polygon)
* ...

**Installation**

```console
pip install git+https://github.com/ollinevalainen/satellitetools.git@develop
```


**Some explanation of the satellite data filtering process used:**

The satellite data is filtered based on the information available in the scene classification band (SCL) of the Sentinel-2 Level 2A products. 
In the SCL band every pixel is classified into one of the following classes: 

[
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

The percentage of each of those classes within the polygon/area-of-interest(AOI) is calculated from SCL data (get_quality_info). This information is/can be used to to filter the dates for which the actual data is requested from GEE or AWS. There's two parameters, qi_threhold and qi_filter, which are used to filter the data. The qi_filter is a list of the (unwanted) classes whose percentages within the AOI is summed. Currently, qi_filter is by default this:

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

These classes are considered unwanted and user can change the classes included in that list by the qi_filter parameter. If the sum of the percentages of the classes defined in the qi_filter is less than the qi_threshold, the data is considered good/acceptable. I have been using the 0.02 which I have empirically tested and found to be ok for our purposes.

The SCL band provided in the S2 level-2A products is automatically generated and it has errors/misclassifications in it. This affects also the accuracy of the filtering processs and some bad data acquisition dates usually passes the filter.
