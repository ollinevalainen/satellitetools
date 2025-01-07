# Changelog

## satellitetools 2.1.0

### Feature
* Set CI-CD workflow using Github Actions and Python Semantic Release

### Fix
* Fixed printing/logging bugs

## satellitetools 2.0.0 🛰️ 
This is the second release of satellitetools package :tada:

The package ended-up being quite useful in multiple research projects, so development of the package has been continuous. This release is a milestone for years of small changes and recent major refactoring of code during past two months. But now develop branch was finally merged to master.

In the future, the plan is to publish releases more frequently when possible developments in the develop branch gets merged to master branch.

Too much has changed compared to the first release, but here's summary of changes respect to the previous develop branch before major refactoring done during Sep-Dec 2024:

 **Breaking changes:**:warning:

* Renamed xrtools.py to timeseries.py
* RequestParams class is now Sentinel2RequestParams in sentinel2 submodule
* Restructured files and removed nesting. For example gee imported as satellitetools.gee instead of satellitetools.gee.gee
* The quality information dataframe doesn't have anymore `"Date"` column, but instead index named `"acquisition_time"` as UTC aware timestamp (from `pd.to_datetime`).

**New features:** 🔧 

* Refactored and made codes more object-oriented and modular:
* There's now parent classes `Sentinel2DataCollection`, `Sentinel2Item` and `Sentinel2Metadata` which have datasource specific child classes:
    - `GEESentinel2DataCollection`
    - `AWSSentinel2DataCollection`, `AWSSentinel2Item`, `AWSSentinel2Metadata`
* The parent classes have methods that are common for both data sources and the child classes have methods that are specific to the data source.
* Improved handling of Sentinel-2 bands and scene classification classes with `S2Band` and `SCLClass` classes
* Biophysical processor is now a class `SNAPBiophysProcessor`, also the biophysical variables and vegetation indices are now Enum classes.
* Enabled easier importing and access of classes and submodules. For example, you can define the data source with `satellitetools.DataSource.GEE` and the bands with `satellitetools.S2Band.B4`.
* Added tests for pytest.
* Improved docstrings, added examples and updated README.
* Documentation using sphinx
* Testing for python versions 3.10, 3.11 and 3.12 using nox
* Changed from print statements to logging.

## 1.0.0
This release is currently (7 February 2022) in use in Field Observatory (https://www.fieldobservatory.org/) and cited in a scientific article "Towards agricultural soil carbon monitoring, reporting and verification through Field Observatory Network (FiON)" that describes Field Observatory and has been accepted for publishing in EGU's Geoscientific Instrumentation, Methods and Data Systems journal (https://doi.org/10.5194/gi-2021-21).