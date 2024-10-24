from pathlib import Path

import xarray as xr

import satellitetools as sattools
from satellitetools.biophys.biophys import BiophysVariable, VegetationIndex


class TestSNAPBiophysicalProcessor:
    def test_lai_with_aws_data(self):
        variable = BiophysVariable.LAI
        ds = xr.open_dataset(Path(__file__).parent / "test_data" / "test_aws_20m.nc")
        ds = sattools.biophys.run_snap_biophys(ds, variable)
        assert ds.lai is not None

    def test_lai_with_gee_data(self):
        variable = BiophysVariable.LAI
        ds = xr.open_dataset(Path(__file__).parent / "test_data" / "test_gee_20m.nc")
        ds = sattools.biophys.run_snap_biophys(ds, variable)
        assert ds.lai is not None


class TestVegetationIndices:
    def test_ndvi(self):
        ds = xr.open_dataset(Path(__file__).parent / "test_data" / "test_aws_20m.nc")
        vi = VegetationIndex.NDVI
        ds = sattools.biophys.compute_vegetation_index(ds, vi)
        assert ds.ndvi is not None

    def test_gcc(self):
        ds = xr.open_dataset(Path(__file__).parent / "test_data" / "test_aws_20m.nc")
        vi = VegetationIndex.GCC
        ds = sattools.biophys.compute_vegetation_index(ds, vi)
        assert ds.gcc is not None

    def test_ci_red_edge(self):
        ds = xr.open_dataset(Path(__file__).parent / "test_data" / "test_aws_20m.nc")
        vi = VegetationIndex.CI_RED_EDGE
        ds = sattools.biophys.compute_vegetation_index(ds, vi)
        assert ds.ci_red_edge is not None
