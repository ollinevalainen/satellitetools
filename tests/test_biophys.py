from pathlib import Path

import xarray as xr

import satellitetools as sattools
from satellitetools.biophys.biophys import BiophysVariable


class TestSNAPBiophysicalProcessor:
    def test_lai_with_aws_data(self):
        variable = BiophysVariable.LAI
        ds = xr.open_dataset(Path(__file__).parent / "test_data" / "test_aws_20m.nc")
        ds = sattools.biophys.run_snap_biophys(ds, variable)
        assert ds.lai is not None
        return ds

    def test_lai_with_gee_data(self):
        variable = BiophysVariable.LAI
        ds = xr.open_dataset(Path(__file__).parent / "test_data" / "test_gee_20m.nc")
        ds = sattools.biophys.run_snap_biophys(ds, variable)
        assert ds.lai is not None
        return ds
