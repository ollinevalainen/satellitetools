from pathlib import Path

import xarray as xr

import satellitetools as sattools
from satellitetools.biophys.biophys import BiophysVariable, VegetationIndex


class TestTimeseries:
    def test_timeseries(self):
        variable = BiophysVariable.LAI
        ds = xr.open_dataset(Path(__file__).parent / "test_data" / "test_aws_20m.nc")
        ds = sattools.biophys.run_snap_biophys(ds, variable)
        ds = sattools.biophys.compute_ndvi(ds)
        df = sattools.timeseries.xr_dataset_to_timeseries(
            ds,
            [BiophysVariable.LAI, VegetationIndex.NDVI],
            add_uncertainty=True,
            add_confidence_intervals=True,
            confidence_level="95",
        )
        assert not df.empty
