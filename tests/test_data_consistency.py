from pathlib import Path

import xarray as xr


class TestDataConsistency:
    def test_gee_aws_10m(self):
        ds_gee_10m = xr.open_dataset(
            Path(__file__).parent / "test_data/test_gee_10m.nc"
        )
        ds_aws_10m = xr.open_dataset(
            Path(__file__).parent / "test_data/test_aws_10m.nc"
        )
        assert ds_gee_10m.dims == ds_aws_10m.dims
        assert ds_gee_10m.attrs.keys() == ds_aws_10m.attrs.keys()
        # Check that data is roughly the same
        assert all(
            ds_aws_10m.band_data.mean(dim=["x", "y", "time"]).round(2)
            == ds_gee_10m.band_data.mean(dim=["x", "y", "time"]).round(2)
        )

    def test_gee_aws_20m(self):

        ds_gee_20m = xr.open_dataset(
            Path(__file__).parent / "test_data/test_gee_20m.nc"
        )

        ds_aws_20m = xr.open_dataset(
            Path(__file__).parent / "test_data/test_aws_20m.nc"
        )

        assert ds_gee_20m.dims == ds_aws_20m.dims
        assert ds_gee_20m.attrs.keys() == ds_aws_20m.attrs.keys()
        # Check that data is roughly the same
        assert all(
            ds_aws_20m.band_data.mean(dim=["x", "y", "time"]).round(2)
            == ds_gee_20m.band_data.mean(dim=["x", "y", "time"]).round(2)
        )

    def test_gee_aws_2022(self):

        ds_gee = xr.open_dataset(Path(__file__).parent / "test_data/test_gee_2022.nc")

        ds_aws = xr.open_dataset(Path(__file__).parent / "test_data/test_aws_2022.nc")

        assert ds_gee.dims == ds_aws.dims
        assert ds_gee.attrs.keys() == ds_aws.attrs.keys()
        # Check that data is roughly the same
        assert all(
            ds_aws.band_data.mean(dim=["x", "y", "time"]).round(2)
            == ds_gee.band_data.mean(dim=["x", "y", "time"]).round(2)
        )
