# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:34:08 2020

@author: Olli Nevalainen (olli.nevalainen@fmi.fi),
 Finnish Meteorological Institute)

Olli's python implementation of ESA SNAP s2toolbox biophysical processor and
computation of vegetation indices.
See ATBD at https://step.esa.int/docs/extra/ATBD_S2ToolBox_L2B_V1.1.pdf
And java source code at
https://github.com/senbox-org/s2tbx/tree/master/s2tbx-biophysical/src/main/java/org/esa/s2tbx/biophysical 
SNAP auxdata from:
https://github.com/senbox-org/s2tbx/tree/master/s2tbx-biophysical/src/main/resources/auxdata/2_1

Caveats
Currently changes out of bounds inputs and outputs to nan (or min or max value
if output wihtin tolerance). Maybe output flagging information as well ( i.e.
diffferent flags input and output out of bounds).

Convex hull input checking currently disabled. It's computationally slow and
 not sure of its benefits. Better to filter out bad data based on L2A quality
 info/classification and hope averaging removes some bad pixels. 
"""  # noqa
import os
from enum import Enum

import numpy as np
import xarray as xr

from satellitetools.common.sentinel2 import S2Band


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

    @classmethod
    def get_by_name(cls, name: str) -> "BiophysVariable":
        """Get BiophysVariable by name.

        Parameters
        ----------
        name : str
            Name of the biophysical variable.

        Returns
        -------
        Union[None, "BiophysVariable"]
            BiophysVariable or raise Valuer error if not found.

        """
        for var in BiophysVariable:
            if var.name == name:
                return var
        raise ValueError(f"BiophysVariable {name} not found.")


SNAP_BIO_BANDS = [
    S2Band.B3.value,
    S2Band.B4.value,
    S2Band.B5.value,
    S2Band.B6.value,
    S2Band.B7.value,
    S2Band.B8A.value,
    S2Band.B11.value,
    S2Band.B12.value,
]
SNAP_BIO_RMSE = {
    BiophysVariable.FAPAR: 0.05,
    BiophysVariable.FCOVER: 0.04,
    BiophysVariable.LAI: 0.89,
    BiophysVariable.LAI_Cab: 56,
    BiophysVariable.LAI_Cw: 0.03,
}


class SNAPBiophysProcessor:
    """SNAP Biophysical Processor.

    Attributes
    ----------
    data_cube : xr.DataArray
        Data cube with Sentinel-2 10-20m bands and observation geometry layers
        biophysical bands.
    variable : BiophysVariable
        Biophysical variable to compute.
    """

    def __init__(self, data_cube: xr.DataArray, variable: BiophysVariable):
        """Initialize SNAPBiophysProcessor.

        Parameters
        ----------
        data_cube : xr.DataArray
            Data cube with Sentinel-2 10-20m bands and observation geometry layers
            biophysical bands.
        variable : BiophysVariable
            Biophysical variable to compute.
        """
        self.data_cube = data_cube
        self.variable = variable
        self.nn_params = self._read_snap_nn_params(variable)

    def _read_snap_nn_params(self, variable: BiophysVariable) -> dict:
        # path_to_s2tbx_biophysical
        snap_bio_path = os.path.join(
            os.path.dirname(__file__), "snap-auxdata/biophysical/2_1/"
        )
        var_name_tuple = (variable.name, variable.name)
        norm_minmax = np.loadtxt(
            snap_bio_path + "%s/%s_Normalisation" % var_name_tuple, delimiter=","
        )
        denorm_minmax = np.loadtxt(
            snap_bio_path + "%s/%s_Denormalisation" % var_name_tuple, delimiter=","
        )
        layer1_weights = np.loadtxt(
            snap_bio_path + "%s/%s_Weights_Layer1_Neurons" % var_name_tuple,
            delimiter=",",
        )
        layer1_bias = np.loadtxt(
            snap_bio_path + "%s/%s_Weights_Layer1_Bias" % var_name_tuple,
            delimiter=",",
        ).reshape(-1, 1)
        layer2_weights = np.loadtxt(
            snap_bio_path + "%s/%s_Weights_Layer2_Neurons" % var_name_tuple,
            delimiter=",",
        ).reshape(1, -1)
        layer2_bias = np.loadtxt(
            snap_bio_path + "%s/%s_Weights_Layer2_Bias" % var_name_tuple,
            delimiter=",",
        ).reshape(1, -1)
        extreme_cases = np.loadtxt(
            snap_bio_path + "%s/%s_ExtremeCases" % var_name_tuple, delimiter=","
        )

        defdom_min = np.loadtxt(
            snap_bio_path + "%s/%s_DefinitionDomain_MinMax" % var_name_tuple,
            delimiter=",",
        )[0, :].reshape(-1, 1)
        defdom_max = np.loadtxt(
            snap_bio_path + "%s/%s_DefinitionDomain_MinMax" % var_name_tuple,
            delimiter=",",
        )[1, :].reshape(-1, 1)
        defdom_grid = np.loadtxt(
            snap_bio_path + "%s/%s_DefinitionDomain_Grid" % var_name_tuple,
            delimiter=",",
        )
        nn_params = {
            "norm_minmax": norm_minmax,
            "denorm_minmax": denorm_minmax,
            "layer1_weights": layer1_weights,
            "layer1_bias": layer1_bias,
            "layer2_weights": layer2_weights,
            "layer2_bias": layer2_bias,
            "defdom_min": defdom_min,
            "defdom_max": defdom_max,
            "defdom_grid": defdom_grid,
            "extreme_cases": extreme_cases,
        }
        return nn_params

    def _normalization(self, x):
        x_min = self.nn_params["norm_minmax"][:, 0].reshape(-1, 1)
        x_max = self.nn_params["norm_minmax"][:, 1].reshape(-1, 1)
        x_norm = 2 * (x - x_min) / (x_max - x_min) - 1
        return x_norm

    def _denormalization(self, x):
        x_min = self.nn_params["denorm_minmax"][0]
        x_max = self.nn_params["denorm_minmax"][1]
        x_denorm = 0.5 * (x + 1) * (x_max - x_min)
        return x_denorm

    def _input_ouf_of_range(self, x):
        x_copy = x.copy()
        x_bands = x_copy[:8, :]

        # check min max domain
        defdom_min = self.nn_params["defdom_min"][:, 0].reshape(-1, 1)
        defdom_max = self.nn_params["defdom_max"][:, 0].reshape(-1, 1)
        bad_input_mask = (x_bands < defdom_min) | (x_bands > defdom_max)
        bad_vector = np.any(bad_input_mask, axis=0)
        x_bands[:, bad_vector] = np.nan

        # convex hull check, currently disabled due to time consumption vs benefit
        # gridProject = lambda v: np.floor(
        #     10 * (v - defdom_min) / (defdom_max - defdom_min) + 1
        # ).astype(int)
        # x_bands = gridProject(x_bands)
        # isInGrid = lambda v: any((v == x).all() for x in nn_params[variable]["defdom_grid"]) # noqa
        # notInGrid = ~np.array([isInGrid(v) for v in x_bands.T])
        # x[:, notInGrid | bad_vector] = np.nan

        x_copy[:, bad_vector] = np.nan
        return x_copy

    def _output_ouf_of_range(self, x):
        y = np.copy(x)
        tolerance = self.nn_params["extreme_cases"][0]
        output_min = self.nn_params["extreme_cases"][1]
        output_max = self.nn_params["extreme_cases"][2]

        y[x < (output_min + tolerance)] = np.nan
        y[(x > (output_min + tolerance)) & (x < output_min)] = output_min
        y[(x < (output_max - tolerance)) & (x > output_max)] = output_max
        y[x > (output_max - tolerance)] = np.nan
        return y

    def _layer1(self, x):
        return np.tanh(
            self.nn_params["layer1_weights"].dot(x) + self.nn_params["layer1_bias"]
        )

    def _layer2(self, x):
        return self.nn_params["layer2_weights"].dot(x) + self.nn_params["layer2_bias"]

    def _run_neural_network(self, x):
        y = self._input_ouf_of_range(x)
        y = self._normalization(y)
        y = self._layer1(y)
        y = self._layer2(y)
        y = self._denormalization(y)
        y = self._output_ouf_of_range(y)
        y = y.reshape(1, np.shape(x)[1])
        return y

    def run(self) -> xr.DataArray:
        """Run SNAP Biophysical Processor.

        Returns
        -------
        xr.DataArray
            Array of the computed biophysical variable.

        """

        # Flatten spatial dimensions to 1D for computation
        data_stacked = self.data_cube.stack(yx=("y", "x"))
        # Compute variable and unstack to original 2D shape
        variable_arr = xr.apply_ufunc(
            self._run_neural_network,
            data_stacked,
            input_core_dims=[["band", "yx"]],
            output_core_dims=[["yx"]],
            vectorize=True,
        ).unstack()
        return variable_arr


def create_data_cube(ds: xr.Dataset) -> xr.DataArray:
    """Create data cube for SNAP biophysical processor.

    Parameters
    ----------
    ds : xr.Dataset
        xarray dataset.

    Returns
    -------
    xr.DataArray
        Data cube with Sentinel-2 10-20m bands and observation geometry layers
        biophysical bands.

    """
    band_data = ds.band_data.sel(band=SNAP_BIO_BANDS)
    # generate view angle bands/layers
    vz = (
        np.ones_like(band_data[:, 0, :, :]).T
        * np.cos(np.radians(ds.view_zenith)).values
    )
    vz = vz[..., np.newaxis]
    vzarr = xr.DataArray(
        vz,
        coords=[ds.x, ds.y, ds.time, ["view_zenith"]],
        dims=["x", "y", "time", "band"],
    )

    sz = (
        np.ones_like(band_data[:, 0, :, :]).T * np.cos(np.radians(ds.sun_zenith)).values
    )
    sz = sz[..., np.newaxis]
    szarr = xr.DataArray(
        sz,
        coords=[ds.x, ds.y, ds.time, ["sun_zenith"]],
        dims=["x", "y", "time", "band"],
    )

    raz = (
        np.ones_like(band_data[:, 0, :, :]).T
        * np.cos(np.radians(ds.sun_azimuth - ds.view_azimuth)).values
    )
    raz = raz[..., np.newaxis]
    razarr = xr.DataArray(
        raz,
        coords=[ds.x, ds.y, ds.time, ["relative_azimuth"]],
        dims=["x", "y", "time", "band"],
    )
    data_cube = xr.concat([band_data, vzarr, szarr, razarr], dim="band")
    return data_cube


def run_snap_biophys(ds: xr.Dataset, variable: BiophysVariable) -> xr.Dataset:
    """Compute specified variable using the SNAP algorithm.

    See ATBD at https://step.esa.int/docs/extra/ATBD_S2ToolBox_L2B_V1.1.pdf

    Parameters
    ----------
    ds : xr dataset
        xarray dataset.
    variable : Union[str, BiophysVariable]
        Options 'FAPAR', 'FCOVER', 'LAI', 'LAI_Cab' or 'LAI_Cw'

    Returns
    -------
    xarray dataset
        Adds the specified variable array to the input dataset.

    """
    # If variable is string, convert to BiophysVariable
    if variable in [v.name for v in BiophysVariable]:
        variable = BiophysVariable.get_by_name(variable)

    data_cube = create_data_cube(ds)
    snap_bio_processor = SNAPBiophysProcessor(data_cube, variable)
    output_array = snap_bio_processor.run()
    return ds.assign({variable.value: output_array})


# =============================================================================
# def multidim_intersect(arr1, arr2):
#     arr1_view = arr1.view([('',arr1.dtype)]*arr1.shape[1])
#     arr2_view = arr2.view([('',arr2.dtype)]*arr2.shape[1])
#     isin = np.isin(arr1_view, arr2_view)
#     return isin
# =============================================================================


def compute_ndvi(ds: xr.Dataset) -> xr.Dataset:
    """Compute NDVI

    Parameters
    ----------
    ds : xarray dataset

    Returns
    -------
    xarray dataset
        Adds 'ndvi' xr array to xr dataset.

    """
    b4 = ds.band_data.sel(band=S2Band.B4.value)
    b8 = ds.band_data.sel(band=S2Band.B8A.value)
    ndvi = (b8 - b4) / (b8 + b4)
    return ds.assign({VegetationIndex.NDVI.value: ndvi})


def compute_ci_red_edge(ds: xr.Dataset) -> xr.Dataset:
    """Compute CI_Red_Edge vegetation index.

    Parameters
    ----------
    ds : xarray dataset

    Returns
    -------
    xarray dataset
        Adds 'ci_red_edge' xr array to xr dataset.

    """
    b5 = ds.band_data.sel(band=S2Band.B5.value)
    b7 = ds.band_data.sel(band=S2Band.B7.value)
    ci_red_edge = (b7 / b5) - 1
    return ds.assign({VegetationIndex.CI_RED_EDGE.value: ci_red_edge})


def compute_gcc(ds: xr.Dataset) -> xr.Dataset:
    """Compute GCC vegetation index.

    Parameters
    ----------
    ds : xarray dataset

    Returns
    -------
    xarray dataset
        Adds 'gcc' xr array to xr dataset.

    """
    b2 = ds.band_data.sel(band=S2Band.B2.value)
    b3 = ds.band_data.sel(band=S2Band.B3.value)
    b4 = ds.band_data.sel(band=S2Band.B4.value)
    gcc = b3 / (b2 + b3 + b4)
    return ds.assign({VegetationIndex.GCC.value: gcc})


def compute_vegetation_index(ds: xr.Dataset, vi: VegetationIndex) -> xr.Dataset:
    """Compute vegetation index.

    Parameters
    ----------
    ds : xarray dataset
    vi : VegetationIndex
        Vegetation index to compute.

    Returns
    -------
    xarray dataset
        Adds vegetation index array to xr dataset.

    """
    if vi == VegetationIndex.NDVI:
        return compute_ndvi(ds)
    if vi == VegetationIndex.CI_RED_EDGE:
        return compute_ci_red_edge(ds)
    if vi == VegetationIndex.GCC:
        return compute_gcc(ds)
    raise ValueError(f"Vegetation index {vi} not found.")
