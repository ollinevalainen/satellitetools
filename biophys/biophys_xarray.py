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

Caveats
Currently changes out of bounds inputs and outputs to nan (or min or max value
if output wihtin tolerance). Maybe output flagging information as well ( i.e.
diffferent flags input and output out of bounds).

Convex hull input checking currently disabled. It's computationally slow and
 not sure of its benefits. Better to filter out bad data based on L2A quality
 info/classification\
    and hope averaging removes some bad pixels.
"""
import os
import numpy as np
import xarray as xr

# Read SNAP Biophysical processor neural network parameters

SNAP_BIO_BANDS = ["B3", "B4", "B5", "B6", "B7", "B8A", "B11", "B12"]
SNAP_BIO_RMSE = {
    "fapar": 0.05,
    "fcover": 0.04,
    "lai": 0.89,
    "lai_cab": 56,
    "lai_cw": 0.03,
}
# path_to_s2tbx_biophysical
snap_bio_path = os.path.join(os.path.dirname(__file__), "snap-auxdata/biophysical/2_1/")
nn_params = {}
for var in ["FAPAR", "FCOVER", "LAI", "LAI_Cab", "LAI_Cw"]:
    norm_minmax = np.loadtxt(
        snap_bio_path + "%s/%s_Normalisation" % (var, var), delimiter=","
    )
    denorm_minmax = np.loadtxt(
        snap_bio_path + "%s/%s_Denormalisation" % (var, var), delimiter=","
    )
    layer1_weights = np.loadtxt(
        snap_bio_path + "%s/%s_Weights_Layer1_Neurons" % (var, var), delimiter=","
    )
    layer1_bias = np.loadtxt(
        snap_bio_path + "%s/%s_Weights_Layer1_Bias" % (var, var), delimiter=","
    ).reshape(-1, 1)
    layer2_weights = np.loadtxt(
        snap_bio_path + "%s/%s_Weights_Layer2_Neurons" % (var, var), delimiter=","
    ).reshape(1, -1)
    layer2_bias = np.loadtxt(
        snap_bio_path + "%s/%s_Weights_Layer2_Bias" % (var, var), delimiter=","
    ).reshape(1, -1)
    extreme_cases = np.loadtxt(
        snap_bio_path + "%s/%s_ExtremeCases" % (var, var), delimiter=","
    )

    defdom_min = np.loadtxt(
        snap_bio_path + "%s/%s_DefinitionDomain_MinMax" % (var, var), delimiter=","
    )[0, :].reshape(-1, 1)
    defdom_max = np.loadtxt(
        snap_bio_path + "%s/%s_DefinitionDomain_MinMax" % (var, var), delimiter=","
    )[1, :].reshape(-1, 1)
    defdom_grid = np.loadtxt(
        snap_bio_path + "%s/%s_DefinitionDomain_Grid" % (var, var), delimiter=","
    )
    nn_params[var] = {
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


def _normalization(x, x_min, x_max):
    x_norm = 2 * (x - x_min) / (x_max - x_min) - 1
    return x_norm


def _denormalization(y_norm, y_min, y_max):
    y = 0.5 * (y_norm + 1) * (y_max - y_min)
    return y


# LAI tests
# x = np.array([0.057979, 0.0078856, 0.093585, 0.2585, 0.28253, 0.30874, 0.1708, 0.069808, 0.98434, 0.40581, -0.55142]).reshape(-1,1)
# x1 =   np.array([0.056024,0.012462	,0.088543,	0.41626	,0.49575,	0.51452	,0.14425,	0.043583,	0.99367,	0.90957	,-0.99999]).reshape(-1,1)

# =============================================================================
# def multidim_intersect(arr1, arr2):
#     arr1_view = arr1.view([('',arr1.dtype)]*arr1.shape[1])
#     arr2_view = arr2.view([('',arr2.dtype)]*arr2.shape[1])
#     isin = np.isin(arr1_view, arr2_view)
#     return isin
# =============================================================================


def _input_ouf_of_range(x, variable):
    x_copy = x.copy()
    x_bands = x_copy[:8, :]

    # check min max domain
    defdom_min = nn_params[variable]["defdom_min"][:, 0].reshape(-1, 1)
    defdom_max = nn_params[variable]["defdom_max"][:, 0].reshape(-1, 1)
    bad_input_mask = (x_bands < defdom_min) | (x_bands > defdom_max)
    bad_vector = np.any(bad_input_mask, axis=0)
    x_bands[:, bad_vector] = np.nan

    # convex hull check, currently disabled due to time consumption vs benefit
    # gridProject = lambda v: np.floor(10 * (v - defdom_min) / (defdom_max - defdom_min) + 1 ).astype(int)
    # x_bands = gridProject(x_bands)
    # isInGrid = lambda v: any((v == x).all() for x in nn_params[variable]['defdom_grid'])
    # notInGrid = ~np.array([isInGrid(v) for v in x_bands.T])
    # x[:,notInGrid | bad_vector] = np.nan

    x_copy[:, bad_vector] = np.nan
    return x_copy


def _output_ouf_of_range(output, variable):
    new_output = np.copy(output)
    tolerance = nn_params[variable]["extreme_cases"][0]
    output_min = nn_params[variable]["extreme_cases"][1]
    output_max = nn_params[variable]["extreme_cases"][2]

    new_output[output < (output_min + tolerance)] = np.nan
    new_output[(output > (output_min + tolerance)) & (output < output_min)] = output_min
    new_output[(output < (output_max - tolerance)) & (output > output_max)] = output_max
    new_output[output > (output_max - tolerance)] = np.nan
    return new_output


def _compute_variable(x, variable):

    x_norm = np.zeros_like(x)
    x = _input_ouf_of_range(x, variable)
    x_norm = _normalization(
        x,
        nn_params[variable]["norm_minmax"][:, 0].reshape(-1, 1),
        nn_params[variable]["norm_minmax"][:, 1].reshape(-1, 1),
    )

    out_layer1 = np.tanh(
        nn_params[variable]["layer1_weights"].dot(x_norm)
        + nn_params[variable]["layer1_bias"]
    )
    out_layer2 = (
        nn_params[variable]["layer2_weights"].dot(out_layer1)
        + nn_params[variable]["layer2_bias"]
    )
    output = _denormalization(
        out_layer2,
        nn_params[variable]["denorm_minmax"][0],
        nn_params[variable]["denorm_minmax"][1],
    )[0]
    output = _output_ouf_of_range(output, variable)
    output = output.reshape(1, np.shape(x)[1])
    return output


def _s2_lists_to_pixel_vectors(single_date_dict):
    band_list = ["B3", "B4", "B5", "B6", "B7", "B8A", "B11", "B12"]
    pixel_vector = np.zeros(shape=(11, len(single_date_dict[band_list[0]])))

    for i, b in enumerate(band_list):
        pixel_vector[i, :] = np.array(single_date_dict[b]) / 10000.0

    pixel_vector[8, :] = np.cos(np.radians(single_date_dict["view_zenith"]))
    pixel_vector[9, :] = np.cos(np.radians(single_date_dict["sun_zenith"]))
    pixel_vector[10, :] = np.cos(
        np.radians(single_date_dict["sun_azimuth"] - single_date_dict["view_azimuth"])
    )
    return pixel_vector


def compute_ndvi(dataset):
    """Compute NDVI

    Parameters
    ----------
    dataset : xarray dataset

    Returns
    -------
    xarray dataset
        Adds 'ndvi' xr array to xr dataset.

    """
    b4 = dataset.band_data.sel(band="B4")
    b8 = dataset.band_data.sel(band="B8A")
    ndvi = (b8 - b4) / (b8 + b4)
    return dataset.assign({"ndvi": ndvi})


def compute_ci_red_edge(dataset):
    """Compute CI_Red_Edge vegetation index.

    Parameters
    ----------
    dataset : xarray dataset

    Returns
    -------
    xarray dataset
        Adds 'ci_red_edge' xr array to xr dataset.

    """
    b5 = dataset.band_data.sel(band="B5")
    b7 = dataset.band_data.sel(band="B7")
    ci_red_edge = (b7 / b5) - 1
    return dataset.assign({"ci_red_edge": ci_red_edge})


def compute_gcc(dataset):
    """Compute GCC vegetation index.

    Parameters
    ----------
    dataset : xarray dataset

    Returns
    -------
    xarray dataset
        Adds 'gcc' xr array to xr dataset.

    """
    b2 = dataset.band_data.sel(band="B2")
    b3 = dataset.band_data.sel(band="B3")
    b4 = dataset.band_data.sel(band="B4")
    gcc = b3 / (b2 + b3 + b4)
    return dataset.assign({"gcc": gcc})


def run_snap_biophys(dataset, variable):
    """Compute specified variable using the SNAP algorithm.

    See ATBD at https://step.esa.int/docs/extra/ATBD_S2ToolBox_L2B_V1.1.pdf

    Parameters
    ----------
    dataset : xr dataset
        xarray dataset.
    variable : str
        Options 'FAPAR', 'FCOVER', 'LAI', 'LAI_Cab' or 'LAI_Cw'

    Returns
    -------
    xarray dataset
        Adds the specified variable array to dataset (variable name in
        lowercase).

    """
    band_data = dataset.band_data.sel(band=SNAP_BIO_BANDS)
    # generate view angle bands/layers
    vz = (
        np.ones_like(band_data[:, 0, :, :]).T
        * np.cos(np.radians(dataset.view_zenith)).values
    )
    vz = vz[..., np.newaxis]
    vzarr = xr.DataArray(
        vz,
        coords=[dataset.x, dataset.y, dataset.time, ["view_zenith"]],
        dims=["x", "y", "time", "band"],
    )

    sz = (
        np.ones_like(band_data[:, 0, :, :]).T
        * np.cos(np.radians(dataset.sun_zenith)).values
    )
    sz = sz[..., np.newaxis]
    szarr = xr.DataArray(
        sz,
        coords=[dataset.x, dataset.y, dataset.time, ["sun_zenith"]],
        dims=["x", "y", "time", "band"],
    )

    raz = (
        np.ones_like(band_data[:, 0, :, :]).T
        * np.cos(np.radians(dataset.sun_azimuth - dataset.view_azimuth)).values
    )
    raz = raz[..., np.newaxis]
    razarr = xr.DataArray(
        raz,
        coords=[dataset.x, dataset.y, dataset.time, ["relative_azimuth"]],
        dims=["x", "y", "time", "band"],
    )

    newarr = xr.concat([band_data, vzarr, szarr, razarr], dim="band")
    newarr = newarr.stack(yx=("y", "x"))
    arr = xr.apply_ufunc(
        _compute_variable,
        newarr,
        input_core_dims=[["band", "yx"]],
        output_core_dims=[["yx"]],
        kwargs={"variable": variable},
        vectorize=True,
    ).unstack()
    return dataset.assign({variable.lower(): arr})


def estimate_gpp_vi_lue(vi, daily_par, model_name):
    """Estimate GPP using simple vegetation index based models and PAR.

    This function has not been properly tested (i.e.used for a while)

    Parameters
    ----------
    vi : float
        Vegetation index values.
    daily_par : float
        Daily PAR as MJ/s/m².
    model_name : str, optional
        Name of the model (see biophys_xarray.GPP_LUE_models).

    Returns
    -------
    gpp : float
        Estimated gross primary productivity.

    """
    vi_name = "_".join(model_name.split("_")[:-1])
    gpp = GPP_LUE_MODELS[vi_name][model_name]["model"](vi, daily_par)
    return gpp


# GPP estimation models
GPP_LUE_MODELS = {
    "ci_red_edge": {
        "ci_red_edge_1": {
            "model": lambda vi, par: 4.80 * np.log(vi * par * 1e3) - 37.93,
            "species": "soybean",
            "reference": "Peng & Gitelson, 2012",
        },
        "ci_red_edge_2": {
            "model": lambda vi, par: 0.31 * (vi * par) - 0.1,
            "species": "grass",
            "reference": "Huang et al. 2019",
        },
    },
    "ci_green": {
        "ci_green_1": {
            "model": lambda vi, par: 5.13 * np.log(vi * par * 1e3) - 46.92,
            "species": "soybean",
            "reference": "Peng & Gitelson, 2012",
        },
        "ci_green_2": {
            "model": lambda vi, par: 14.7 * np.log(vi * par * 1e3 + 27900.61) - 154,
            "species": "maize",
            "reference": "Peng & Gitelson, 2012",
        },
    },
    "NDVI": {
        "NDVI_1": {
            "model": lambda vi, par: 2.07 * (vi * par) - 6.19,
            "species": "soybean",
            "reference": "Gitelson et al., 2012",
        },
        "NDVI_2": {
            "model": lambda vi, par: 3.11 * (vi * par) - 9.22,
            "species": "maize",
            "reference": "Gitelson et al., 2012",
        },
        "NDVI_3": {
            "model": lambda vi, par: (
                -3.26 * 1e-8 * (vi * par * 1e3) ** 2
                + 1.7 * 1e-3 * (vi * par * 1e3)
                - 2.17
            ),
            "species": "soybean",
            "reference": "Peng & Gitelson, 2012",
        },
        "NDVI_4": {
            "model": lambda vi, par: 1.94e-3 * (vi * par * 1e3) - 2.59,
            "species": "maize",
            "reference": "Peng & Gitelson, 2012",
        },
    },
    "gndvi": {
        "gndvi_1": {
            "model": lambda vi, par: 2.86 * (vi * par) - 11.9,
            "species": "soybean",
            "reference": "Gitelson et al., 2012",
        },
        "gndvi_2": {
            "model": lambda vi, par: 4 * (vi * par) - 15.4,
            "species": "maize",
            "reference": "Gitelson et al., 2012",
        },
    },
    "evi": {
        "evi_1": {
            "model": lambda vi, par: (2.26 * (vi * par) - 3.73),
            "species": "soybean",
            "reference": "Peng et al., 2013",
        },
        "evi_2": {
            "model": lambda vi, par: (3.49 * (vi * par) - 4.92),
            "species": "maize",
            "reference": "Peng et al., 2013",
        },
    },
    "reNDVI": {
        "reNDVI_1": {
            "model": lambda vi, par: 1.61 * (vi * par) - 1.75,
            "species": "mixed",
            "reference": "Wolanin et al., 2019",
        },
        "reNDVI_2": {
            "model": lambda vi, par: (
                -1.19 * 1e-7 * (vi * par * 1e3) ** 2
                + 3 * 1e-3 * (vi * par * 1e3)
                - 2.70
            ),
            "species": "soybean",
            "reference": "Peng & Gitelson, 2012",
        },
        "reNDVI_3": {
            "model": lambda vi, par: (
                -3.41 * 1e-8 * (vi * par * 1e3) ** 2
                + 2.77 * 1e-3 * (vi * par * 1e3)
                - 2.06
            ),
            "species": "maize",
            "reference": "Peng & Gitelson, 2012",
        },
    },
    "fapar": {
        "fapar_1": {
            "model": lambda vi, par: 1.10 * (vi * par),
            "species": "grass",
            "reference": "Olli Qvidja test fapar*x*PAR",
        }
    },
}
