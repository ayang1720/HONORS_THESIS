
"""
Alex Yang — Honors Thesis

Project directory:
    /users/PAS3252/ayang1720/HONORS_THESIS

Input data:
    /fs/ess/PAS2635/LandAir_Predictability

Output data:
    /fs/scratch/PAS3252/yang/HONORS_THESIS
"""

from datetime import timedelta
from pathlib import Path

import numpy as np
from netCDF4 import Dataset
from scipy.stats import beta, false_discovery_control


# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

NUMBER_OF_ENSEMBLES = 500
NUMBER_OF_CELLS = 40962
HOURS_PER_DAY = 24
OUTPUT_INTERVAL_HOURS = 3
INTERVALS_PER_DAY = HOURS_PER_DAY // OUTPUT_INTERVAL_HOURS

RESPONSE_VARIABLES = ["height_500hPa"]

ENSEMBLE_DIRECTORY = Path(
    "/fs/ess/PAS2635/LandAir_Predictability"
)


# ---------------------------------------------------------------------
# Construct daily-mean response arrays
# ---------------------------------------------------------------------

def construct_array(
    responses,
    start_time,
    number_of_ensembles,
    ncells,
):
    """
    Construct daily-mean response variables for every ensemble member.

    Parameters
    ----------
    responses : list[str]
        NetCDF variable names to extract.
    start_time : datetime
        Starting time of the daily averaging period.
    number_of_ensembles : int
        Number of ensemble members.
    ncells : int
        Number of horizontal grid cells.

    Returns
    -------
    ndarray
        Array with shape:
        (number of responses, number of ensembles, number of cells).

        Each value is the mean of eight 3-hourly outputs over one day.
    """
    day = start_time.strftime("%m/%d")

    daily_means = np.zeros(
        (len(responses), number_of_ensembles, ncells)
    )

    # The eight output times spanning one day.
    output_times = [
        start_time + timedelta(hours=step * OUTPUT_INTERVAL_HOURS)
        for step in range(INTERVALS_PER_DAY)
    ]

    for ensemble_index in range(number_of_ensembles):
        member_directory = (
            ENSEMBLE_DIRECTORY
            / f"member_{ensemble_index + 1:05d}"
        )

        # Accumulate the eight 3-hourly values for each response.
        for output_time in output_times:
            timestamp = output_time.strftime("%Y-%m-%d_%H.%M.%S")

            input_file = member_directory / f"diag.{timestamp}.nc"

            with Dataset(input_file) as nc:
                for response_index, response_name in enumerate(responses):
                    response_values = np.squeeze(
                        np.asarray(nc[response_name])
                    )

                    daily_means[response_index, ensemble_index, :] += (
                        response_values / INTERVALS_PER_DAY
                    )

        print(
            f"Day: {day} | "
            f"Ensemble: {ensemble_index + 1}/{number_of_ensembles}"
        )

    return daily_means


# ---------------------------------------------------------------------
# Correlation across ensemble members
# ---------------------------------------------------------------------

def cor_pvalues(soil_moisture, response_array):
    """
    Calculate correlations across ensembles, independently at each cell.

    Parameters
    ----------
    soil_moisture : ndarray, shape (ensembles,)
        Soil moisture values for the ensemble members.
    response_array : ndarray, shape (cells, ensembles)
        Response values for each cell and ensemble member.

    Returns
    -------
    ndarray, shape (cells,)
        Pearson correlation coefficients.
    """
    number_of_ensembles = len(soil_moisture)

    soil_moisture_perturbations = (
        soil_moisture - np.mean(soil_moisture)
    )

    response_perturbations = (
        response_array
        - np.mean(response_array, axis=1, keepdims=True)
    )

    covariance = np.sum(
        soil_moisture_perturbations[None, :]
        * response_perturbations,
        axis=1,
    ) / (number_of_ensembles - 1)

    soil_moisture_std = np.std(soil_moisture, ddof=1)

    response_std = np.std(
        response_array,
        axis=1,
        ddof=1,
    )

    with np.errstate(divide="ignore", invalid="ignore"):
        correlations = covariance / (
            soil_moisture_std * response_std
        )

    return correlations.astype(float)


# ---------------------------------------------------------------------
# Calculate p-values and apply false-discovery-rate correction
# ---------------------------------------------------------------------

def pvalue(correlations, number_of_ensembles):
    """
    Calculate two-sided correlation p-values and apply the
    Benjamini-Yekutieli false-discovery-rate correction.

    Parameters
    ----------
    correlations : ndarray
        Correlation coefficients, typically shaped
        (responses, cells).
    number_of_ensembles : int
        Number of ensemble members used to estimate each correlation.

    Returns
    -------
    ndarray
        Corrected p-values, with the same shape as correlations.
    """
    degrees_distribution = beta(
        number_of_ensembles / 2 - 1,
        number_of_ensembles / 2 - 1,
        loc=-1,
        scale=2,
    )

    raw_pvalues = 2 * degrees_distribution.cdf(
        -np.abs(correlations)
    )

    return false_discovery_control(
        raw_pvalues,
        axis=0,
        method="by",
    )


# ---------------------------------------------------------------------
# GLACE-ESA correlation analysis
# ---------------------------------------------------------------------

def glace_esa_pvalues(
    glaced_smois_array,
    response_variables,
    time,
    number_of_ensembles,
    ncells,
):
    """
    Calculate corrected p-values relating GLACE soil moisture to
    daily-mean response variables at a specified time.

    Parameters
    ----------
    glaced_smois_array : ndarray, shape (ensembles,)
        GLACE soil moisture values from the initial simulation.
    response_variables : list[str]
        NetCDF response variable names.
    time : datetime
        Starting time of the response averaging period.
    number_of_ensembles : int
        Number of ensemble members.
    ncells : int
        Number of horizontal grid cells.

    Returns
    -------
    ndarray, shape (responses, cells)
        False-discovery-rate-corrected p-values.
    """
    global_arrays = construct_array(
        responses=response_variables,
        start_time=time,
        number_of_ensembles=number_of_ensembles,
        ncells=ncells,
    )

    correlations = np.empty(
        (len(response_variables), ncells)
    )

    for response_index, _ in enumerate(response_variables):
        # Original shape: (ensembles, cells).
        response_data = global_arrays[response_index]

        # Correlation function expects: (cells, ensembles).
        response_data = response_data.T

        correlations[response_index] = cor_pvalues(
            glaced_smois_array,
            response_data,
        )

    return pvalue(
        correlations,
        number_of_ensembles,
    )