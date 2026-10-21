
from pathlib import Path
import pickle

import numpy as np
from netCDF4 import Dataset


# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

LAYERS = {
    "top": 0,
    "bottom": -1,
}

REGIONS = {
    "gp": {
        "name": "Great Plains",
        "latitude": (37.0, 44.0),
        "longitude": (360 - 104.0, 360 - 97.0),
    },
    "ca": {
        "name": "California",
        "latitude": (32.5, 42.0),
        "longitude": (360 - 124.5, 360 - 114.1),
    },
    "fl": {
        "name": "Florida",
        "latitude": (24.4, 31.0),
        "longitude": (360 - 87.6, 360 - 80.0),
    },
}

RUSSIA_BOUNDS = {
    "latitude": (70.0, 90.0),
    "longitude": (90.0, 120.0),
}

NUMBER_OF_ENSEMBLES = 500

INITIAL_CONDITIONS_FILE = (
    "/fs/ess/PAS2635/Generalized_Predictability_MPAS/"
    "120km_uniform/history.2025-10-01_06.00.00.nc"
)

ENSEMBLE_DIRECTORY = Path(
    "/fs/ess/PAS2635/LandAir_Predictability"
)

OUTPUT_FILE = Path(
    "/fs/scratch/PAS3252/yang/HONORS_THESIS/xs_dict.pkl"
)


# ---------------------------------------------------------------------
# Coordinates and geographic masks
# ---------------------------------------------------------------------

with Dataset(INITIAL_CONDITIONS_FILE) as nc:
    latitudes = np.rad2deg(np.asarray(nc["latCell"]))
    longitudes = np.rad2deg(np.asarray(nc["lonCell"]))


def create_geographic_mask(latitude_bounds, longitude_bounds):
    """Return a Boolean mask for a latitude-longitude bounding box."""
    lat_min, lat_max = latitude_bounds
    lon_min, lon_max = longitude_bounds

    return (
        (latitudes >= lat_min)
        & (latitudes <= lat_max)
        & (longitudes >= lon_min)
        & (longitudes <= lon_max)
    )


region_masks = {
    region_key: create_geographic_mask(
        region["latitude"],
        region["longitude"],
    )
    for region_key, region in REGIONS.items()
}

russia_mask = create_geographic_mask(
    RUSSIA_BOUNDS["latitude"],
    RUSSIA_BOUNDS["longitude"],
)


# ---------------------------------------------------------------------
# Ensemble-mean soil moisture and Russian geopotential height
# ---------------------------------------------------------------------

def create_smois_arrays(layer_index, region_mask):
    """
    Calculate ensemble-mean soil moisture for a region and layer.

    Also calculate ensemble-mean 500-hPa height over Russia.

    Returns
    -------
    regional_soil_moisture : ndarray
        One regional mean per ensemble member.
    russian_height : ndarray
        One Russian mean 500-hPa height per ensemble member.
    """
    regional_soil_moisture = np.empty(NUMBER_OF_ENSEMBLES)
    russian_height = np.empty(NUMBER_OF_ENSEMBLES)

    for ensemble_index in range(NUMBER_OF_ENSEMBLES):
        member_number = ensemble_index + 1

        input_file = (
            ENSEMBLE_DIRECTORY
            / f"member_{member_number:05d}"
            / "diag.2021-07-14_00.00.00.nc"
        )

        with Dataset(input_file) as nc:
            soil_moisture = np.squeeze(
                np.asarray(nc["smois"])
            )[:, layer_index]

            height_500hpa = np.squeeze(
                np.asarray(nc["height_500hPa"])
            )

        # Exclude ocean cells, where soil moisture equals 1.
        regional_values = soil_moisture[region_mask]
        regional_values = regional_values[regional_values != 1]

        regional_soil_moisture[ensemble_index] = np.mean(
            regional_values
        )

        russian_height[ensemble_index] = np.mean(
            height_500hpa[russia_mask]
        )

    return regional_soil_moisture, russian_height


# ---------------------------------------------------------------------
# Build output dictionary
# ---------------------------------------------------------------------

xs_dict = {}

for region_key, region in REGIONS.items():
    print(f"Processing {region['name']}...")

    for layer_name, layer_index in LAYERS.items():
        print(f"  Soil layer: {layer_name}")

        soil_moisture, russian_height = create_smois_arrays(
            layer_index,
            region_masks[region_key],
        )

        output_key = f"x_{'t' if layer_name == 'top' else 'b'}_{region_key}"
        xs_dict[output_key] = soil_moisture

        # Preserve the original x2 key and its ensemble ordering.
        if "x2" not in xs_dict:
            xs_dict["x2"] = russian_height


# ---------------------------------------------------------------------
# Save results
# ---------------------------------------------------------------------

OUTPUT_FILE.parent.mkdir(parents=True, exist_ok=True)

with OUTPUT_FILE.open("wb") as output_file:
    pickle.dump(xs_dict, output_file)

print(f"Saved results to {OUTPUT_FILE}")