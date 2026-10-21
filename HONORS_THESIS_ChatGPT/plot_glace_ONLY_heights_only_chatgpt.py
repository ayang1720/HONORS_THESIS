
"""
Alex Yang — Honors Thesis
MPAS GLACE-ESA significance time series.

NOTE:
Update the output filename if creating a different type of plot.
"""

import pickle
import sys
from datetime import datetime
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np

import glace_esa_heights_only


# ---------------------------------------------------------------------
# Settings
# ---------------------------------------------------------------------

BASE_DIRECTORY = Path(
    "/fs/scratch/PAS3252/yang/HONORS_THESIS"
)

PLOT_DIRECTORY = (
    BASE_DIRECTORY / "timeseries" / "smois_6_cases"
)

PICKLE_DIRECTORY = BASE_DIRECTORY / "statgood_pickle_files"
SOIL_MOISTURE_FILE = BASE_DIRECTORY / "xs_dict.pkl"
ESA_FILE = Path("ESA.pkl")

NUMBER_OF_ENSEMBLES = 500
NUMBER_OF_CELLS = 40962

RESPONSE_VARIABLE = "height_500hPa"

START_TIME = datetime.strptime(
    "20210714210000", "%Y%m%d%H%M%S"
)

END_TIME = datetime.strptime(
    "20210801210000", "%Y%m%d%H%M%S"
)

# Requested analysis time, supplied as a command-line argument.
analysis_time = datetime.strptime(
    sys.argv[1], "%Y%m%d%H%M%S"
)


# ---------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------

with SOIL_MOISTURE_FILE.open("rb") as file:
    soil_moisture_data = pickle.load(file)

with ESA_FILE.open("rb") as file:
    esa_data = pickle.load(file)

latitudes = esa_data["lat_1d"]
longitudes = esa_data["lon_1d"]


# ---------------------------------------------------------------------
# Define the six cases
# ---------------------------------------------------------------------

# Each case contains its soil-moisture key and subplot position.
CASES = {
    "case_1": {
        "soil_moisture_key": "x_t_gp",
        "subplot": (0, 0),
        "label": "Top soil — Great Plains",
    },
    "case_2": {
        "soil_moisture_key": "x_b_gp",
        "subplot": (0, 1),
        "label": "Bottom soil — Great Plains",
    },
    "case_3": {
        "soil_moisture_key": "x_t_ca",
        "subplot": (1, 0),
        "label": "Top soil — California",
    },
    "case_4": {
        "soil_moisture_key": "x_b_ca",
        "subplot": (1, 1),
        "label": "Bottom soil — California",
    },
    "case_5": {
        "soil_moisture_key": "x_t_fl",
        "subplot": (2, 0),
        "label": "Top soil — Florida",
    },
    "case_6": {
        "soil_moisture_key": "x_b_fl",
        "subplot": (2, 1),
        "label": "Bottom soil — Florida",
    },
}


# ---------------------------------------------------------------------
# Create figure and map axes
# ---------------------------------------------------------------------

fig, axes = plt.subplots(
    nrows=3,
    ncols=2,
    figsize=(24, 15),
    constrained_layout=True,
    subplot_kw={"projection": ccrs.PlateCarree()},
)

for ax in axes.flat:
    ax.set_extent(
        [-180, 180, -90, 90],
        crs=ccrs.PlateCarree(),
    )

    gridlines = ax.gridlines(
        draw_labels=True,
        color="black",
        linewidth=1,
        linestyle=":",
    )

    gridlines.xlabel_style = {"fontsize": 14}
    gridlines.ylabel_style = {"fontsize": 14}

    ax.add_feature(
        cfeature.LAND,
        edgecolor="black",
        linewidth=1.0,
        facecolor="none",
        zorder=100,
    )


# ---------------------------------------------------------------------
# Calculate p-values and plot each case
# ---------------------------------------------------------------------

statgood_dict = {}

for case_name, case_info in CASES.items():
    print(f"Processing {case_name}: {case_info['label']}")

    soil_moisture = soil_moisture_data[
        case_info["soil_moisture_key"]
    ]

    # The analysis function expects a list of response variable names.
    pvalues = glace_esa_heights_only.glace_esa_pvalues(
        soil_moisture,
        [RESPONSE_VARIABLE],
        analysis_time,
        NUMBER_OF_ENSEMBLES,
        NUMBER_OF_CELLS,
    )

    # Store the cell-wise p-values for later analysis.
    # The function returns shape (responses, cells).
    statgood = pvalues[0]
    statgood_dict[case_name] = statgood

    row, column = case_info["subplot"]
    ax = axes[row, column]

    ax.set_title(
        f"Soil moisture vs. {RESPONSE_VARIABLE}\n"
        f"{analysis_time:%m/%d} — {case_info['label']}",
        fontsize=16,
    )

    contour = ax.tricontourf(
        longitudes,
        latitudes,
        statgood,
        levels=np.linspace(0, 0.1, 11),
        cmap="plasma",
        transform=ccrs.PlateCarree(),
        extend="both",
    )

    fig.colorbar(contour, ax=ax)


# ---------------------------------------------------------------------
# Save figure and p-value data
# ---------------------------------------------------------------------

PLOT_DIRECTORY.mkdir(parents=True, exist_ok=True)
PICKLE_DIRECTORY.mkdir(parents=True, exist_ok=True)

date_string = analysis_time.strftime("%Y-%m-%d")

plot_file = PLOT_DIRECTORY / f"{date_string}.png"
pickle_file = PICKLE_DIRECTORY / (
    f"{analysis_time:%Y-%m-%d}.pkl"
)

fig.savefig(plot_file, dpi=300)
plt.close(fig)

with pickle_file.open("wb") as file:
    pickle.dump(statgood_dict, file)

print(f"Saved figure: {plot_file}")
print(f"Saved p-values: {pickle_file}")