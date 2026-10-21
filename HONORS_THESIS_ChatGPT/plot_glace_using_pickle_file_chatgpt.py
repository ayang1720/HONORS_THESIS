
"""
Alex Yang — Honors Thesis
Generate daily six-panel maps of GLACE-ESA p-values.
"""

import pickle
from datetime import datetime, timedelta
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np


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

PLOT_LEVELS = np.linspace(0, 0.1, 11)


# ---------------------------------------------------------------------
# Case definitions
# ---------------------------------------------------------------------

CASES = {
    "case_1": {
        "subplot": (0, 0),
        "label": "Top soil — Great Plains",
    },
    "case_2": {
        "subplot": (0, 1),
        "label": "Bottom soil — Great Plains",
    },
    "case_3": {
        "subplot": (1, 0),
        "label": "Top soil — California",
    },
    "case_4": {
        "subplot": (1, 1),
        "label": "Bottom soil — California",
    },
    "case_5": {
        "subplot": (2, 0),
        "label": "Top soil — Florida",
    },
    "case_6": {
        "subplot": (2, 1),
        "label": "Bottom soil — Florida",
    },
}


# ---------------------------------------------------------------------
# Load geographic coordinates
# ---------------------------------------------------------------------

with ESA_FILE.open("rb") as file:
    esa_data = pickle.load(file)

latitudes = esa_data["lat_1d"]
longitudes = esa_data["lon_1d"]


# ---------------------------------------------------------------------
# Plotting functions
# ---------------------------------------------------------------------

def create_map_axes():
    """Create and format the six map panels."""
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
            linewidth=1,
            facecolor="none",
            zorder=100,
        )

    return fig, axes


def plot_case(ax, pvalues, case_name, case_info, analysis_time):
    """Plot p-values for one case on the specified map axis."""
    ax.set_title(
        f"Soil moisture on 07/14 vs. {RESPONSE_VARIABLE}\n"
        f"{analysis_time:%m/%d} — {case_info['label']}",
        fontsize=16,
    )

    contour = ax.tricontourf(
        longitudes,
        latitudes,
        pvalues,
        levels=PLOT_LEVELS,
        cmap="plasma",
        transform=ccrs.PlateCarree(),
        extend="both",
    )

    return contour


def plot_day(analysis_time):
    """Load and plot all six cases for one day."""
    date_string = analysis_time.strftime("%Y-%m-%d")
    pickle_file = PICKLE_DIRECTORY / f"{date_string}.pkl"
    plot_file = PLOT_DIRECTORY / f"{date_string}.png"

    with pickle_file.open("rb") as file:
        statgood_dict = pickle.load(file)

    fig, axes = create_map_axes()

    for case_name, case_info in CASES.items():
        print(f"Plotting {case_name} for {date_string}")

        row, column = case_info["subplot"]
        ax = axes[row, column]

        contour = plot_case(
            ax=ax,
            pvalues=statgood_dict[case_name],
            case_name=case_name,
            case_info=case_info,
            analysis_time=analysis_time,
        )

        fig.colorbar(contour, ax=ax)

    fig.savefig(plot_file, dpi=300)
    plt.close(fig)

    print(f"Saved: {plot_file}")


# ---------------------------------------------------------------------
# Generate the daily time series
# ---------------------------------------------------------------------

PLOT_DIRECTORY.mkdir(parents=True, exist_ok=True)

current_time = START_TIME

while current_time <= END_TIME:
    plot_day(current_time)
    current_time += timedelta(days=1)