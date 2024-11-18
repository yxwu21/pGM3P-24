from pathlib import Path
from src.result_analysis_func import (
    extract_dc_values_to_df,
    extract_dc_data,
)
from temp_result_analysis import read_expt_data
from matplotlib import pyplot as plt
import numpy as np
import pandas as pd


def analyze_diffusion_constants_with_temp(
    md_data_path, analysis_folder, case_name, data_file_name="DC.dat"
):
    results = {}
    md_path = Path(md_data_path)
    subdir = md_path / case_name

    folder_name = subdir.name
    analysis_path = subdir / "analysis" / analysis_folder
    file_path = analysis_path / data_file_name

    if file_path.exists():  # Check if the file exists before proceeding
        df = extract_dc_values_to_df(file_path)
        dc_results = extract_dc_data(df)
        results[folder_name] = dc_results
    else:
        results[folder_name] = "Data file not found"

    return results


def prepare_data(file_path):
    # Read experimental data
    data = read_expt_data(file_path)

    # Convert temperature to numeric
    data["Temperature"] = data["Temperature"].astype(float)

    # Check if temperature is in Celsius and convert to Kelvin if necessary
    if "celsius" in file_path.lower():
        temperatures = data["Temperature"] + 273.15
    else:
        temperatures = data["Temperature"]

    diffusions = data["Diffusion"]
    return temperatures, diffusions


def plot_data(file_path, linestyle, color, label, marker):
    temperatures, diffusions = prepare_data(file_path)
    plt.plot(
        temperatures,
        diffusions,
        linestyle=linestyle,
        color=color,
        label=label,
        linewidth=2,
        marker=marker,
    )


if __name__ == "__main__":
    temperatures2 = set(range(240, 381, 10))
    temperatures3 = set(range(240, 381, 5))

    interval = "5"
    # interval = "10"

    if interval == "5":
        temperatures = list(set(range(240, 381, 5)))
    else:
        temperatures = temperatures3 - temperatures2
    dc_data = {}
    for temp in temperatures:
        md_data_path = Path(
            f"/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-{temp}/nvt-npt-dipole"
        )

        analysis_folder = "diffusion_-11_-1"

        case_name = "vital-cosmos-120_340"
        results = analyze_diffusion_constants_with_temp(
            md_data_path,
            analysis_folder,
            case_name,
            data_file_name="DC.dat",
        )
        dc_data[temp] = results

    # Initialize lists to store the temperatures and corresponding 'r' diffusion constants
    temperatures = []
    diffusion_constants_r = []

    # Extract data from dc_data
    for temp, data in dc_data.items():
        r_diffusion_constant = data["vital-cosmos-120_340"]["Diffusion Constant"]["r"]
        temperatures.append(temp)
        diffusion_constants_r.append(r_diffusion_constant)

    # Sort the data by temperature
    sorted_temps_indices = sorted(
        range(len(temperatures)), key=lambda i: temperatures[i]
    )
    sorted_temperatures = [temperatures[i] for i in sorted_temps_indices]
    sorted_diffusion_constants_r = [
        diffusion_constants_r[i] for i in sorted_temps_indices
    ]

    # Plotting
    plt.figure(figsize=(8, 6))  # Set the figure size (optional)
    plt.plot(
        sorted_temperatures,
        sorted_diffusion_constants_r,
        marker="o",
        linestyle="-",
        color="blue",
        label=case_name,
        linewidth=2,
    )

    data_name_to_path = {
        "expt-Easteal_1989": "datasets/expt/expt_temp_diffusion-Easteal_1989.kelvin",
        "expt-Gillien_1972": "datasets/expt/expt_temp_diffusion-Gillien_1972.kelvin",
        "TIP3P": "datasets/ref_benchmark/bench_temp_diffusion.kelvin.TIP3P",
        "TIP4P": "datasets/ref_benchmark/bench_temp_diffusion.kelvin.TIP4P",
        "TIP5P": "datasets/ref_benchmark/bench_temp_diffusion.kelvin.TIP5P",
    }

    data_name_to_color = {
        "expt-Easteal_1989": "black",
        "expt-Gillien_1972": "purple",
        "TIP3P": "red",
        "TIP4P": "orange",
        "TIP5P": "gray",
    }

    data_name_to_marker = {
        "expt-Easteal_1989": None,
        "expt-Gillien_1972": None,
        "TIP3P": "^",
        "TIP4P": "d",
        "TIP5P": ">",
    }

    data_name_to_linestyle = {
        "expt-Easteal_1989": "dashdot",
        "expt-Gillien_1972": "dashed",
        "TIP3P": "--",
        "TIP4P": "--",
        "TIP5P": "--",
    }

    for data_name in data_name_to_path.keys():
        plot_data(
            file_path=data_name_to_path[data_name],
            linestyle=data_name_to_linestyle[data_name],
            color=data_name_to_color[data_name],
            label=data_name,
            marker=data_name_to_marker[data_name],
        )

    plt.xlabel("Temperature (K)", fontsize=18)
    plt.ylabel(r"$D(10^{-5} \, \text{cm}^2/\text{s})$", fontsize=18)
    plt.legend(fontsize=14)
    plt.tick_params(axis="both", which="major", labelsize=14)
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(f"outputs/fig/dc_vs_temperature_plot_{case_name}_{interval}_bench.png")
    plt.close()
