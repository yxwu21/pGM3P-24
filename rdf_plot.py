import re
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from diffusion_plot import merge_dictionaries
from scipy.optimize import curve_fit


def extract_rdf_data(file_path):
    rdf_data = {}
    with open(file_path, "r") as file:
        for line in file:
            # Skip comments and empty lines
            if line.startswith("#") or line.strip() == "":
                continue
            # Split the line into distance and value
            distance, value = map(float, line.split())
            rdf_data[distance] = value
    return rdf_data


def load_rdf_data(md_folder, cases, target_folder_pattern):
    temp_rdf = {case: {} for case in cases}

    md_folder = Path(md_folder)
    if not md_folder.exists():
        print(f"The directory {md_folder} does not exist.")
        return temp_rdf

    temp_folders = os.listdir(md_folder)

    for temp_folder in temp_folders:
        temp = temp_folder.split("-")[-1]
        for case in cases:
            target_folder = md_folder / temp_folder / "*" / case / target_folder_pattern
            output_files = glob.glob(f"{target_folder}/*.out")
            tmp_rdf = []
            for output_file in output_files:
                tmp_rdf_value = extract_rdf_data(output_file)
                tmp_rdf.append(tmp_rdf_value)
                # print(output_file)
                # print(tmp_rdf_value)

            if tmp_rdf:  # Ensure there is data before assigning
                temp_rdf[case][temp] = tmp_rdf_value

    return temp_rdf


def load_experimental_data(file_paths_with_names):
    expt_temp_rdf = {}
    for temp, file_path in file_paths_with_names.items():
        expt_temp_rdf[temp] = extract_rdf_data(file_path)
    return expt_temp_rdf


def plot_rdf_data(temp_rdf, cases, temperatures, save_folder, labels_markers):
    # Ensure the save folder exists
    save_folder_path = Path(save_folder)
    save_folder_path.mkdir(parents=True, exist_ok=True)

    for temp in temperatures:
        plt.figure(figsize=(8, 6))
        legend_handles = {}

        for case in cases:
            if temp in temp_rdf.get(case, {}):
                rdf_data = temp_rdf[case][temp]
                distances = list(rdf_data.keys())
                values = list(rdf_data.values())

                label, marker, color, linestyle, linewidth, marker_size = (
                    labels_markers.get(case, (case, "o", None, "-", 1, 12))
                )

                (handle,) = plt.plot(
                    distances,
                    values,
                    marker=marker,
                    color=color,
                    linestyle=linestyle,
                    linewidth=linewidth,
                    label=label,
                    markersize=marker_size,
                    markerfacecolor="none" if not color else color,
                    markeredgewidth=linewidth,
                )
                legend_handles[label] = handle

        plt.xlabel(r"Distance $(\AA)$", fontsize=18)
        plt.ylabel(r"$g(r_{\mathrm{OO}})$", fontsize=18)

        ordered_labels = [
            labels_markers[case][0] for case in cases if case in labels_markers
        ]
        ordered_handles = [
            legend_handles[label] for label in ordered_labels if label in legend_handles
        ]

        plt.legend(ordered_handles, ordered_labels, fontsize=12)
        plt.tick_params(axis="both", which="major", labelsize=12)
        plt.tight_layout()

        plot_path = save_folder_path / f"rdf_vs_distance_{temp}K.pdf"
        plt.savefig(plot_path)
        plt.close()


def calculate_rOO_for_all_models(merged_dict):
    rOO_results = {}

    # Loop over each model in merged_dict
    for model in merged_dict:
        if "298" in merged_dict[model]:
            rdf_data = merged_dict[model]["298"]  # Get RDF data at 298 K

            # Find the key (r) that corresponds to the maximum value (g(r))
            rOO = max(rdf_data, key=rdf_data.get)

            # Store the result in the dictionary
            rOO_results[model] = rOO
            print(model, ":", rOO)
        else:
            print(f"No data available for model {model} at 298 K")

    return rOO_results


def gaussian(x, a, b, c):
    return a * np.exp(-((x - b) ** 2) / (2 * c**2))


def calculate_precise_rOO(merged_dict, model, temperature="298"):
    if temperature not in merged_dict[model]:
        print(f"No data available for model {model} at {temperature} K")
        return None

    # Get the RDF data for the given model and temperature
    rdf_data = merged_dict[model][temperature]

    # Convert the RDF data to numpy arrays for fitting
    r_values = np.array(list(rdf_data.keys()))
    g_values = np.array(list(rdf_data.values()))

    # Find a rough estimate of the first peak for the initial guess
    rough_rOO = r_values[np.argmax(g_values)]

    # Limit the fitting to a region around the rough estimate of rOO
    # You can adjust the range to capture the peak properly (e.g., 0.5 units around the peak)
    fitting_region = (r_values > rough_rOO - 0.5) & (r_values < rough_rOO + 0.5)
    r_fit = r_values[fitting_region]
    g_fit = g_values[fitting_region]

    # Perform Gaussian fitting
    initial_guess = [
        max(g_fit),
        rough_rOO,
        0.1,
    ]  # Initial guesses for amplitude, mean, and standard deviation
    params, _ = curve_fit(gaussian, r_fit, g_fit, p0=initial_guess)

    # Extract the peak (mean of the Gaussian) as the precise rOO
    rOO_precise = params[1]

    print(model, ":", rOO_precise)

    # # Plot the RDF data and the Gaussian fit
    # plt.figure(figsize=(8, 6))
    # plt.plot(r_values, g_values, "o", label="RDF Data", color="blue")
    # plt.plot(r_fit, gaussian(r_fit, *params), "-", label="Gaussian Fit", color="red")
    # plt.axvline(
    #     x=rOO_precise,
    #     color="green",
    #     linestyle="--",
    #     label=f"Precise rOO = {rOO_precise:.3f}",
    # )
    # plt.xlabel("r (distance)")
    # plt.ylabel("g(r) (RDF)")
    # plt.title(f"RDF and Gaussian Fit for {model} at {temperature} K")
    # plt.legend()
    # plt.grid(True)
    # plt.show()

    return rOO_precise


if __name__ == "__main__":

    rdf_type = "rdf_goo"
    # rdf_type = "rdf_ghh"
    # rdf_type = "rdf_goh"

    cases = [
        # "vital-cosmos-120_320",
        "vital-cosmos-120_340",
        # "vital-cosmos-120_578",
    ]

    md_folder = Path("MD_data-0-MC")
    # md_folder = Path("MD_data-pGM-temp")
    target_folder_pattern = f"analysis/{rdf_type}*"
    temp_rdf = load_rdf_data(md_folder, cases, target_folder_pattern)
    # print(temp_rdf)

    bench_cases = [
        "TIP3P",
        # "TIP4P",
        "TIP5P",
        "OPC",
        "OPC3",
        "SPCE",
        "TIP4PEW",
    ]

    bench_folder = Path("MD_data_Benchmark_reorg")
    bench_target_folder_pattern = target_folder_pattern
    bench_temp_rdf = load_rdf_data(
        bench_folder, bench_cases, bench_target_folder_pattern
    )
    # print(bench_temp_rdf)

    # read experiment data
    file_paths_with_names = {
        "298": f"datasets/expt/rdf/expt_Soper_2013/analysis/{rdf_type}/{rdf_type}.out",
    }
    expt_temp_rdf = load_experimental_data(file_paths_with_names)
    # print(expt_temp_rdf)

    bench_temp_rdf["Expt"] = expt_temp_rdf

    merge_dict_list = [bench_temp_rdf, temp_rdf]
    merged_dict = merge_dictionaries(merge_dict_list)
    print(merged_dict)
    print(merged_dict.keys())
    print(merged_dict["vital-cosmos-120_340"]["298"])
    print("*" * 180)
    calculate_rOO_for_all_models(merged_dict)
    print("*" * 180)
    for model in merged_dict.keys():
        calculate_precise_rOO(merged_dict, model)

    target_case = "vital-cosmos-120_340"
    # save_folder = f"outputs/final_figs/{bench_folder}/{rdf_type}"
    # save_folder = f"outputs/acs_figs/{bench_folder}/{rdf_type}"
    save_folder = f"plots/{bench_folder}/{rdf_type}"

    labels_markers = {
        # "vital-cosmos-120_340": ("Vital Cosmos 120_340", None, "blue", "-", 1.5, 8),
        "TIP3P": ("TIP3P", None, "green", "--", 1.5, 8),
        # "TIP4P": ("TIP4P", None, "red", "--", 1.5, 8),
        "TIP5P": ("TIP5P", None, "purple", "--", 1.5, 8),
        "OPC": ("OPC", None, "y", "--", 1.5, 8),
        "OPC3": ("OPC3", None, "darkorange", "--", 1.5, 8),
        "SPCE": ("SPC/E", None, "deeppink", "--", 1.5, 8),
        "TIP4PEW": ("TIP4P-Ew", None, "red", "--", 1.5, 8),
        "Expt": ("Expt.", None, "black", "-", 1.5, 12),
        # "vital-cosmos-120_340": ("Vital Cosmos 120_340", None, "blue", "-", 1.5, 8),
        "vital-cosmos-120_340": ("pGM3P-24", None, "blue", "-", 1.5, 8),
    }

    temperatures = ["298"]
    plot_rdf_data(
        merged_dict,
        list(labels_markers.keys()),
        temperatures,
        save_folder,
        labels_markers,
    )
