import re
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from src.result_analysis_func import extract_dc_values_to_df


def read_expt_data(file_path):
    expt_temp_density = {}
    convert_to_kelvin = "celsius" in file_path.lower()

    try:
        with open(file_path, "r") as file:
            next(file)  # Skip the header line
            for line in file:
                if line.strip():  # Skip empty lines
                    parts = line.split()
                    temp = float(parts[0])
                    density = float(parts[1])
                    if convert_to_kelvin:
                        temp = temp + 273.15  # Convert Celsius to Kelvin
                    expt_temp_density[temp] = density
    except Exception as e:
        print(f"Error reading data from {file_path}: {e}")
        return None

    return expt_temp_density


def read_multiple_expt_data(file_paths_with_names):
    merged_data = {}

    for name, file_path in file_paths_with_names.items():
        data = read_expt_data(file_path)
        if data is not None:
            merged_data[name] = data

    return merged_data


def extract_dc_data(file_path, labels=["r", "x", "y", "z"]):
    # Read the data into a DataFrame
    df = extract_dc_values_to_df(file_path)

    # Dictionary to store the results for each label under respective categories
    results = {
        "Diffusion Constant": {},
        "Slope": {},
        "Intercept": {},
        "Correlation": {},
    }

    # Pre-determine the relevant columns to avoid repeated searching
    relevant_columns = {
        "Diffusion Constant": next(
            (col for col in df.columns if col.endswith("[D]")), None
        ),
        "Slope": next((col for col in df.columns if col.endswith("[Slope]")), None),
        "Intercept": next(
            (col for col in df.columns if col.endswith("[Intercept]")), None
        ),
        "Correlation": next(
            (col for col in df.columns if col.endswith("[Corr]")), None
        ),
    }

    # Check for missing columns and inform the user
    for key, value in relevant_columns.items():
        if value is None:
            print(f"Warning: No column found for '{key}'. Check data integrity.")

    # Iterate over each label to extract data
    for index, label in enumerate(labels):
        for key, col_name in relevant_columns.items():
            if col_name:  # Ensure the column was found before attempting to access it
                results[key][label] = df.loc[index, col_name]
            else:
                results[key][label] = None  # Default value if column is not present

    return results


# def load_diffusion_data(md_folder, cases, target_folder_pattern):
#     temp_diffusion = {case: {} for case in cases}

#     if not md_folder.exists():
#         print(f"The directory {md_folder} does not exist.")
#         return temp_diffusion

#     temp_folders = os.listdir(md_folder)

#     for temp_folder in temp_folders:
#         temp = temp_folder.split("-")[-1]
#         for case in cases:
#             target_folder = md_folder / temp_folder / "*" / case / target_folder_pattern
#             output_files = glob.glob(f"{target_folder}/*.dat")
#             tmp_diffusion = []
#             for output_file in output_files:
#                 tmp_diffusion_value = extract_dc_data(output_file)[
#                     "Diffusion Constant"
#                 ]["r"]
#                 tmp_diffusion.append(tmp_diffusion_value)
#                 print(output_file)
#                 print(tmp_diffusion_value)

#             temp_diffusion[case][temp] = tmp_diffusion_value
#             # print("-" * 160)

#     return temp_diffusion


def load_diffusion_data(md_folder, cases, target_folder_pattern):
    # Initialize the dictionary to store diffusion data
    temp_diffusion = {case: {} for case in cases}

    # Check if the folder exists
    if not md_folder.exists():
        print(f"The directory {md_folder} does not exist.")
        return temp_diffusion

    # List the subdirectories (temperature folders)
    temp_folders = os.listdir(md_folder)

    # Loop through each temperature folder
    for temp_folder in temp_folders:
        temp = temp_folder.split("-")[-1]  # Extract temperature information
        # print(f"Processing temperature: {temp} (folder: {temp_folder})")  # Debugging

        # Loop through each case
        for case in cases:
            # Create the pattern to match the target folder and files
            target_folder = md_folder / temp_folder / "*" / case / target_folder_pattern
            output_files = glob.glob(f"{target_folder}/*.dat")

            # print(f"Case: {case}, Searching in: {target_folder}")  # Debugging
            # print(f"Found files: {output_files}")  # Debugging

            tmp_diffusion = []  # Temporary list to store diffusion values for this case
            if not output_files:
                print(
                    f"No files found for case {case} at temperature {temp}. Skipping."
                )  # Debugging
                continue

            # Loop through each output file and extract diffusion data
            for output_file in output_files:
                try:
                    tmp_diffusion_value = extract_dc_data(output_file)[
                        "Diffusion Constant"
                    ]["r"]
                    tmp_diffusion.append(tmp_diffusion_value)
                    # print(
                    #     f"File: {output_file}, Diffusion Value: {tmp_diffusion_value}"
                    # )  # Debugging
                except KeyError as e:
                    print(f"Error extracting data from {output_file}: {e}")  # Debugging
                except Exception as e:
                    print(
                        f"Unexpected error processing {output_file}: {e}"
                    )  # Debugging

            # Check if any diffusion value was extracted
            if tmp_diffusion:
                temp_diffusion[case][float(temp)] = tmp_diffusion_value
            else:
                print(
                    f"No valid diffusion values extracted for case {case} at temperature {temp}."
                )  # Debugging

    return temp_diffusion


def merge_dictionaries(dicts, selected_keys=None):
    merged_dict = {}

    # If no selected keys are provided, merge all keys
    if selected_keys is None or not selected_keys:
        selected_keys = set()
        for d in dicts:
            selected_keys.update(d.keys())

    for d in dicts:
        for key in selected_keys:
            if key in d:
                merged_dict[key] = d[key]

    return merged_dict


# def plot_diffusion_vs_temperature(data, labels_markers, save_folder):
#     # Ensure the save folder exists
#     save_folder_path = Path(save_folder)
#     save_folder_path.mkdir(parents=True, exist_ok=True)

#     plt.figure(figsize=(12, 8))

#     # Plot each case with customized label, marker, color, linestyle, and linewidth
#     for case, values in data.items():
#         temperatures = sorted(values.keys(), key=float)
#         densities = [values[temp] for temp in temperatures]

#         # Convert to numpy arrays
#         temperatures = np.array(temperatures, dtype=float)
#         densities = np.array(densities, dtype=float)

#         # Get customization options, with defaults
#         label, marker, color, linestyle, linewidth = labels_markers.get(
#             case, (case, "o", None, "-", 1)
#         )

#         print(
#             f"Plotting {case} with label '{label}', marker '{marker}', color '{color}', linestyle '{linestyle}', linewidth '{linewidth}'"
#         )
#         plt.plot(
#             temperatures,
#             densities,
#             marker=marker,
#             color=color,
#             linestyle=linestyle,
#             linewidth=linewidth,
#             label=label,
#         )

#     # Customize the plot
#     plt.xlabel("Temperature (K)", fontsize=18)
#     plt.ylabel(r"$D(10^{-5} \, \text{cm}^2/\text{s})$", fontsize=18)
#     plt.legend(fontsize=14)
#     plt.tick_params(axis="both", which="major", labelsize=14)
#     plt.tight_layout()
#     plt.grid(True)

#     # Save the plot
#     plot_path = save_folder_path / "diffusion_vs_temperature.pdf"
#     plt.savefig(plot_path)
#     plt.close()


# def plot_diffusion_vs_temperature(
#     data, labels_markers, save_folder, plot_trendline=False, hollow_markers=False
# ):
#     # Ensure the save folder exists
#     save_folder_path = Path(save_folder)
#     save_folder_path.mkdir(parents=True, exist_ok=True)

#     plt.figure(figsize=(8, 6))

#     # Plot each case with customized label, marker, color, linestyle, linewidth, and marker size
#     for case, values in data.items():
#         temperatures = sorted(values.keys(), key=float)
#         densities = [values[temp] for temp in temperatures]

#         # Convert to numpy arrays
#         temperatures = np.array(temperatures, dtype=float)
#         densities = np.array(densities, dtype=float)

#         # Get customization options, with defaults
#         label, marker, color, linestyle, linewidth, marker_size = labels_markers.get(
#             case, (case, "o", None, "-", 1, 12)
#         )

#         print(
#             f"Plotting {case} with label '{label}', marker '{marker}', color '{color}', linestyle '{linestyle}', linewidth '{linewidth}', marker_size '{marker_size}'"
#         )

#         if plot_trendline:
#             # Plot the original data points as scatter plot
#             plt.scatter(
#                 temperatures,
#                 densities,
#                 marker=marker,
#                 edgecolors=color,  # Use edgecolors for scatter plot
#                 facecolors="none" if hollow_markers else color,
#                 label=label,
#                 s=marker_size**2,
#             )
#             # Fit a linear model to the data and plot the trendline
#             z = np.polyfit(temperatures, densities, 3)
#             p = np.poly1d(z)
#             trendline_color = color if color else "black"
#             plt.plot(
#                 temperatures,
#                 p(temperatures),
#                 linestyle="-",
#                 color=trendline_color,
#                 linewidth=linewidth,
#             )
#         else:
#             # Plot the original data as specified in labels_markers
#             plt.plot(
#                 temperatures,
#                 densities,
#                 marker=marker,
#                 color=color,
#                 linestyle=linestyle,
#                 linewidth=linewidth,
#                 label=label,
#                 markersize=marker_size,
#                 markerfacecolor="none" if hollow_markers else color,
#                 markeredgewidth=linewidth,
#             )

#     # Customize the plot
#     plt.xlabel("Temperature (K)", fontsize=18)
#     plt.ylabel(r"$D(10^{-5} \, \text{cm}^2/\text{s})$", fontsize=18)
#     plt.legend(fontsize=14)
#     plt.tick_params(axis="both", which="major", labelsize=14)
#     plt.tight_layout()
#     # plt.grid(True)

#     # Save the plot
#     if plot_trendline:
#         plot_path = save_folder_path / "diffusion_vs_temperature-trend.pdf"
#     else:
#         plot_path = save_folder_path / "diffusion_vs_temperature.pdf"
#     plt.savefig(plot_path)
#     plt.close()


def plot_diffusion_vs_temperature(
    data, labels_markers, save_folder, plot_trendline=False, hollow_markers=False
):
    # Ensure the save folder exists
    save_folder_path = Path(save_folder)
    save_folder_path.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(8, 6))

    legend_handles = {}

    # Plot each case with customized label, marker, color, linestyle, linewidth, and marker size
    for case, values in data.items():
        temperatures = sorted(values.keys(), key=float)
        densities = [values[temp] for temp in temperatures]

        # Convert to numpy arrays
        temperatures = np.array(temperatures, dtype=float)
        densities = np.array(densities, dtype=float)

        # Get customization options, with defaults
        label, marker, color, linestyle, linewidth, marker_size = labels_markers.get(
            case, (case, "o", None, "-", 1, 12)
        )

        print(
            f"Plotting {case} with label '{label}', marker '{marker}', color '{color}', linestyle '{linestyle}', linewidth '{linewidth}', marker_size '{marker_size}'"
        )

        if plot_trendline:
            # Plot the original data points as scatter plot
            handle = plt.scatter(
                temperatures,
                densities,
                marker=marker,
                edgecolors=color,  # Use edgecolors for scatter plot
                facecolors="none" if hollow_markers else color,
                label=label,
                s=marker_size**2,
            )
            legend_handles[label] = handle
            # Fit a linear model to the data and plot the trendline
            z = np.polyfit(temperatures, densities, 3)
            p = np.poly1d(z)
            trendline_color = color if color else "black"
            plt.plot(
                temperatures,
                p(temperatures),
                linestyle="-",
                color=trendline_color,
                linewidth=linewidth,
            )
        else:
            # Plot the original data as specified in labels_markers
            (handle,) = plt.plot(
                temperatures,
                densities,
                marker=marker,
                color=color,
                linestyle=linestyle,
                linewidth=linewidth,
                label=label,
                markersize=marker_size,
                markerfacecolor="none" if hollow_markers else color,
                markeredgewidth=linewidth,
            )
            legend_handles[label] = handle

    # Customize the plot
    plt.xlabel("Temperature (K)", fontsize=18)
    plt.ylabel(r"$D(10^{-5} \, \text{cm}^2/\text{s})$", fontsize=18)

    # Reorder the legend according to labels_markers
    ordered_labels = [
        labels_markers[case][0] for case in labels_markers if case in labels_markers
    ]
    ordered_handles = [
        legend_handles[label] for label in ordered_labels if label in legend_handles
    ]

    plt.legend(ordered_handles, ordered_labels, fontsize=12)
    plt.tick_params(axis="both", which="major", labelsize=12)
    plt.tight_layout()
    # plt.grid(True)

    # Save the plot
    if plot_trendline:
        plot_path = save_folder_path / "diffusion_vs_temperature-trend.pdf"
    else:
        plot_path = save_folder_path / "diffusion_vs_temperature.pdf"
    plt.savefig(plot_path)
    plt.close()


if __name__ == "__main__":

    cases = [
        # "vital-cosmos-120_320",
        "vital-cosmos-120_340",
        # "vital-cosmos-120_578",
    ]

    md_folder = Path("MD_data-0-MC")
    target_folder_pattern = "analysis/diffusion_-11_-1"
    temp_diffusion = load_diffusion_data(md_folder, cases, target_folder_pattern)
    # print(temp_diffusion)

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
    bench_target_folder_pattern = "analysis/diffusion_-10_21"
    bench_temp_diffusion = load_diffusion_data(
        bench_folder, bench_cases, bench_target_folder_pattern
    )
    # print(bench_temp_diffusion)

    # read experiment data
    file_paths_with_names = {
        # "Expt.Gillien-1972": "datasets/expt/expt_temp_diffusion-Gillien_1972.kelvin",
        # "Expt.Easteal-1989": "datasets/expt/expt_temp_diffusion-Easteal_1989.kelvin",
        "Expt": "datasets/expt/expt_temp_diffusion-merge.kelvin",
    }

    expt_temp_diffusion = read_multiple_expt_data(file_paths_with_names)
    # print(expt_temp_diffusion)

    merge_dict_list = [bench_temp_diffusion, expt_temp_diffusion, temp_diffusion]
    merged_dict = merge_dictionaries(merge_dict_list)
    print(merged_dict)
    print(merged_dict.keys())

    target_case = "vital-cosmos-120_340"
    # save_folder = f"outputs/final_figs/{bench_folder}"
    save_folder = f"outputs/acs_figs/{bench_folder}"

    labels_markers = {
        "TIP3P": ("TIP3P", "o", "green", "--", 1.5, 8),
        # "TIP4P": ("TIP4P", "^", "red", "--", 1.5, 8),
        "TIP5P": ("TIP5P", "s", "purple", "--", 1.5, 8),
        "OPC": ("OPC", "<", "y", "--", 1.5, 8),
        "OPC3": ("OPC3", ">", "darkorange", "--", 1.5, 8),
        "SPCE": ("SPC/E", "8", "deeppink", "--", 1.5, 8),
        "TIP4PEW": ("TIP4P-EW", "^", "red", "--", 1.5, 8),
        # "Expt.Gillien-1972": ("Expt.", None, "gray", "--", 1.5, 8),
        # "Expt.Easteal-1989": ("Expt.", None, "black", "--", 1.5, 8),
        "Expt": ("Expt.", None, "black", "-", 1.5, 12),
        # "vital-cosmos-120_340": ("Vital Cosmos 120_340", "D", "blue", "-", 1.5, 8),
        "vital-cosmos-120_340": ("pGM-water", "D", "blue", "-", 1.5, 8),
    }

    # labels_markers = {
    #     "vital-cosmos-120_340": ("Vital Cosmos 120_340", None, "blue", "-", 3, 12),
    #     "TIP3P": ("TIP3P", None, "green", "--", 3, 12),
    #     "TIP4P": ("TIP4P", None, "red", "dotted", 3, 12),
    #     "TIP5P": ("TIP5P", None, "purple", "-.", 3, 12),
    #     "Expt.Gillien-1972": ("Expt.", "D", "gray", "-", 3, 12),
    #     "Expt.Easteal-1989": ("Expt.", "D", "black", "-", 3, 12),
    # }

    plot_diffusion_vs_temperature(
        merged_dict,
        labels_markers,
        save_folder,
        plot_trendline=False,
        hollow_markers=True,
    )
