import re
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from diffusion_plot import merge_dictionaries


def read_expt_data(file_path, name):
    """
    Reads experimental data from a file and returns it as a dictionary.

    Parameters:
    file_path (str): Path to the experimental data file.

    Returns:
    dict: Dictionary with experimental temperatures (in Kelvin if 'celsius' is in the file path) as keys and densities as values.
    """
    expt_temp_density = {}
    convert_to_kelvin = "celsius" in file_path.lower()

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

    temp_dict = {name: expt_temp_density}
    return temp_dict


def extract_density(file_path):
    """
    Extracts the density value from the specified section of the given file.

    Parameters:
    file_path (str): The path to the file from which to extract the density value.

    Returns:
    str: The extracted density value or 'Not found' if the value is not found.
    """
    # Read the content of the file
    with open(file_path, "r") as file:
        file_content = file.read()

    # Define the pattern to match the section " A V E R A G E S   O V E R"
    section_pattern = (
        r"A V E R A G E S   O V E R\s+[\d\sA-Z]+\n(.*?Density\s*=\s*[\d.]+)"
    )

    # Search for the pattern in the file content
    section_match = re.search(section_pattern, file_content, re.DOTALL)

    # Define the pattern to extract the density value
    pattern = r"Density\s*=\s*([\d.]+)"

    # Extract the density value only if it is within the specified section
    if section_match:
        density_match = re.search(pattern, section_match.group(1))
        density_value = density_match.group(1) if density_match else "Not found"
    else:
        density_value = "Not found"

    return density_value


def calculate_average(densities):
    """
    Calculates the average of a list of densities.

    Parameters:
    densities (list): List of density values as strings.

    Returns:
    float: The average density.
    """
    # Convert densities to float and filter out 'Not found'
    numeric_densities = [float(d) for d in densities if d != "Not found"]
    if not numeric_densities:
        return "No valid densities found"
    return sum(numeric_densities) / len(numeric_densities)


def load_data(md_folder, cases, target_folder_pattern):
    """
    Loads the density data from the given MD folder for the specified cases.

    Parameters:
    md_folder (Path): The path to the MD data folder.
    cases (list): List of case names to process.

    Returns:
    dict: A dictionary containing the average densities for each case and temperature.
    """
    temp_density = {case: {} for case in cases}

    if not md_folder.exists():
        print(f"The directory {md_folder} does not exist.")
        return temp_density

    temp_folders = os.listdir(md_folder)

    for temp_folder in temp_folders:
        temp = temp_folder.split("-")[-1]
        for case in cases:
            target_folder = md_folder / temp_folder / "*" / case / target_folder_pattern
            output_files = glob.glob(f"{target_folder}/*.out")
            tmp_density = []
            for output_file in output_files:
                tmp_density_value = extract_density(output_file)
                tmp_density.append(tmp_density_value)
                # print(output_file)
                # print(tmp_density_value)

            tmp_avg_density = calculate_average(tmp_density)
            temp_density[case][float(temp)] = tmp_avg_density
            # print("-" * 160)

    return temp_density


# def plot_density_vs_temperature(data, labels_markers, save_folder):
#     """
#     Plots the density vs. temperature for the given data and saves the plot.

#     Parameters:
#     data (dict): Dictionary containing the density data for different cases.
#     labels_markers (dict): Dictionary containing the label and marker for each case.
#     save_folder (str): Folder path to save the plot.
#     """
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
#     plt.xlabel("Temperature (K)")
#     plt.ylabel("Density (g/cm^3)")
#     plt.title("Density vs. Temperature")
#     plt.legend()
#     plt.grid(True)

#     # Save the plot
#     plot_path = save_folder_path / "density_vs_temperature-MC-Langevin.pdf"
#     plt.savefig(plot_path)
#     plt.close()


def plot_density_vs_temperature(
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
    plt.ylabel(r"Density (g/cm$^3$)", fontsize=18)

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
        plot_path = save_folder_path / "density_vs_temperature-trend.pdf"
    else:
        plot_path = save_folder_path / "density_vs_temperature.pdf"
    plt.savefig(plot_path)
    plt.close()


if __name__ == "__main__":
    cases = [
        # "vital-cosmos-120_320",
        "vital-cosmos-120_340",
        # "vital-cosmos-120_578",
    ]

    md_folder = Path("MD_data-0-MC")
    target_folder_pattern = "MD/1c-200ps_10iter_output_npt"
    temp_density = load_data(md_folder, cases, target_folder_pattern)

    print(temp_density)
    # print(temp_density["vital-cosmos-120_340"]["298"])

    # bench_cases = [
    #     "TIP3P",
    #     # "TIP4P",
    #     "TIP5P",
    #     "OPC",
    #     "OPC3",
    #     "SPCE",
    #     "TIP4PEW",
    # ]

    # bench_folder = Path("MD_data_Benchmark_reorg")
    # bench_target_folder_pattern = "MD/4_Prod"
    # bench_temp_density = load_data(
    #     bench_folder, bench_cases, bench_target_folder_pattern
    # )

    # # print(bench_temp_density)

    # # read experiment data
    # expt_data_file_path = "datasets/expt/expt_temp_density-IPTS_68.celsius.short"
    # expt_temp_density = read_expt_data(expt_data_file_path, name="IPTS_68")
    # # print(expt_temp_density)

    # merge_dict_list = [bench_temp_density, expt_temp_density, temp_density]
    # merged_dict = merge_dictionaries(merge_dict_list)
    # # print(merged_dict)
    # # print(merged_dict.keys())

    # target_case = "vital-cosmos-120_340"
    # # save_folder = f"outputs/final_figs/{bench_folder}"
    # save_folder = f"outputs/acs_figs/{bench_folder}"

    # # labels_markers = {
    # #     "vital-cosmos-120_340": ("Vital Cosmos 120_340", "o", "blue", "-", 2),
    # #     "TIP3P": ("TIP3P", "x", "green", "--", 2),
    # #     "TIP4P": ("TIP4P", "^", "red", "--", 2),
    # #     "TIP5P": ("TIP5P", "s", "purple", "--", 2),
    # #     "IPTS_68": ("Expt", None, "black", "--", 2),
    # # }

    # labels_markers = {
    #     "TIP3P": ("TIP3P", "o", "green", "--", 1.5, 8),
    #     # "TIP4P": ("TIP4P", "^", "red", "--", 1.5, 8),
    #     "TIP5P": ("TIP5P", "s", "purple", "--", 1.5, 8),
    #     "OPC": ("OPC", "<", "y", "--", 1.5, 8),
    #     "OPC3": ("OPC3", ">", "darkorange", "--", 1.5, 8),
    #     "SPCE": ("SPC/E", "8", "deeppink", "--", 1.5, 8),
    #     "TIP4PEW": ("TIP4P-EW", "^", "red", "--", 1.5, 8),
    #     "IPTS_68": ("Expt.", None, "black", "-", 2, 12),
    #     # "vital-cosmos-120_340": ("Vital Cosmos 120_340", "D", "blue", "-", 1.5, 8),
    #     "vital-cosmos-120_340": ("pGM-water", "D", "blue", "-", 1.5, 8),
    # }

    # # plot_density_vs_temperature(merged_dict, labels_markers, save_folder)
    # plot_density_vs_temperature(
    #     merged_dict,
    #     labels_markers,
    #     save_folder,
    #     plot_trendline=False,
    #     hollow_markers=True,
    # )
