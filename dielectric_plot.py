import re
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from pathlib import Path
from diffusion_plot import read_multiple_expt_data, merge_dictionaries


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


def plot_dielectric_vs_temperature(
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
    plt.ylabel("Static Dielectric", fontsize=18)

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
        plot_path = save_folder_path / "dielectric_vs_temperature-pgm-mpi-trend.pdf"
    else:
        plot_path = save_folder_path / "dielectric_vs_temperature-pgm-mpi.pdf"
    plt.savefig(plot_path)
    plt.close()


if __name__ == "__main__":

    # cases = [
    #     # "vital-cosmos-120_320",
    #     "vital-cosmos-120_340",
    #     # "vital-cosmos-120_578",
    # ]

    # md_folder = Path("MD_data-0-MC")
    # target_folder_pattern = "analysis/diffusion_-11_-1"
    # temp_diffusion = load_diffusion_data(md_folder, cases, target_folder_pattern)
    # # print(temp_diffusion)

    temp_file_paths_with_names = {
        # "vital-cosmos-120_340": "datasets/pgm_dielectric/vital-cosmos-120_340.txt",
        "vital-cosmos-120_340": "datasets/pgm_dielectric/pgm-mpi_vital-cosmos-120_340.txt",
    }

    temp_dielectric = read_multiple_expt_data(temp_file_paths_with_names)
    print(temp_dielectric)

    bench_cases = [
        "TIP3P",
        # "TIP4P",
        "TIP5P",
        "OPC",
        "OPC3",
        "SPCE",
        "TIP4PEW",
    ]

    bench_folder = "MD_data_Benchmark_reorg"

    file_paths_with_names = {
        "TIP3P": f"datasets/{bench_folder}/benchmark_analysis_dielectric/TIP3P/tip3p.txt",
        # "TIP4P": f"datasets/{bench_folder}/benchmark_analysis_dielectric/TIP4P/tip4p.txt",
        "TIP5P": f"datasets/{bench_folder}/benchmark_analysis_dielectric/TIP5P/tip5p.txt",
        "TIP4PEW": f"datasets/{bench_folder}/benchmark_analysis_dielectric/TIP4PEW/tip4pew.txt",
        "OPC": f"datasets/{bench_folder}/benchmark_analysis_dielectric/OPC/opc.dat",
        "OPC3": f"datasets/{bench_folder}/benchmark_analysis_dielectric/OPC3/opc3.txt",
        "SPCE": f"datasets/{bench_folder}/benchmark_analysis_dielectric/SPCE/spce.txt",
    }

    bench_temp_dielectric = read_multiple_expt_data(file_paths_with_names)

    print(bench_temp_dielectric)

    # read experiment data
    file_paths_with_names = {
        "Expt": "datasets/expt/expt_temp_dielectric-Maryott_1956.celsius",
    }

    expt_temp_dielectric = read_multiple_expt_data(file_paths_with_names)
    print(expt_temp_dielectric)

    merge_dict_list = [bench_temp_dielectric, expt_temp_dielectric, temp_dielectric]
    # merge_dict_list = [bench_temp_dielectric, expt_temp_dielectric]
    merged_dict = merge_dictionaries(merge_dict_list)
    print(merged_dict)
    print(merged_dict.keys())

    target_case = "vital-cosmos-120_340"
    save_folder = f"outputs/final_figs/{bench_folder}"

    labels_markers = {
        "TIP3P": ("TIP3P", "o", "green", "--", 1.5, 8),
        # "TIP4P": ("TIP4P", "^", "red", "--", 1.5, 8),
        "TIP5P": ("TIP5P", "s", "purple", "--", 1.5, 8),
        "OPC": ("OPC", "<", "y", "--", 1.5, 8),
        "OPC3": ("OPC3", ">", "darkorange", "--", 1.5, 8),
        "SPCE": ("SPC/E", "8", "deeppink", "--", 1.5, 8),
        "TIP4PEW": ("TIP4P-EW", "^", "red", "--", 1.5, 8),
        "Expt": ("Expt.", None, "black", "-", 1.5, 12),
        "vital-cosmos-120_340": ("Vital Cosmos 120_340", "D", "blue", "-", 1.5, 8),
    }

    # labels_markers = {
    #     "vital-cosmos-120_340": ("Vital Cosmos 120_340", None, "blue", "-", 3, 12),
    #     "TIP3P": ("TIP3P", None, "green", "--", 3, 12),
    #     "TIP4P": ("TIP4P", None, "red", "dotted", 3, 12),
    #     "TIP5P": ("TIP5P", None, "purple", "-.", 3, 12),
    #     "Expt.Gillien-1972": ("Expt.", "D", "gray", "-", 3, 12),
    #     "Expt.Easteal-1989": ("Expt.", "D", "black", "-", 3, 12),
    # }

    plot_dielectric_vs_temperature(
        merged_dict,
        labels_markers,
        save_folder,
        plot_trendline=True,
        hollow_markers=True,
    )
