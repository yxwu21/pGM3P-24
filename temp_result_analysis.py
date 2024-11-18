from typing import Literal, Union
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

from matplotlib.lines import lineStyles
from pathlib import Path
from src.result_analysis_func import (
    calculate_averages_single,
    find_matching_paths,
    analyze_diffusion_constants,
    calculate_averages_single_dipole,
)
from sklearn.metrics import mean_squared_error, mean_absolute_error
from result_analysis import extract_dipole_data


def read_expt_data(file_path):
    experimental_data = pd.read_csv(file_path, sep="\s+")
    return experimental_data


def calculate_errors(temperatures, densities, expt_temperatures, expt_densities):
    interpolated_expt_densities = np.interp(
        temperatures, expt_temperatures, expt_densities
    )
    mse = mean_squared_error(densities, interpolated_expt_densities)
    mae = mean_absolute_error(densities, interpolated_expt_densities)
    return mse, mae


def process_md_data(
    md_storage_path,
    md_data_pattern,
    input_dir_name,
    keys_to_average,
    output_pattern,
    dipole_folder_name: Union[str, None] = None,
):
    md_paths = list(md_storage_path.glob(md_data_pattern))
    temp_density = {}

    for md_path in md_paths:
        try:
            paths = find_matching_paths(md_path, input_dir_name)
            if dipole_folder_name is None:
                average_data = calculate_averages_single(
                    paths,
                    input_dir_name,
                    keys_to_average,
                    output_pattern,
                )
            else:
                # also read dipole
                average_data = calculate_averages_single_dipole(
                    paths,
                    input_dir_name,
                    keys_to_average,
                    output_pattern,
                    dipole_folder_name,
                )

            temp = int(md_path.parts[-2].split("-")[-1])
            temp_density[temp] = average_data
        except Exception as e:
            print(f"An error occurred: {e}")

    return temp_density


def process_md_data_with_dipole(
    md_storage_path, md_data_pattern, input_dir_name, keys_to_average, output_pattern
):
    md_paths = list(md_storage_path.glob(md_data_pattern))
    temp_dipole = {}

    for md_path in md_paths:
        try:
            paths = find_matching_paths(md_path, input_dir_name)
            average_data = calculate_averages_single_dipole(
                paths, input_dir_name, keys_to_average, output_pattern
            )
            temp = int(md_path.parts[-2].split("-")[-1])
            temp_dipole[temp] = average_data
        except Exception as e:
            print(f"An error occurred: {e}")

    return temp_dipole


def filter_temp_data(temp_density, excluded_temps):
    filtered_temp_data = {
        temp: {case: values for case, values in cases.items() if "Density" in values}
        for temp, cases in temp_density.items()
        if temp not in excluded_temps
    }
    return filtered_temp_data


def extract_data_to_tuple_with_key(
    filtered_temp_data, key: Literal["Density", "Dipole"] = "Density"
):
    case_data = {}
    for temp, cases in filtered_temp_data.items():
        for case, values in cases.items():
            if case not in case_data:
                case_data[case] = []
            case_data[case].append((int(temp), values[key]))
    return case_data


def plot_data(case_data: dict, expt_temperatures: dict, expt_densities: dict):
    plt.figure(figsize=(12, 8))
    plt.plot(
        expt_temperatures,
        expt_densities,
        linestyle="--",
        color="black",
        label="expt-Kell_1975",
    )

    error_metrics = {}
    for case, data in case_data.items():
        data_sorted = sorted(data)
        temperatures = [x[0] for x in data_sorted]
        densities = [x[1] for x in data_sorted]
        plt.plot(temperatures, densities, marker="o", label=case)
        mse, mae = calculate_errors(
            temperatures, densities, expt_temperatures, expt_densities
        )
        error_metrics[case] = {"MSE": mse, "MAE": mae}

    plt.xlabel("Temperature (K)")
    plt.ylabel("Density (g/cm³)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("outputs/fig/density_vs_temperature_plot.png")
    plt.close()

    return error_metrics


def filter_cases(error_metrics, metric, min_value, max_value):
    filtered_cases = {
        case: metrics
        for case, metrics in error_metrics.items()
        if min_value <= metrics[metric] <= max_value
    }
    return filtered_cases


def plot_filtered_data(
    case_data: dict,
    valid_case_key_to_metric: dict,
    expt_temperatures: Union[pd.Series, None],
    expt_densities: Union[pd.Series, None],
    metric: str,
    plot_type: Literal["Density", "Dipole"] = "Density",
):
    plot_type2ylabel = {"Density": "Density (g/cm³)", "Dipole": "Dipole"}

    plt.figure(figsize=(12, 8))
    if expt_temperatures is None or expt_densities is None:
        pass
    else:
        plt.plot(
            expt_temperatures,
            expt_densities,
            linestyle="--",
            color="black",
            label="expt-Kell_1975",
        )

    for case_key in valid_case_key_to_metric:
        data = case_data[case_key]
        data_sorted = sorted(data)
        temperatures = [x[0] for x in data_sorted]
        densities = [x[1] for x in data_sorted]
        plt.plot(temperatures, densities, marker="o", label=case_key)

    plt.xlabel("Temperature (K)")
    plt.ylabel(plot_type2ylabel[plot_type])
    # plt.ylim(0.9, 1.06)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"outputs/fig/filtered_{plot_type}_vs_temperature_plot_{metric}.png")
    plt.close()


def save_to_excel(
    case_data,
    error_metrics,
    original_temp_density,
    tuple_type: Literal["Density", "Dipole"] = "Density",
):
    data_for_df = {}
    for case, data in case_data.items():
        if case not in data_for_df:
            data_for_df[case] = {}
        for temp, density in sorted(data):
            data_for_df[case][f"T = {temp}K"] = density

    df = pd.DataFrame.from_dict(data_for_df, orient="index")
    df.reset_index(inplace=True)
    df.rename(columns={"index": "Case Name"}, inplace=True)
    df["MSE"] = df["Case Name"].map(lambda case: error_metrics[case]["MSE"])
    df["MAE"] = df["Case Name"].map(lambda case: error_metrics[case]["MAE"])

    # Append T=298K data
    if 298 in original_temp_density:
        for case in df["Case Name"]:
            df.loc[df["Case Name"] == case, f"T = 298K"] = (
                original_temp_density[298].get(case, {}).get(tuple_type, np.nan)
            )

    df.to_excel(f"outputs/excel/{tuple_type}_Temperature.xlsx", index=False)

    best_case_mse = min(error_metrics, key=lambda x: error_metrics[x]["MSE"])
    best_case_mae = min(error_metrics, key=lambda x: error_metrics[x]["MAE"])

    print(
        f"Best case based on MSE: {best_case_mse} with MSE = {error_metrics[best_case_mse]['MSE']}"
    )
    print(
        f"Best case based on MAE: {best_case_mae} with MAE = {error_metrics[best_case_mae]['MAE']}"
    )

    error_df = pd.DataFrame.from_dict(error_metrics, orient="index")
    error_df.reset_index(inplace=True)
    error_df.rename(columns={"index": "Case Name"}, inplace=True)
    error_df.to_excel(f"outputs/excel/{tuple_type}_Error_Metrics.xlsx", index=False)


def plot_trendlines_for_cases(
    case_data, cases_to_plot, expt_temperatures, expt_densities
):
    plt.figure(figsize=(12, 8))
    plt.plot(
        expt_temperatures,
        expt_densities,
        linestyle="--",
        color="black",
        label="expt-Kell_1975",
    )

    for case in cases_to_plot:
        if case in case_data:
            data = case_data[case]
            data_sorted = sorted(data)
            temperatures = [x[0] for x in data_sorted]
            densities = [x[1] for x in data_sorted]
            plt.plot(temperatures, densities, marker="o", label=case)

            # Fit a polynomial (degree 2 in this example)
            poly_coeffs = np.polyfit(temperatures, densities, 2)
            poly_fit = np.poly1d(poly_coeffs)
            trendline = poly_fit(temperatures)
            plt.plot(temperatures, trendline, linestyle="--", label=f"{case} trendline")

    plt.xlabel("Temperature (K)")
    plt.ylabel("Density (g/cm³)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("outputs/fig/trendlines_density_vs_temperature_plot.png")
    plt.close()


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

    densities = data["Density"]
    return temperatures, densities


def plot_trendlines_for_cases_with_benchmark(
    case_data, cases_to_plot, file_path, linestyle, color, label, marker
):
    temperatures, densities = prepare_data(file_path)
    plt.figure(figsize=(12, 8))
    plt.plot(
        temperatures,
        densities,
        linestyle=linestyle,
        color=color,
        label=label,
        linewidth=2,
        marker=marker,
    )

    for case in cases_to_plot:
        if case in case_data:
            data = case_data[case]
            data_sorted = sorted(data)
            temperatures = [x[0] for x in data_sorted]
            densities = [x[1] for x in data_sorted]
            plt.plot(
                temperatures,
                densities,
                marker="o",
                label=case,
                color="blue",
                linestyle="None",
            )

            # Fit a polynomial (degree 2 in this example)
            poly_coeffs = np.polyfit(temperatures, densities, 2)
            poly_fit = np.poly1d(poly_coeffs)
            trendline = poly_fit(temperatures)
            plt.plot(
                temperatures,
                trendline,
                linestyle="-",
                label=f"{case} trendline",
                color="blue",
                linewidth=2,
            )

    plt.xlabel("Temperature (K)", fontsize=18)
    plt.ylabel("Density (g/cm³)", fontsize=18)
    # plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.legend(fontsize=14)
    plt.tick_params(axis="both", which="major", labelsize=12)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("outputs/fig/trendlines_density_vs_temperature_plot_with_bench.png")
    plt.close()


# def plot_trendlines_for_cases_with_benchmark(
#     case_data,
#     cases_to_plot,
#     expt_temperatures,
#     expt_densities,
#     spc_temperatures,
#     spc_densities,
#     tip3p_temperatures,
#     tip3p_densities,
#     tip4p_temperatures,
#     tip4p_densities,
# ):
#     plt.figure(figsize=(8, 6))
#     plt.plot(
#         expt_temperatures,
#         expt_densities,
#         linestyle="dashdot",
#         color="black",
#         label="expt-Kell_1975",
#         linewidth=2,
#     )

#     plt.plot(
#         spc_temperatures,
#         spc_densities,
#         linestyle="--",
#         color="green",
#         label="SPC",
#         linewidth=2,
#         marker="s",
#     )

#     plt.plot(
#         tip3p_temperatures,
#         tip3p_densities,
#         linestyle="--",
#         color="red",
#         label="TIP3P",
#         linewidth=2,
#         marker="^",
#     )

#     plt.plot(
#         tip4p_temperatures,
#         tip4p_densities,
#         linestyle="--",
#         color="orange",
#         label="TIP4P",
#         linewidth=2,
#         marker="d",
#     )

#     for case in cases_to_plot:
#         if case in case_data:
#             data = case_data[case]
#             data_sorted = sorted(data)
#             temperatures = [x[0] for x in data_sorted]
#             densities = [x[1] for x in data_sorted]
#             plt.plot(
#                 temperatures,
#                 densities,
#                 marker="o",
#                 label=case,
#                 color="blue",
#                 linestyle="None",
#             )

#             # Fit a polynomial (degree 2 in this example)
#             poly_coeffs = np.polyfit(temperatures, densities, 2)
#             poly_fit = np.poly1d(poly_coeffs)
#             trendline = poly_fit(temperatures)
#             plt.plot(
#                 temperatures,
#                 trendline,
#                 linestyle="-",
#                 label=f"{case} trendline",
#                 color="blue",
#                 linewidth=2,
#             )

# plt.xlabel("Temperature (K)", fontsize=18)
# plt.ylabel("Density (g/cm³)", fontsize=18)
# # plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
# plt.legend(fontsize=14)
# plt.tick_params(axis="both", which="major", labelsize=12)
# plt.grid(True)
# plt.tight_layout()
# plt.savefig("outputs/fig/trendlines_density_vs_temperature_plot_with_bench.png")
# plt.close()


def extract_and_export_data(case_data, md_data_path, dipole_folder_name):
    diffusion_results = analyze_diffusion_constants(md_data_path, "diffusion_-11_-1")

    dipole_data = {}
    diffusion_constants = {}
    for case in case_data.keys():
        dipole_file_path = md_data_path / case / "MD" / dipole_folder_name / "fort.200"
        dipole_data[case] = extract_dipole_data(str(dipole_file_path))["t"]
        if case in diffusion_results:
            diffusion_constants[case] = diffusion_results[case]["Diffusion Constant"][
                "r"
            ]
        else:
            diffusion_constants[case] = np.nan

    return dipole_data, diffusion_constants


def plot_filtered_data_with_benchmark(
    case_data,
    error_metrics,
    expt_temperatures,
    expt_densities,
    metric,
    min_value,
    max_value,
    bench_case_data,
):
    filtered_cases = filter_cases(error_metrics, metric, min_value, max_value)

    plt.figure(figsize=(12, 8))
    plt.plot(
        expt_temperatures,
        expt_densities,
        linestyle="--",
        color="black",
        label="expt-Kell_1975",
    )

    for case in filtered_cases:
        data = case_data[case]
        data_sorted = sorted(data)
        temperatures = [x[0] for x in data_sorted]
        densities = [x[1] for x in data_sorted]
        plt.plot(temperatures, densities, marker="o", label=case)

    for case, data in bench_case_data.items():
        data = bench_case_data[case]
        data_sorted = sorted(data)
        temperatures = [x[0] for x in data_sorted]
        densities = [x[1] for x in data_sorted]
        plt.plot(temperatures, densities, linestyle="--", label=case, marker="^")

    plt.xlabel("Temperature (K)")
    plt.ylabel("Density (g/cm³)")
    # plt.ylim(0.9, 1.06)
    plt.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(
        f"outputs/fig/benchmark_filtered_density_vs_temperature_plot_{metric}.png"
    )
    plt.close()


def main():
    # read my md results
    md_storage_path = Path("/home8/yxwu/pGM_water_model/MD_data")
    md_data_pattern = "langevin_Berendsen-*/nvt-npt-dipole"
    input_dir_name = "1c-200ps_10iter_output_npt"
    keys_to_average = ["Density", "Etot"]
    output_pattern = "*.10.out"
    dipole_folder_name = "dipole_1c-1step_output_npt"

    temperatures2 = set(range(240, 381, 10))
    temperatures3 = set(range(240, 381, 5))

    # exclude_temperatures = temperatures3 - temperatures2
    # excluded_temps = [298] + list(exclude_temperatures)
    excluded_temps = [298]

    temp_density = process_md_data(
        md_storage_path,
        md_data_pattern,
        input_dir_name,
        keys_to_average,
        output_pattern,
        dipole_folder_name,
    )
    filtered_temp_data = filter_temp_data(temp_density, excluded_temps)
    case_data_of_density = extract_data_to_tuple_with_key(
        filtered_temp_data, key="Density"
    )
    case_data_of_dipole = extract_data_to_tuple_with_key(
        filtered_temp_data, key="Dipole"
    )

    # read benchmark data
    wat_num = "512"
    bench_md_storage_path = Path("/home8/yxwu/pGM_water_model/MD_data_Benchmark_reorg")
    bench_md_data_pattern = f"Benchmark-*/{wat_num}"
    bench_input_dir_name = "4_Prod"
    bench_temp_density = process_md_data(
        bench_md_storage_path,
        bench_md_data_pattern,
        bench_input_dir_name,
        keys_to_average,
        output_pattern,
    )
    bench_filtered_temp_data = filter_temp_data(bench_temp_density, excluded_temps)
    bench_case_data = extract_data_to_tuple_with_key(bench_filtered_temp_data)

    # read experiment data
    expt_data_file_path = "datasets/expt/expt_temp_density-IPTS_68.celsius.full"
    expt_data = read_expt_data(expt_data_file_path)
    if "celsius" in expt_data_file_path.lower():
        expt_temperatures = expt_data["Temperature"] + 273.15
    else:
        expt_temperatures = expt_data["Temperature"]
    expt_densities = expt_data["Density"]

    # # read benchmark data
    # spc_data_file_path = "datasets/ref_benchmark/bench_temp_density.celsius.SPC"
    # spc_data = read_expt_data(spc_data_file_path)
    # if "celsius" in spc_data_file_path.lower():
    #     spc_temperatures = spc_data["Temperature"] + 273.15
    # else:
    #     spc_temperatures = spc_data["Temperature"]
    # spc_densities = spc_data["Density"]

    # tip3p_data_file_path = "datasets/ref_benchmark/bench_temp_density.celsius.TIP3P"
    # tip3p_data = read_expt_data(tip3p_data_file_path)
    # if "celsius" in tip3p_data_file_path.lower():
    #     tip3p_temperatures = tip3p_data["Temperature"] + 273.15
    # else:
    #     tip3p_temperatures = tip3p_data["Temperature"]
    # tip3p_densities = tip3p_data["Density"]

    # tip4p_data_file_path = "datasets/ref_benchmark/bench_temp_density.celsius.TIP4P"
    # tip4p_data = read_expt_data(tip4p_data_file_path)
    # if "celsius" in tip4p_data_file_path.lower():
    #     tip4p_temperatures = tip4p_data["Temperature"] + 273.15
    # else:
    #     tip4p_temperatures = tip4p_data["Temperature"]
    # tip4p_densities = tip4p_data["Density"]

    # plot all lines
    # TODO: change it to plot dipole benchmark
    error_metrics = plot_data(case_data_of_density, expt_temperatures, expt_densities)
    save_to_excel(case_data_of_density, error_metrics, temp_density)
    save_to_excel(case_data_of_dipole, error_metrics, temp_density, tuple_type="Dipole")

    # Filter case data based on the density measurement
    # Define the range for filtering
    metric = "MAE"  # or "MAE"
    min_value = 0.0
    max_value = 0.01
    valid_case_key_to_metric = filter_cases(error_metrics, metric, min_value, max_value)

    # Plot filtered data by chosen metric
    plot_filtered_data(
        case_data_of_density,
        valid_case_key_to_metric,
        expt_temperatures,
        expt_densities,
        metric,
        plot_type="Density",
    )

    plot_filtered_data(
        case_data_of_dipole,
        valid_case_key_to_metric,
        expt_temperatures=None,
        expt_densities=None,
        metric=metric,
        plot_type="Dipole",
    )

    # plot_filtered_data_with_benchmark(
    #     case_data_of_density,
    #     error_metrics,
    #     expt_temperatures,
    #     expt_densities,
    #     metric,
    #     min_value,
    #     max_value,
    #     bench_case_data,
    # )

    # List of cases to plot with trendlines
    cases_to_plot = ["vital-cosmos-120_340"]
    plot_trendlines_for_cases(
        case_data_of_density, cases_to_plot, expt_temperatures, expt_densities
    )

    data_name_to_path = {
        "expt-Kell_1975": "datasets/expt/expt_temp_density-IPTS_68.celsius.full",
        "SPC": "datasets/ref_benchmark/bench_temp_density.celsius.SPC",
        "TIP3P": "datasets/ref_benchmark/bench_temp_density.celsius.TIP3P",
        "TIP4P": "datasets/ref_benchmark/bench_temp_density.celsius.TIP4P",
    }

    data_name_to_color = {
        "expt-Kell_1975": "black",
        "SPC": "green",
        "TIP3P": "red",
        "TIP4P": "orange",
    }

    data_name_to_linestyle = {
        "expt-Kell_1975": "--",
        "SPC": "--",
        "TIP3P": "--",
        "TIP4P": "--",
    }

    data_name_to_marker = {
        "expt-Kell_1975": None,
        "SPC": "s",
        "TIP3P": "^",
        "TIP4P": "d",
    }

    for data_name in data_name_to_path.keys():
        plot_trendlines_for_cases_with_benchmark(
            case_data_of_density,
            cases_to_plot,
            file_path=data_name_to_path[data_name],
            linestyle=data_name_to_linestyle[data_name],
            color=data_name_to_color[data_name],
            label=data_name,
            marker=data_name_to_marker[data_name],
        )

    temp_dipole = process_md_data_with_dipole(
        md_storage_path,
        md_data_pattern,
        input_dir_name,
        keys_to_average,
        output_pattern,
    )

    filtered_temp_data = filter_temp_data(temp_dipole, excluded_temps)
    dipole_case_data = extract_data_to_tuple_with_key(filtered_temp_data, key="Dipole")
    print(dipole_case_data)

    # plot_filtered_data_with_dipole(
    #     dipole_case_data,
    #     error_metrics,
    #     expt_temperatures,
    #     expt_densities,
    #     metric,
    #     min_value,
    #     max_value,
    # )


if __name__ == "__main__":
    main()
