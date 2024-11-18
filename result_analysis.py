import pandas as pd
import json
import numpy as np
import glob
import re
import os

from pathlib import Path
from src.result_analysis_func import (
    analyze_diffusion_constants,
    extract_output_data,
    find_matching_paths,
    calculate_averages,
)


def load_json_params(json_filepath):
    """Load JSON data from a file."""
    try:
        with open(json_filepath, "r") as file:
            return json.load(file)
    except FileNotFoundError:
        print(f"Warning: No JSON file found at {json_filepath}")
        return {}


def find_dc_in_range(
    data, lower_bound=2.3, upper_bound=2.4, label="Diffusion Constant", set="r"
):
    results = {}
    # Iterate over each key and sub-dictionary in the main dictionary
    for key, value in data.items():
        # Check if the 'Diffusion Constant' for 'r' is within the specified range
        if lower_bound <= value[label][set] <= upper_bound:
            results[key] = value

    return results


def extract_dipole_data(file_path):
    with open(file_path, "r") as file:
        x_values, y_values, z_values, t_values = [], [], [], []

        # Process each line in the file
        for line in file:
            # Match lines containing the dipole moments
            if "total moment (x/y/z/t)" in line:
                # Extract numbers using regex
                numbers = list(
                    map(float, re.findall(r"-?\d+\.?\d*[Ee]?[+-]?\d*", line))
                )
                if (
                    len(numbers) == 5
                ):  # This should skip the index and use the four dipole values
                    _, x, y, z, t = numbers
                    x_values.append(x)
                    y_values.append(y)
                    z_values.append(z)
                    t_values.append(t)

    # Compute averages
    x_avg = sum(x_values) / len(x_values) if x_values else None
    y_avg = sum(y_values) / len(y_values) if y_values else None
    z_avg = sum(z_values) / len(z_values) if z_values else None
    t_avg = sum(t_values) / len(t_values) if t_values else None

    average_data = {"x": x_avg, "y": y_avg, "z": z_avg, "t": t_avg}

    return average_data


def export_data_to_excel(
    results,
    average_data,
    plot_lists,
    md_data_path,
    dipole_folder_name,
    output_excel_path,
):
    data = []
    for rdf_name in plot_lists:
        path = md_data_path / rdf_name
        json_file_path = path / "model_params.json"  # Adjust path as needed
        model_params = load_json_params(json_file_path)

        dipole_file_path = path / "MD" / dipole_folder_name / "fort.200"
        dipole_data = extract_dipole_data(str(dipole_file_path))
        avg_total_dipole = dipole_data["t"]

        diffusion_constant = (
            results.get(rdf_name, {}).get("Diffusion Constant", {}).get("r", "N/A")
        )
        avg_density = average_data.get(rdf_name, {}).get("Density", "N/A")
        avg_etot = average_data.get(rdf_name, {}).get("Etot", "N/A")

        # Combine all data into one dictionary for the current case
        case_data = {
            "filter_case": rdf_name,
            "Diffusion_constant": diffusion_constant,
            "Density": avg_density,
            "Etot": avg_etot,
            "Dipole": avg_total_dipole,
        }
        case_data.update(model_params)  # Add model parameters to the case data

        data.append(case_data)

    df = pd.DataFrame(data)
    df.to_excel(output_excel_path, index=False, engine="openpyxl")
    print(f"Data has been written to {output_excel_path}")


def organize_simulation_data(filter_results, md_data_path):
    output_dir = Path("pgm_wat_candidates")
    for filter_case in filter_results:
        candidate_path = output_dir / filter_case
        candidate_path.mkdir(exist_ok=True)
        prmtop_path = (
            md_data_path
            / filter_case
            / "MD"
            / "prep_files"
            / "case_1_pgm_512wat.prmtop"
        )
        rst_path = (
            md_data_path
            / filter_case
            / "MD"
            / "1a-10ps_output_nvt"
            / "case_1_pgm_512wat.nvt.pmemd-pgm.rst"
        )

        os.system(f"cp {prmtop_path} {candidate_path}")
        os.system(f"cp {rst_path} {candidate_path}")


if __name__ == "__main__":
    # md_data_path = Path(
    #     "/home/yxwu/pGM_water_model/MD_data/langevin_Berendsen_full-298/nvt-npt-dipole"
    # )

    # input_dir_name = "1c-200ps_10iter_output_npt"
    # keys_to_average = ["Density", "Etot"]

    # paths = find_matching_paths(md_data_path, input_dir_name)
    # average_data = calculate_averages(paths, input_dir_name, keys_to_average)

    # analysis_folder = "diffusion_-11_-1"

    # results = analyze_diffusion_constants(md_data_path, analysis_folder)
    # filter_results = find_dc_in_range(
    #     results,
    #     # lower_bound=2.3,
    #     # upper_bound=2.41,
    #     lower_bound=2.2,
    #     upper_bound=2.3,
    # )
    # print(results)
    # for filter_case in filter_results:
    #     diffusion_constant = filter_results[filter_case]["Diffusion Constant"]["r"]
    #     avg_density = average_data[filter_case]["Density"]
    #     avg_etot = average_data[filter_case]["Etot"]

    # #     print(f"{filter_case}:", filter_results[filter_case]["Diffusion Constant"]["r"])

    # print("Total:", len(filter_results.keys()))
    # # print(filter_results.keys())

    # categorized_filters = {}

    # for case in filter_results.keys():
    #     key, value = case.rsplit("_", 1)
    #     if key not in categorized_filters:
    #         categorized_filters[key] = []
    #     categorized_filters[key].append(int(value))

    # print(categorized_filters)

    # # # print(average_data)

    # """
    # export data to excel
    # """
    # dipole_folder_name = "dipole_1c-1step_output_npt"
    # output_excel_path = "outputs/excel/sec_filter_cases.xlsx"
    # plot_lists = filter_results.keys()
    # export_data_to_excel(
    #     results,
    #     average_data,
    #     plot_lists,
    #     md_data_path,
    #     dipole_folder_name,
    #     output_excel_path,
    # )

    # avg = extract_dipole_data("/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-250/nvt-npt-dipole/vital-cosmos-120_340/MD/dipole_1d-200ps_output_npt/fort.200")
    avg = extract_dipole_data("/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-250/nvt-npt-dipole/vital-cosmos-120_340/MD/dipole_1c-1step_output_npt/fort.200")
    print(avg)
