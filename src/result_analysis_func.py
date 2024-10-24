import re
import pandas as pd
import glob
import numpy as np

from pathlib import Path
from matplotlib import pyplot as plt


def find_matching_paths(rst_origin, input_dir_name):
    rst_origin_path = Path(rst_origin)
    pattern = f"*/MD/{input_dir_name}"
    matching_paths = list(rst_origin_path.glob(pattern))
    if not matching_paths:
        raise ValueError(
            f"No matching paths found for the pattern {rst_origin_path}/{pattern}"
        )
    return matching_paths


def extract_subpath_after_key(full_path: Path, key_directory: str) -> str:
    path_parts = full_path.parts
    try:
        index = path_parts.index(key_directory)
    except ValueError:
        raise ValueError(f"Key directory '{key_directory}' not found in path")
    relevant_path_parts = path_parts[index + 1 :]
    relevant_path = "/".join(relevant_path_parts)
    return relevant_path


def extract_dc_values_to_df(file_path):
    # Read the data into a DataFrame, skipping comments and using whitespace as delimiter
    df = pd.read_csv(file_path, sep="\s+")
    return df


def extract_dc_data(df, labels=["r", "x", "y", "z"]):
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


def analyze_diffusion_constants(md_data_path, analysis_folder, data_file_name="DC.dat"):
    results = {}
    md_path = Path(md_data_path)
    subdirs = [d for d in md_path.iterdir() if d.is_dir()]

    for subdir in subdirs:
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


def extract_output_data(filepath):
    # Dictionary to store the extracted values
    data = {}

    # Start and end markers for the data section
    start_marker = "A V E R A G E S   O V E R"
    end_marker = (
        "------------------------------------------------------------------------------"
    )

    # Flag to determine if we are in the relevant section
    in_section = False

    with open(filepath, "r") as file:
        for line in file:
            # Check if the line is the start of the relevant section
            if start_marker in line:
                in_section = True
                continue

            # Check if the line is the end of the relevant section
            if end_marker in line and in_section:
                break

            # If in the section, process the line
            if in_section:
                # Normalize the line to ensure spaces around '=' for consistent splitting
                line = line.replace("=", " = ")
                parts = line.split()
                key = None
                for i, part in enumerate(parts):
                    if part == "=" and i > 0 and i + 1 < len(parts):
                        key = parts[i - 1]
                        value = parts[i + 1]
                        # Remove any trailing commas which might be present in the values
                        value = value.replace(",", "")
                        try:
                            # Try to convert the value to float if possible
                            data[key] = float(value)
                        except ValueError:
                            data[key] = value
    return data


def calculate_averages(paths, input_dir_name, keys_to_average):
    avg_output_data = {}

    for path in paths:
        path = Path(path)  # Ensure `path` is a Path object
        data_file_pattern = str(path / "*.out")
        data_paths = glob.glob(data_file_pattern)
        values_dict = {
            key: [] for key in keys_to_average
        }  # Dictionary to store lists of values for each key

        for data_path in data_paths:
            data = extract_output_data(data_path)
            for key in keys_to_average:
                if key in data:
                    values_dict[key].append(data[key])

        case_name = path.parents[1].name
        if case_name not in avg_output_data:
            avg_output_data[case_name] = {}

        for key, values in values_dict.items():
            if values:  # Check if the list is not empty
                avg_value = sum(values) / len(values)
                avg_output_data[case_name][key] = avg_value

    return avg_output_data


def calculate_averages_single(paths, input_dir_name, keys_to_average, output_pattern):
    avg_output_data = {}

    for path in paths:
        path = Path(path)  # Ensure `path` is a Path object
        data_file_pattern = str(path / output_pattern)
        data_paths = glob.glob(data_file_pattern)
        values_dict = {
            key: [] for key in keys_to_average
        }  # Dictionary to store lists of values for each key

        for data_path in data_paths:
            data = extract_output_data(data_path)
            for key in keys_to_average:
                if key in data:
                    values_dict[key].append(data[key])

        case_name = path.parents[1].name
        if case_name not in avg_output_data:
            avg_output_data[case_name] = {}

        for key, values in values_dict.items():
            if values:  # Check if the list is not empty
                avg_value = sum(values) / len(values)
                avg_output_data[case_name][key] = avg_value

    return avg_output_data


def calculate_averages_single_dipole(
    paths, input_dir_name, keys_to_average, output_pattern, dipole_folder_name
):
    avg_output_data = {}

    for path in paths:
        path = Path(path)  # Ensure `path` is a Path object
        data_file_pattern = str(path / output_pattern)
        data_paths = glob.glob(data_file_pattern)
        values_dict = {
            key: [] for key in keys_to_average
        }  # Dictionary to store lists of values for each key

        for data_path in data_paths:
            data = extract_output_data(data_path)
            for key in keys_to_average:
                if key in data:
                    values_dict[key].append(data[key])

        case_name = path.parents[1].name
        if case_name not in avg_output_data:
            avg_output_data[case_name] = {}

        for key, values in values_dict.items():
            if values:  # Check if the list is not empty
                avg_value = sum(values) / len(values)
                avg_output_data[case_name][key] = avg_value

        dipole_file_path = path.parent / dipole_folder_name / "fort.200"
        dipole_data = extract_dipole_data(str(dipole_file_path))
        avg_total_dipole = dipole_data["t"]
        avg_output_data[case_name]["Dipole"] = avg_total_dipole

    return avg_output_data


def plot_rdf_single(data: dict, rdf_type: str, rdf_name: str, fig_folder: Path):
    ylabels = {
        "rdf_goo": "goo",
        "rdf_goh": "goh",
        "rdf_ghh": "ghh",
    }
    ylimits = {
        "rdf_goo": (-0.1, 4),
        "rdf_goh": (-0.1, 2),
        "rdf_ghh": (-0.1, 1.75),
    }
    fig, ax = plt.subplots(figsize=(8, 6))
    distance_cutoff = 9

    # Define default styles
    plot_styles = {
        "rdf": {"label": rdf_name, "linestyle": "-", "color": "red"},
        "exp_rdf": {"label": "expt_Soper_2013", "linestyle": "--", "color": "blue"},
    }

    # Customize for specific rdf_type
    if rdf_type == "rdf_goo":
        plot_styles["exp_rdf"].update({"label": "expt_Skinner_2013", "color": "black"})
        plot_styles["exp_rdf_soper"] = {
            "label": "expt_Soper_2013",
            "linestyle": "--",
            "color": "blue",
        }

    # Convert and check data
    for key in list(data):
        data[key] = np.array(data[key]) if isinstance(data[key], list) else data[key]
        if data[key].ndim != 2 or data[key].shape[1] < 2:
            print(f"Data for {key} is not 2-dimensional or lacks sufficient columns.")
            data.pop(key)

    # Plotting
    for key, style in plot_styles.items():
        if key in data:
            dataset = data[key]
            filter_rdf = dataset[:, 0] < distance_cutoff
            y = (
                dataset[filter_rdf, 1] / 2
                if rdf_type == "rdf_ghh" and key == "rdf"
                else dataset[filter_rdf, 1]
            )
            ax.plot(
                dataset[filter_rdf, 0],
                y,
                label=style["label"],
                linestyle=style["linestyle"],
                color=style["color"],
            )
        else:
            print(f"Missing data for {key}.")

    ax.legend()
    ax.set_xlabel("#Distance (Ang)")
    ax.set_ylabel(ylabels.get(rdf_type, "Unknown Type"))
    if rdf_type in ylimits:
        ax.set_ylim(ylimits[rdf_type])

    output_path = fig_folder / f"{rdf_type}.pdf"
    fig.savefig(output_path)
    plt.close(fig)


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
