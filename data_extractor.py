import re
import glob
import numpy as np
import pickle
import os
import math
import MDAnalysis as mda

from pathlib import Path
from diffusion_plot import extract_dc_data, load_diffusion_data
from density_plot import load_data
from dipole_calc import calculate_average_dipole_moment_per_molecule


class DataExtractor:

    def __init__(self, proj_path, case, wat_num, file_pattern):
        self.proj_path = Path(proj_path)
        self.case = case
        self.wat_num = wat_num
        self.file_pattern = file_pattern
        self.temp_data = {}

    def extract_number(self, file_path):
        """Extract a number from a file. This method should be overridden by subclasses."""
        raise NotImplementedError("Subclasses should implement this method.")

    def get_file_paths(self, folder):
        """Get paths for each temperature folder."""
        return glob.glob(f"{folder}/{self.case}/{self.wat_num}/*")

    def process_temp_folder(self, temp_folder):
        """Process a temperature folder and extract relevant data."""
        temp_folder_path = Path(temp_folder)
        temp = temp_folder_path.name.split("-")[-1]
        prod_folder = temp_folder_path

        file_paths = glob.glob(f"{prod_folder}/{self.file_pattern}")

        # Avoid multiple calls to extract_number for the same file
        data_values = []
        for file_path in file_paths:
            extracted_value = self.extract_number(file_path)
            # debug flag
            print(file_path)
            print(extracted_value)
            if extracted_value is not None:
                data_values.append(extracted_value)

        return temp, data_values

    def process_folder(self, folder):
        """Process each folder and extract data from temperature subfolders."""
        folder_name = Path(folder).name
        self.initialize_nested_dict(folder_name)

        temp_folders = self.get_file_paths(folder)
        for temp_folder in temp_folders:
            temp, data_values = self.process_temp_folder(temp_folder)
            if data_values:
                avg_temp_value = np.average(data_values)
                self.temp_data[folder_name][self.case][float(temp)] = avg_temp_value
            else:
                print(f"No valid data found in {temp_folder}.")

    def initialize_nested_dict(self, folder_name):
        """Initialize nested dictionary for storing results."""
        if folder_name not in self.temp_data:
            self.temp_data[folder_name] = {}
        if self.case not in self.temp_data[folder_name]:
            self.temp_data[folder_name][self.case] = {}

    def run(self):
        """Run the data extraction for all folders matching the pattern."""
        folders = glob.glob(str(self.proj_path))
        for folder in folders:
            self.process_folder(folder)
        return self.temp_data


class DielectricExtractor(DataExtractor):

    def extract_number(self, file_path):
        """Extract the last number from a dielectric file, dividing by the average temperature and volume."""
        try:
            prod_folder_path = Path(file_path).parents[0]
            file_num = Path(file_path).name.split("_")[1]

            output_file = glob.glob(
                f"{prod_folder_path}/*{self.wat_num}wat.*{file_num}*.out"
            )[0]

            with open(output_file, "r") as file:
                file_content = file.read()

            section_pattern = r"A V E R A G E S   O V E R\s+[\d\sA-Z]+\n(.*?TEMP\(K\)\s*=\s*[\d.]+.*?VOLUME\s*=\s*[\d.]+)"
            section_match = re.search(section_pattern, file_content, re.DOTALL)

            if section_match:
                section_text = section_match.group(1)

                # Extract the average temperature
                temp_pattern = r"TEMP\(K\)\s*=\s*([\d.]+)"
                temp_match = re.search(temp_pattern, section_text)
                avg_temp = float(temp_match.group(1)) if temp_match else None

                # Extract the average volume
                volume_pattern = r"VOLUME\s*=\s*([\d.]+)"
                volume_match = re.search(volume_pattern, section_text)
                avg_volume = float(volume_match.group(1)) if volume_match else None

                # Print temperature and volume for debugging
                # print(output_file)
                # print(f"Temperature: {avg_temp}, Volume: {avg_volume}")

                if avg_temp is not None and avg_volume is not None:
                    # Extract the last number from the dielectric file
                    with open(file_path, "r") as file:
                        lines = file.readlines()
                        if lines:
                            last_line = lines[-1].strip()
                            last_number_str = last_line.split()[-1]
                            last_number = float(last_number_str)

                            # Return the last number divided by the product of temperature and volume
                            return last_number / (avg_temp * avg_volume)
                else:
                    print("Error: Temperature or Volume not found.")
            else:
                print("Error: Section pattern not found in the output file.")

        except (ValueError, FileNotFoundError, IndexError) as e:
            print(f"Error processing file {file_path}: {e}")

        return None


class DiffusionExtractor(DataExtractor):

    def extract_number(self, file_path):
        """Extract diffusion constant 'r' value."""
        try:
            diffusion_data = extract_dc_data(file_path)
            diffusion_constant_r = diffusion_data["Diffusion Constant"].get("r", None)
            return (
                float(diffusion_constant_r)
                if diffusion_constant_r is not None
                else None
            )
        except (ValueError, KeyError, FileNotFoundError) as e:
            print(f"Error processing file {file_path}: {e}")
        return None


class DensityExtractor(DataExtractor):

    def extract_number(self, file_path):
        """Extract density value from a file."""
        try:
            with open(file_path, "r") as file:
                file_content = file.read()

            section_pattern = (
                r"A V E R A G E S   O V E R\s+[\d\sA-Z]+\n(.*?Density\s*=\s*[\d.]+)"
            )
            section_match = re.search(section_pattern, file_content, re.DOTALL)
            pattern = r"Density\s*=\s*([\d.]+)"

            if section_match:
                density_match = re.search(pattern, section_match.group(1))
                return float(density_match.group(1)) if density_match else None
        except (ValueError, FileNotFoundError) as e:
            print(f"Error reading file {file_path}: {e}")
        return None


class EtotExtractor(DataExtractor):

    def extract_number(self, file_path):
        """Extract Etot value from a file."""
        try:
            with open(file_path, "r") as file:
                file_content = file.read()

            # Refined pattern to capture the entire section
            section_pattern = (
                r"A V E R A G E S\s+O V E R\s+\d+\s+S T E P S\s+(.+?)(?:\n\s*\n|\Z)"
            )
            section_match = re.search(section_pattern, file_content, re.DOTALL)

            if section_match:
                section_content = section_match.group(1)

                # Pattern to capture the Etot value
                etot_pattern = r"Etot\s*=\s*(-?[\d.]+)"
                value_match = re.search(etot_pattern, section_content)

                if value_match:
                    return float(value_match.group(1))

        except (ValueError, FileNotFoundError) as e:
            print(f"Error reading file {file_path}: {e}")


class EPtotExtractor(DataExtractor):

    def extract_number(self, file_path):
        """Extract EPtot value from a file."""
        try:
            with open(file_path, "r") as file:
                file_content = file.read()

            # Refined pattern to capture the entire section
            section_pattern = (
                r"A V E R A G E S\s+O V E R\s+\d+\s+S T E P S\s+(.+?)(?:\n\s*\n|\Z)"
            )
            section_match = re.search(section_pattern, file_content, re.DOTALL)

            if section_match:
                section_content = section_match.group(1)

                # Pattern to capture the Etot value
                eptot_pattern = r"EPtot\s*=\s*(-?[\d.]+)"
                value_match = re.search(eptot_pattern, section_content)

                if value_match:
                    return float(value_match.group(1))

        except (ValueError, FileNotFoundError) as e:
            print(f"Error reading file {file_path}: {e}")


class VolumeExtractor(DataExtractor):

    def extract_number(self, file_path):
        """Extract Volume value from a file."""
        try:
            with open(file_path, "r") as file:
                file_content = file.read()

            # Refined pattern to capture the entire section
            section_pattern = (
                r"A V E R A G E S\s+O V E R\s+\d+\s+S T E P S\s+(.+?)(?:\n\s*\n|\Z)"
            )
            section_match = re.search(section_pattern, file_content, re.DOTALL)

            if section_match:
                section_content = section_match.group(1)

                # Pattern to capture the Etot value
                volume_pattern = r"VOLUME\s*=\s*(-?[\d.]+)"
                value_match = re.search(volume_pattern, section_content)

                if value_match:
                    return float(value_match.group(1))

        except (ValueError, FileNotFoundError) as e:
            print(f"Error reading file {file_path}: {e}")


class VolumeflucExtractor(DataExtractor):

    def extract_number(self, file_path):
        """Extract Volume Fluctuation value from a file."""
        try:
            with open(file_path, "r") as file:
                file_content = file.read()

            # Refined pattern to capture the entire section
            section_pattern = r"R M S\s+F L U C T U A T I O N S\s*(.+?)(?:\n\s*\n|\Z)"
            section_match = re.search(section_pattern, file_content, re.DOTALL)

            if section_match:
                section_content = section_match.group(1)

                # Pattern to capture the Etot value
                volume_pattern = r"VOLUME\s*=\s*(-?[\d.]+)"
                value_match = re.search(volume_pattern, section_content)

                if value_match:
                    return float(value_match.group(1))

        except (ValueError, FileNotFoundError) as e:
            print(f"Error reading file {file_path}: {e}")


class DipoleCalculator(DataExtractor):

    def extract_number(self, file_paths):
        # print(file_paths)
        trajectory_files = sorted(file_paths)
        if file_paths:
            md_folder = Path(file_paths[0]).parents[1]
            # if "-298" in str(md_folder):
            prep_folder = Path(md_folder) / "prep_files"

            topology_file = glob.glob(f"{prep_folder}/*.prmtop")[0]

            universe = mda.Universe(topology_file, trajectory_files)
            mu_liquid = calculate_average_dipole_moment_per_molecule(universe)
            print(f"Average dipole moment in the liquid phase: {mu_liquid} Debye")

            return float(mu_liquid)

    def process_temp_folder(self, temp_folder):
        """Process a temperature folder and extract relevant data."""
        temp_folder_path = Path(temp_folder)
        temp = temp_folder_path.name.split("-")[-1]
        prod_folder = temp_folder_path

        print(temp, ":")

        file_paths = glob.glob(f"{prod_folder}/{self.file_pattern}")

        # Avoid multiple calls to extract_number for the same file
        data_values = []

        extracted_value = self.extract_number(file_paths)
        # debug flag
        # print(file_path)
        # print(extracted_value)
        if extracted_value is not None:
            data_values.append(extracted_value)

        return temp, data_values


# def calculate_kappa_T(avg_data):

#     # Constants
#     boltzmann_constant = 1.380649

#     # Initialize the property_data dictionary
#     property_data = {}

#     # Iterate over cases in avg_data
#     for case in avg_data:
#         if "Expt" not in case:
#             property_data[case] = {"kappa_T": {}}

#             # Iterate over temperatures for the current case
#             for temperature in avg_data[case]["volume"]:
#                 volume = avg_data[case]["volume"][temperature]
#                 volume_fluc = avg_data[case]["volume_fluctuation"][temperature]

#                 # Calculate heat capacity using the provided formula
#                 kappa_T = (volume_fluc * 1e4) / (
#                     volume * boltzmann_constant * temperature
#                 )

#                 # Save the result in property_data
#                 property_data[case]["kappa_T"][temperature] = kappa_T

#     return property_data


def calculate_kappa_T(avg_data):

    # Initialize the property_data dictionary
    property_data = {}

    # Iterate over cases in avg_data
    for case in avg_data:
        if "Expt" not in case:
            property_data[case] = {"kappa_T": {}}

            # Iterate over temperatures for the current case
            for temperature in avg_data[case]["volume"]:
                volume = avg_data[case]["volume"][temperature]
                volume_fluc = avg_data[case]["volume_fluctuation"][temperature]

                # Convert volume to m^3
                V = volume * 1e-30  # 1 Å^3 = 1e-30 m^3

                # Calculate <V^2> - <V>^2
                V_squared_diff = volume_fluc * 1e-60  # Convert to m^6

                # Boltzmann constant
                kB = 1.380649e-23  # J/K

                # Calculate κT in Pa^-1
                kappa_T = V_squared_diff / (kB * temperature * V)

                # Convert Pa^-1 to [10^-6 bar^-1]
                # 1 Pa^-1 = 10^5 bar^-1 = 10^11 [10^-6 bar^-1]
                kappa_T_bar = kappa_T * 1e11

                # Save the result in property_data
                property_data[case]["kappa_T"][temperature] = kappa_T_bar

    return property_data


# def calculate_kappa_T(avg_data):

#     # Constants
#     boltzmann_constant = 1.380649e-23  # in J/K (converted from kcal to J)
#     # kcal_to_joule = 4184  # 1 kcal = 4184 J
#     pa_to_bar = 1e-5  # 1 Pa = 10^-5 bar
#     angstrom_to_meter = 1e-10  # 1 Å = 10^-10 m
#     conversion_factor = angstrom_to_meter**3 * pa_to_bar * 1e6  # Å^3 to 10^-6 bar^-1

#     # Initialize the property_data dictionary
#     property_data = {}

#     # Iterate over cases in avg_data
#     for case in avg_data:
#         if "Expt" not in case:
#             property_data[case] = {"kappa_T": {}}

#             # Iterate over temperatures for the current case
#             for temperature in avg_data[case]["volume"]:
#                 volume = avg_data[case]["volume"][temperature]
#                 volume_fluc = avg_data[case]["volume_fluctuation"][temperature]

#                 # Convert volume to m^3 and calculate kappa_T
#                 kappa_T = (volume_fluc * 1e4) / (
#                     volume * boltzmann_constant * temperature
#                 )

#                 # Apply the unit conversion to 10^-6 bar^-1
#                 kappa_T *= conversion_factor

#                 # Save the result in property_data
#                 property_data[case]["kappa_T"][temperature] = kappa_T

#     return property_data


def average_property_data(all_data, property_key, case, test_keys):
    avg_data = {}
    std_dev = {}

    # Loop through all temperatures in the first test to get all temperature keys
    for temperature in all_data[property_key][test_keys[0]][case].keys():
        values = []

        try:
            # Collect values for each temperature across all test datasets
            for test_key in test_keys:
                if temperature in all_data[property_key][test_key][case]:
                    values.append(all_data[property_key][test_key][case][temperature])
        except KeyError as e:
            # Handle missing temperature key gracefully and print the error message
            print(
                f"KeyError: {e} - {test_key} or {case} or {temperature} may be missing in the data."
            )
            continue
        except Exception as e:
            # Catch any other exceptions and print the error message
            print(f"An error occurred: {e}")
            continue

        # Compute the average if values were found
        if values:
            avg_data[temperature] = np.mean(values)
            std_dev[temperature] = np.std(values, ddof=1) if len(values) > 1 else 0

    return avg_data, std_dev


def read_expt_data(file_path, name):
    """
    Reads experimental data from a file and returns it as a dictionary.
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


def calculate_thermal_expansion_coefficients(density_dict, P=1):
    temperatures = sorted(density_dict.keys())
    coefficients = {}

    for i in range(1, len(temperatures) - 1):
        T1 = temperatures[i - 1]
        T2 = temperatures[i + 1]
        T = temperatures[i]

        rho1 = density_dict[T1]
        rho2 = density_dict[T2]

        alpha_p = -(math.log(rho2) - math.log(rho1)) / (T2 - T1)
        coefficients[T] = alpha_p * 1e4  # Convert to 10^-4 K^-1

    return coefficients


if __name__ == "__main__":

    proj_path = "/home8/yxwu/pGM_water_model/MD_data-pGM-temp*"
    case = "vital-cosmos-120_340"
    wat_num = "512"

    # Initialize a unified dictionary to store all data
    all_data = {
        case: {
            "dielectric": {},
            "density": {},
            "diffusion": {},
            "etot": {},
            "eptot": {},
            "volume": {},
            "volume_fluctuation": {},
            "alpha_p": {},
            "dipole": {},
        }
    }

    # Dielectric extraction for the main case
    dielectric_file_pattern = "MD/3_Prod/fort.100_*"
    dielectric_extractor = DielectricExtractor(
        proj_path, case, wat_num, dielectric_file_pattern
    )
    dielectric_data = dielectric_extractor.run()
    all_data[case]["dielectric"] = dielectric_data

    # Density extraction for the main case
    density_file_pattern = "MD/3_Prod/*.out"
    density_extractor = DensityExtractor(proj_path, case, wat_num, density_file_pattern)
    density_data = density_extractor.run()
    all_data[case]["density"] = density_data

    # Diffusion extraction for the main case
    diffusion_file_pattern = "analysis/diffusion-last_10/*.dat"
    diffusion_extractor = DiffusionExtractor(
        proj_path, case, wat_num, diffusion_file_pattern
    )
    diffusion_data = diffusion_extractor.run()
    all_data[case]["diffusion"] = diffusion_data

    # Etot extraction for the main case
    etot_file_pattern = "MD/3_Prod/*.out"
    etot_extractor = EtotExtractor(proj_path, case, wat_num, etot_file_pattern)
    etot_data = etot_extractor.run()
    all_data[case]["etot"] = etot_data

    # EPtot extraction for the main case
    eptot_file_pattern = "MD/3_Prod/*.out"
    eptot_extractor = EPtotExtractor(proj_path, case, wat_num, eptot_file_pattern)
    eptot_data = eptot_extractor.run()
    all_data[case]["eptot"] = eptot_data

    # Volume extraction for the main case
    volume_file_pattern = "MD/3_Prod/*.out"
    volume_extractor = VolumeExtractor(proj_path, case, wat_num, volume_file_pattern)
    volume_data = volume_extractor.run()
    all_data[case]["volume"] = volume_data

    # Volume fluctuation extraction for the main case
    volumefluc_file_pattern = "MD/3_Prod/*.out"
    volumefluc_extractor = VolumeflucExtractor(
        proj_path, case, wat_num, volumefluc_file_pattern
    )
    volumefluc_data = volumefluc_extractor.run()
    all_data[case]["volume_fluctuation"] = volumefluc_data

    # Dipole calculation for the main case
    dipole_file_pattern = "MD/3_Prod/*.nc"
    dipole_extractor = DipoleCalculator(proj_path, case, wat_num, dipole_file_pattern)
    dipole_data = dipole_extractor.run()
    all_data[case]["dipole"] = dipole_data

    # Calculate thermal expansion coefficients
    for model in all_data.keys():
        for test_key in all_data[model]["density"].keys():
            density_dict = all_data[model]["density"][test_key][model]
            coefficients = calculate_thermal_expansion_coefficients(density_dict)
            all_data[model]["alpha_p"][test_key] = {}
            all_data[model]["alpha_p"][test_key][model] = coefficients

    # List of benchmark cases
    bench_cases = [
        "TIP3P",
        "TIP4P",
        "TIP5P",
        "OPC",
        "OPC3",
        "SPCE",
        "TIP4PEW",
    ]

    # bench_folder_pattern = "MD_data_Benchamrk-pgm_mpi*"
    bench_folder_pattern = "MD_data_Benchamrk-long-pgm_mpi-*"
    bench_proj_path = f"/home8/yxwu/pGM_water_model/{bench_folder_pattern}"

    # Loop through all benchmark cases for density and diffusion extraction
    for bench_case in bench_cases:
        # Initialize the benchmark entry
        all_data[bench_case] = {
            "density": {},
            "diffusion": {},
            "volume": {},
            "volume_fluctuation": {},
            "alpha_p": {},
            "dipole": {},
        }

        # Extract density data for the benchmark case
        bench_density_extractor = DensityExtractor(
            bench_proj_path, bench_case, wat_num, density_file_pattern
        )
        bench_density_data = bench_density_extractor.run()
        all_data[bench_case]["density"] = bench_density_data

        # Extract diffusion data for the benchmark case
        bench_diffusion_extractor = DiffusionExtractor(
            bench_proj_path, bench_case, wat_num, diffusion_file_pattern
        )
        bench_diffusion_data = bench_diffusion_extractor.run()
        all_data[bench_case]["diffusion"] = bench_diffusion_data

        # Extract volume data for the benchmark case
        bench_volume_extractor = VolumeExtractor(
            bench_proj_path, bench_case, wat_num, volume_file_pattern
        )
        bench_volume_data = bench_volume_extractor.run()
        all_data[bench_case]["volume"] = bench_volume_data

        # Extract volume fluctuation data for the benchmark case
        bench_volumefluc_extractor = VolumeflucExtractor(
            bench_proj_path, bench_case, wat_num, volumefluc_file_pattern
        )
        bench_volumefluc_data = bench_volumefluc_extractor.run()
        all_data[bench_case]["volume_fluctuation"] = bench_volumefluc_data

        # Dipole calculation for the main case
        bench_dipole_extractor = DipoleCalculator(
            bench_proj_path, bench_case, wat_num, dipole_file_pattern
        )
        bench_dipole_data = bench_dipole_extractor.run()
        all_data[bench_case]["dipole"] = bench_dipole_data

    # Calculate thermal expansion coefficients
    for model in bench_cases:
        for test_key in all_data[model]["density"].keys():
            density_dict = all_data[model]["density"][test_key][model]
            coefficients = calculate_thermal_expansion_coefficients(density_dict)
            all_data[model]["alpha_p"][test_key] = {}
            all_data[model]["alpha_p"][test_key][model] = coefficients

    # Print all data
    print("All Data:", all_data)
    print("*" * 180)

    with open("all_data.pkl", "wb") as f:
        pickle.dump(all_data, f)

    # with open("long_all_data.pkl", "wb") as f:
    #     pickle.dump(all_data, f)

    # with open("all_data.pkl", "rb") as f:
    #     all_data = pickle.load(f)

    avg_data = {case: {}}
    std_dev_data = {case: {}}

    test_keys = ["MD_data-pGM-temp", "MD_data-pGM-temp-1", "MD_data-pGM-temp-3"]
    bench_test_keys = ["MD_data_Benchamrk-pgm_mpi"]
    bench_test_keys = [
        "MD_data_Benchamrk-pgm_mpi",
        "MD_data_Benchamrk-pgm_mpi-1",
        "MD_data_Benchamrk-pgm_mpi-3",
    ]

    for property_key in all_data[case].keys():
        # Create the structure for the main case in avg_data
        avg_data[case][property_key], std_dev_data[case][property_key] = (
            average_property_data(all_data[case], property_key, case, test_keys)
        )

    for bench_case in bench_cases:
        avg_data[bench_case] = {}
        std_dev_data[bench_case] = {}
        for property_key in all_data[bench_case].keys():
            # Create the structure for the main case in avg_data
            (
                avg_data[bench_case][property_key],
                std_dev_data[bench_case][property_key],
            ) = average_property_data(
                all_data[bench_case],
                property_key,
                bench_case,
                bench_test_keys,
            )

    expt_path_with_names = {
        "density": "datasets/expt/expt_temp_density-IPTS_68.celsius.short",
        "diffusion": "datasets/expt/expt_temp_diffusion-merge.kelvin",
        "dielectric": "datasets/expt/expt_temp_dielectric-Maryott_1956.celsius",
        "alpha_p": "datasets/expt/expt_temp_thermal_expansion_coefficients-TIP4PEw.txt",
    }

    # read experiment data
    for property in expt_path_with_names.keys():

        file_path = expt_path_with_names[property]
        new_data = read_expt_data(file_path, name=property)

        if "Expt" not in avg_data:
            avg_data["Expt"] = {}

        avg_data["Expt"].update(new_data)

    dielectric_benchmark_ref = {
        "TIP3P": "datasets/MD_data_Benchmark_reorg/benchmark_analysis_dielectric/TIP3P/tip3p.txt",
        "TIP4P": "datasets/MD_data_Benchmark_reorg/benchmark_analysis_dielectric/TIP4P/tip4p.txt",
        "TIP5P": "datasets/MD_data_Benchmark_reorg/benchmark_analysis_dielectric/TIP5P/tip5p.txt",
        "OPC": "datasets/MD_data_Benchmark_reorg/benchmark_analysis_dielectric/OPC/opc.dat",
        "OPC3": "datasets/MD_data_Benchmark_reorg/benchmark_analysis_dielectric/OPC3/opc3.txt",
        "SPCE": "datasets/MD_data_Benchmark_reorg/benchmark_analysis_dielectric/SPCE/spce.txt",
        "TIP4PEW": "datasets/MD_data_Benchmark_reorg/benchmark_analysis_dielectric/TIP4PEW/tip4pew.txt",
    }

    for bench_case in dielectric_benchmark_ref.keys():

        file_path = dielectric_benchmark_ref[bench_case]
        new_data = read_expt_data(file_path, name="dielectric")

        avg_data[bench_case].update(new_data)

    print("Averaged Data:", avg_data)
    print("*" * 180)

    print(all_data.keys())
    print(avg_data.keys())

    # for property_key in all_data[case].keys():
    #     for test_key in test_keys:
    #         print(f"{test_key}:", all_data[case][property_key][test_key][case][298])

    with open("avg_data.pkl", "wb") as f:
        pickle.dump(avg_data, f)

    with open("std_dev_data.pkl", "wb") as f:
        pickle.dump(std_dev_data, f)

    # with open("long_avg_data.pkl", "wb") as f:
    #     pickle.dump(avg_data, f)

    # with open("long_std_dev_data.pkl", "wb") as f:
    #     pickle.dump(std_dev_data, f)

    # print(all_data)
    print(avg_data)
    print("*" * 180)

    berendsen_all_data = {case: {"dielectric": {}, "density": {}, "diffusion": {}}}

    cases = ["vital-cosmos-120_340"]
    berendsen_folders = ["MD_data-0", "MD_data-1", "MD_data-2", "MD_data-3"]
    for case in cases:
        for berendsen_folder in berendsen_folders:
            md_folder = Path(berendsen_folder)
            density_target_folder_pattern = "MD/1c-200ps_10iter_output_npt"
            temp_density = load_data(md_folder, cases, density_target_folder_pattern)
            berendsen_all_data[case]["density"][berendsen_folder] = temp_density

            diffusion_target_folder_pattern = "analysis/diffusion_-11_-1"
            temp_diffusion = load_diffusion_data(
                md_folder, cases, diffusion_target_folder_pattern
            )
            berendsen_all_data[case]["diffusion"][berendsen_folder] = temp_diffusion

    berendsen_avg_data = {case: {}}

    berendsen_test_keys = ["MD_data-0", "MD_data-1", "MD_data-2"]

    for property_key in berendsen_all_data[case].keys():
        if berendsen_all_data[case][property_key]:
            # Create the structure for the main case in berendsen_avg_data
            berendsen_avg_data[case][property_key], _ = average_property_data(
                berendsen_all_data[case], property_key, case, berendsen_test_keys
            )

    # print("All Data (Berendsen):", berendsen_all_data)
    # print("Averaged Data (Berendsen):", berendsen_avg_data)

    # print(avg_data["vital-cosmos-120_340"]["dielectric"])
    # print(f"{'Key':<10}{'Value':<20}")
    # print("-" * 30)
    # for key, value in avg_data["vital-cosmos-120_340"]["dielectric"].items():
    #     print(f"{key:<10}{value:<20}")

    avg_property_data = calculate_kappa_T(avg_data)
    print(avg_property_data)
    print("-" * 180)
    print(avg_data[case])
    print(avg_property_data[case])

    with open("avg_property_data.pkl", "wb") as f:
        pickle.dump(avg_property_data, f)

    print(avg_data["vital-cosmos-120_340"]["density"])
    print(all_data)
    print("*" * 180)
    print(std_dev_data)
