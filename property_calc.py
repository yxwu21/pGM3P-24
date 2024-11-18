import re
import glob
import numpy as np
import pickle

from pGM_analysis import extract_data, calculate_properties
from data_extractor import DataExtractor
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.misc import derivative


class CpStepExtractor(DataExtractor):

    def extract_number(self, file_path):
        eptot_values = []
        volume_values = []
        volume_regex = re.compile(r"VOLUME\s+=\s+([-+]?\d*\.\d+|\d+)")
        eptot_regex = re.compile(r"EPtot\s+=\s+([-+]?\d*\.\d+|\d+)")
        averages_flag = re.compile(
            r"A V E R A G E S\s+O V E R\s+\d+\s+S T E P S\s+(.+?)(?:\n\s*\n|\Z)"
        )

        with open(file_path, "r") as file:
            for line in file:
                if averages_flag.search(line):
                    break
                eptot_match = eptot_regex.search(line)
                volume_match = volume_regex.search(line)
                if eptot_match:
                    eptot_values.append(float(eptot_match.group(1)))
                if volume_match:
                    volume_values.append(float(volume_match.group(1)))

        # Limiting to the first 60 values
        return eptot_values[:60], volume_values[:60]

    def process_temp_folder(self, temp_folder):
        """Process a temperature folder and extract relevant data."""
        temp_folder_path = Path(temp_folder)
        temp = temp_folder_path.name.split("-")[-1]
        prod_folder = temp_folder_path

        file_paths = glob.glob(f"{prod_folder}/{self.file_pattern}")

        eptot_data = []
        volume_data = []
        for file_path in file_paths:
            eptot_values, volume_values = self.extract_number(file_path)
            if eptot_values:  # Ensure eptot_values is not empty
                eptot_data.extend(eptot_values)  # Use extend to add all values
            if volume_values:  # Ensure volume_values is not empty
                volume_data.extend(volume_values)  # Use extend to add all values

        return temp, eptot_data, volume_data

    def process_folder(self, folder):
        """Process each folder and extract data from temperature subfolders."""
        folder_name = Path(folder).name
        self.initialize_nested_dict(folder_name)

        temp_folders = self.get_file_paths(folder)
        for temp_folder in temp_folders:
            temp, eptot_data, volume_data = self.process_temp_folder(temp_folder)
            temp = float(temp)  # Convert temperature to float

            # Initialize if not already done
            if temp not in self.temp_data[folder_name][self.case]:
                self.temp_data[folder_name][self.case][temp] = {
                    "eptot": [],
                    "volume": [],
                }

            # Append data to the existing lists
            if eptot_data:
                self.temp_data[folder_name][self.case][temp]["eptot"].extend(eptot_data)

            if volume_data:
                self.temp_data[folder_name][self.case][temp]["volume"].extend(
                    volume_data
                )
            else:
                print(f"No valid data found in {temp_folder}.")


def calculate_quantum_vibrational_heat_capacity(T):
    # Constants
    h = 6.62607015e-34  # J⋅s
    k_B = 1.380649e-23  # J⋅K^-1
    c = 299792458  # m/s (speed of light)
    N_A = 6.022e23  # Avogadro's number

    # Vibrational frequencies in cm^-1
    v_i = np.array([3490, 3280, 1645])

    # Convert frequencies to s^-1
    v_i = v_i * c * 100

    # Compute total heat capacity contribution
    C_v_total = np.sum(
        [
            ((h * v) ** 2 / (k_B * T**2))
            * (np.exp(h * v / (k_B * T)) / (np.exp(h * v / (k_B * T)) - 1) ** 2)
            for v in v_i
        ]
    )

    # Convert to kcal/mol-K
    C_v_total_kcal_mol_K = C_v_total * N_A / 4184

    return C_v_total_kcal_mol_K


def calculate_corrected_heat_capacity(T, correction_term, energies_kcal_per_mol):
    # Convert the list to a numpy array if it's not already
    energies_kcal_per_mol = np.array(energies_kcal_per_mol)

    # Calculate mean energy and mean squared energy
    mean_energy = np.mean(energies_kcal_per_mol)
    mean_energy_squared = np.mean(energies_kcal_per_mol**2)

    # Calculate the fluctuation term
    fluctuation_term = mean_energy_squared - mean_energy**2

    # Boltzmann constant in kcal/(mol·K)
    k_B = 1.9872041e-3

    # Calculate heat capacity at constant pressure (C_P)
    C_P = fluctuation_term / (k_B * T**2) / 512 * 1000

    # Apply the correction term
    C_P_c = C_P + correction_term * 1000

    return C_P_c


def calculate_isothermal_compressibility(volumes, temperature, num_molecules=512):
    k_B = 0.0019872041  # kcal/(mol*K)
    k_B_joules = k_B * 4184  # Convert to J/(mol*K)
    N_A = 6.02214076e23

    volumes = np.array(volumes)
    volumes_m3 = volumes * 1e-30  # Å^3 to m^3
    mean_volume = np.mean(volumes_m3)
    mean_square_volume = np.mean(volumes_m3**2)

    kappa_T = (
        (mean_square_volume - mean_volume**2)
        * N_A
        / (k_B_joules * temperature * mean_volume)
    )  # Pa^-1 per molecule
    kappa_T *= 1e-5  # Convert Pa^-1 to bar^-1
    kappa_T_per_mol = kappa_T / num_molecules  # Convert per molecule to per mole

    kappa_T_per_mol *= 1e6  # Convert to [10^-6 bar^-1]

    return kappa_T_per_mol


def process_case_data(extractor, case, test_keys, temps, multi_run_data):
    """
    Process and extract EPtot and volume data for the given case and test keys.
    """
    step_data = extractor.run()

    for temp in temps:
        multi_run_data[float(temp)] = {
            "eptot": [],
            "volume": [],
        }  # Initialize the lists

        for test_key in test_keys:
            if test_key in step_data and case in step_data[test_key]:
                tmp_eptot = step_data[test_key][case][float(temp)]["eptot"]
                tmp_volume = step_data[test_key][case][float(temp)]["volume"]

                # Extend the lists to concatenate data from all test_keys
                multi_run_data[float(temp)]["eptot"].extend(tmp_eptot)
                multi_run_data[float(temp)]["volume"].extend(tmp_volume)

    return multi_run_data


def calculate_properties_for_temps(temps, multi_run_data, property_data, case):
    """
    Calculate heat capacity and isothermal compressibility for each temperature
    and store them in the property_data dictionary.
    """
    for temp in temps:
        # Perform the calculation: eptot_data + volume_data * 0.0144
        eptot_data = multi_run_data[float(temp)]["eptot"]
        volume_data = multi_run_data[float(temp)]["volume"]
        calculated_data = [e + v * 0.0144 for e, v in zip(eptot_data, volume_data)]

        qm_correction = calculate_quantum_vibrational_heat_capacity(temp)
        C_P_c = calculate_corrected_heat_capacity(temp, qm_correction, calculated_data)

        print(f"{case} - {temp}: Corrected Heat Capacity: {C_P_c} cal/(mol·K)")

        kappa_T = calculate_isothermal_compressibility(volume_data, temp)
        print(
            f"{case} - {temp}: Isothermal compressibility (κ_T) per mole = {kappa_T} [10^-6 bar^-1]"
        )

        # Save the calculated properties
        property_data[case]["Cp"][float(temp)] = C_P_c
        property_data[case]["Kt"][float(temp)] = kappa_T


if __name__ == "__main__":
    proj_path = "/home8/yxwu/pGM_water_model/MD_data-pGM-temp*"
    case = "vital-cosmos-120_340"
    wat_num = "512"

    eptot_file_pattern = "MD/3_Prod/*.out"
    eptot_extractor = CpStepExtractor(proj_path, case, wat_num, eptot_file_pattern)

    # Set temperature range and test keys
    temperatures = set(range(240, 381, 5))
    temps = [298] + list(temperatures)
    test_keys = ["MD_data-pGM-temp", "MD_data-pGM-temp-1", "MD_data-pGM-temp-3"]

    # Initialize multi_run_data for the main case
    multi_run_data = {}
    multi_run_data = process_case_data(
        eptot_extractor, case, test_keys, temps, multi_run_data
    )

    # Initialize property_data for the main case
    property_data = {case: {"Cp": {}, "Kt": {}}}

    # Calculate properties for the main case
    calculate_properties_for_temps(temps, multi_run_data, property_data, case)

    # List of benchmark cases and folder path
    bench_cases = ["TIP3P", "TIP4P", "TIP5P", "OPC", "OPC3", "SPCE", "TIP4PEW"]
    bench_folder_pattern = "MD_data_Benchamrk-pgm_mpi"
    bench_proj_path = f"/home8/yxwu/pGM_water_model/{bench_folder_pattern}"

    for bench_case in bench_cases:
        # Initialize property_data for the benchmark cases if not already done
        if bench_case not in property_data:
            property_data[bench_case] = {"Cp": {}, "Kt": {}}

        # EPtot extraction for the benchmark case
        bench_eptot_extractor = CpStepExtractor(
            bench_proj_path, bench_case, wat_num, eptot_file_pattern
        )

        # Initialize multi_run_data for the benchmark case
        bench_multi_run_data = {}
        bench_multi_run_data = process_case_data(
            bench_eptot_extractor, bench_case, test_keys, temps, bench_multi_run_data
        )

        # Calculate properties for the benchmark case
        calculate_properties_for_temps(
            temps, bench_multi_run_data, property_data, bench_case
        )

    print(property_data)
    # Save property_data to a file using pickle
    with open("property_data.pkl", "wb") as f:
        pickle.dump(property_data, f)
