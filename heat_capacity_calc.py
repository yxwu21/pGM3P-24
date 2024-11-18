import re
import glob
import numpy as np
import pickle
import math

from pGM_analysis import extract_data, calculate_properties
from data_extractor import DataExtractor
from pathlib import Path
from scipy.optimize import curve_fit
from scipy.misc import derivative
from data_extractor import read_expt_data
from scipy.integrate import simpson


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

    # Apply the correction term # in cal/(K*mol)
    C_P_c = C_P + correction_term * 1000

    return C_P_c


# def calculate_isothermal_compressibility(volumes, temperature, num_molecules=512):
#     k_B = 0.0019872041  # kcal/(mol*K)
#     k_B_joules = k_B * 4184  # Convert to J/(mol*K)
#     N_A = 6.02214076e23

#     volumes = np.array(volumes)
#     volumes_m3 = volumes * 1e-30  # Å^3 to m^3
#     mean_volume = np.mean(volumes_m3)
#     mean_square_volume = np.mean(volumes_m3**2)

#     kappa_T = (
#         (mean_square_volume - mean_volume**2)
#         * N_A
#         / (k_B_joules * temperature * mean_volume)
#     )  # Pa^-1 per molecule
#     kappa_T *= 1e-5  # Convert Pa^-1 to bar^-1
#     kappa_T_per_mol = kappa_T / num_molecules  # Convert per molecule to per mole

#     kappa_T_per_mol *= 1e6  # Convert to [10^-6 bar^-1]

#     return kappa_T_per_mol


def calculate_isothermal_compressibility(volumes, temperature, num_molecules=512):
    kB = 1.380649e-23  # J/K

    volumes = np.array(volumes)
    volumes_m3 = volumes * 1e-30  # Å^3 to m^3
    mean_volume = np.mean(volumes_m3)
    mean_square_volume = np.mean(volumes_m3**2)

    volume_fluc = mean_square_volume - mean_volume**2

    # Convert volume to m^3
    V = mean_volume

    # Calculate <V^2> - <V>^2
    V_squared_diff = volume_fluc

    # Calculate κT in Pa^-1
    kappa_T = V_squared_diff / (kB * temperature * V)

    # Convert Pa^-1 to [10^-6 bar^-1]
    # 1 Pa^-1 = 10^5 bar^-1 = 10^11 [10^-6 bar^-1]
    kappa_T_bar = kappa_T * 1e11

    return kappa_T_bar


# def calculate_delta_H_vap(T, U_liquid, V_liquid, pgm=False):
#     # Constants
#     N = 512  # Number of water molecules
#     p = 1  # Pressure in bar
#     R = 0.001987  # Gas constant in kcal/(mol·K)
#     if pgm:
#         U_gas = -976.4144
#     else:
#         U_gas = 0

#     # Calculate Cvib and Cni
#     Cvib = 0.002044 * T - 0.704656
#     Cni = -0.000569 * T + 0.133952

#     if pgm:
#         C = Cvib
#     else:
#         C = Cvib + Cni

#     # Convert V_liquid from Angstrom^3 to L/mol
#     V_liquid_L_mol = V_liquid * 1e-27 * 6.022e23 / 1000 / N

#     # Assume U_liquid is already in kcal for the whole system
#     # Convert to kcal/mol by dividing by N
#     U_liquid_per_mol = U_liquid / N

#     # Convert p*V term to kcal/mol
#     # 1 L·bar = 0.0241838 kcal/mol
#     pV_term = p * V_liquid_L_mol * 0.0241838
#     # pV_term = p * V_liquid * 1.44 * 10e-5

#     E_vib = calculate_quantum_vibrational_heat_capacity(T)
#     print("E_vib:", E_vib)
#     # E_vib = 0

#     # Calculate ΔHvap
#     # delta_H_vap = -U_liquid_per_mol + R * T - pV_term + C + E_vib
#     delta_H_vap = U_gas - U_liquid_per_mol + R * T - pV_term + C + E_vib

#     return delta_H_vap


# def calculate_delta_H_vap(T, U_liquid, V_liquid, kappa_T, alpha_p, pgm=False):
#     # Constants
#     N = 512  # Number of water molecules
#     p_ext = 1  # External pressure in bar
#     R = 0.001987  # Gas constant in kcal/(mol*K)

#     # Antoine equation constants for water (T in °C, p in mm Hg)
#     A = 8.07131
#     B = 1730.63
#     C = 233.426

#     # Convert temperature to Celsius for the Antoine equation
#     T_C = T - 273.15

#     # Calculate p_vap using the Antoine equation (mm Hg to bar conversion)
#     p_vap_mmHg = 10 ** (A - B / (T_C + C))
#     p_vap = p_vap_mmHg * 0.00133322  # Convert mm Hg to bar

#     if pgm:
#         U_gas = -976.4144
#     else:
#         U_gas = 0

#     # Calculate Cvib and Cni
#     Cvib = 0.002044 * T * 0.704656
#     Cni = -0.000569 * T + 0.133952

#     if pgm:
#         C = Cvib
#     else:
#         C = Cvib + Cni

#     # Convert V_liquid from Angstrom^3 to L/mol
#     V_liquid_L_mol = V_liquid * 1e-27 * 6.022e23 / 1000 / N

#     # Assume U_liquid is already in kcal for the whole system
#     # Convert to kcal/mol by dividing by N
#     U_liquid_per_mol = U_liquid / N

#     # Convert p*V term to kcal/mol
#     # 1 L·bar = 0.0241838 kcal/mol
#     pV_term = p_ext * V_liquid_L_mol * 0.0241838

#     # Calculate E_vib (placeholder function, assuming it's defined elsewhere)
#     E_vib = calculate_quantum_vibrational_heat_capacity(T)
#     print("E_vib:", E_vib)

#     # Calculate Cx correction
#     C_x = (
#         V_liquid_L_mol * (1 - (p_vap - p_ext) * kappa_T) - T * V_liquid_L_mol * alpha_p
#     )

#     # Calculate ΔHvap including Cx correction
#     delta_H_vap = U_gas - U_liquid_per_mol + R * T - pV_term + C + E_vib + C_x

#     return delta_H_vap


def calculate_delta_H_vap(T, U_liquid, V_liquid, kappa_T, alpha_p, pgm=False):
    # Constants
    N = 512  # Number of water molecules
    p_ext = 1  # External pressure in bar
    R = 0.001987  # Gas constant in kcal/(mol*K)

    # Antoine equation constants for water (T in °C, p in mm Hg)
    A = 8.07131
    B = 1730.63
    C = 233.426

    # Convert temperature to Celsius for the Antoine equation
    T_C = T - 273.15

    # Calculate p_vap using the Antoine equation (mm Hg to bar conversion)
    p_vap_mmHg = 10 ** (A - B / (T_C + C))
    p_vap = p_vap_mmHg * 0.00133322  # Convert mm Hg to bar

    if pgm:
        U_gas = -976.4144
    else:
        U_gas = 0

    # Calculate Cvib and Cni
    Cvib = 0.002044 * T * 0.704656
    Cni = -0.000569 * T + 0.133952

    if pgm:
        C = Cvib
    else:
        C = Cvib + Cni

    # Convert V_liquid from Angstrom^3 to L/mol
    V_liquid_L_mol = V_liquid * 1e-27 * 6.022e23 / 1000 / N

    # Convert p*V term to kcal/mol
    # 1 L·bar = 0.0241838 kcal/mol
    pV_term = p_ext * V_liquid_L_mol * 0.0241838

    # Calculate E_vib (placeholder function, assuming it's defined elsewhere)
    E_vib = calculate_quantum_vibrational_heat_capacity(T)
    print("E_vib:", E_vib)

    # Calculate Cx correction
    C_x = (
        V_liquid_L_mol * (1 - (p_vap - p_ext) * kappa_T) - T * V_liquid_L_mol * alpha_p
    )

    # Long-range Lennard-Jones correction (placeholder values for sigma and epsilon)
    r_cutoff = 9.0  # Angstrom, cutoff distance
    sigma = 3.1815599649844395  # Angstrom, for water (example value)
    epsilon = 0.1447267928698027  # kcal/mol, for water (example value)

    # Calculate number density (molecules per Angstrom^3)
    rho = N / V_liquid

    # Calculate the long-range correction for Lennard-Jones interactions
    lj_correction = (8 * 3.14159 * rho * epsilon * sigma**6) / (3 * r_cutoff**3) - (
        8 * 3.14159 * rho * epsilon * sigma**12
    ) / (9 * r_cutoff**9)

    print(lj_correction)

    # Convert to kcal/mol by dividing by N
    U_liquid_per_mol = U_liquid / N + lj_correction

    # Calculate ΔHvap including Cx correction and long-range correction
    delta_H_vap = U_gas - U_liquid_per_mol + R * T - pV_term + C + E_vib + C_x

    return delta_H_vap


def calculate_epol(mu, alpha_gas=1.45):
    mu_gas = 0.3853  # Experimental gas-phase dipole moment of water in e*Angstrom (1.85 Debye)

    # Calculate Epol
    epol = (mu - mu_gas) ** 2 / (2 * alpha_gas)

    # Convert from e^2/Angstrom to kcal/mol
    # 1 e^2/Angstrom = 332.06371 kcal/mol
    epol_kcal_mol = epol * 332.06371

    return epol_kcal_mol


def calculate_delta_H_vap_KB(T, rdf, r, density, U_config):
    """
    Calculate heat of vaporization using Kirkwood-Buff integral method.

    :param T: Temperature in Kelvin
    :param rdf: Radial distribution function g(r)
    :param r: Radial distance array corresponding to rdf
    :param density: Number density of the liquid in molecules/Å³
    :param U_config: Configurational energy per molecule in kcal/mol
    :return: Heat of vaporization in kcal/mol
    """
    R = 0.001987  # Gas constant in kcal/(mol·K)

    # Calculate Kirkwood-Buff integral
    KB_integral = simps(4 * np.pi * r**2 * (rdf - 1), r)

    # Calculate isothermal compressibility
    kT = 1 / (density * R * T * (1 + density * KB_integral))

    # Calculate internal pressure
    P_int = density * R * T * (1 - density * KB_integral / 3)

    # Calculate heat of vaporization
    delta_H_vap = R * T - U_config + P_int / density

    return delta_H_vap


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

    # EPtot extraction for the main case
    eptot_file_pattern = "MD/3_Prod/*.out"
    eptot_extractor = CpStepExtractor(proj_path, case, wat_num, eptot_file_pattern)
    step_data = eptot_extractor.run()

    temperatures = set(range(240, 381, 5))
    temps = [298] + list(temperatures)

    test_keys = ["MD_data-pGM-temp", "MD_data-pGM-temp-1", "MD_data-pGM-temp-3"]
    bench_test_keys = ["MD_data_Benchamrk-pgm_mpi"]

    multi_run_data = {}
    for temp in temps:
        multi_run_data[float(temp)] = {
            "eptot": [],
            "volume": [],
        }  # Initialize the lists

        for test_key in test_keys:
            # Get the data for the current test_key
            tmp_eptot = step_data[test_key][case][float(temp)]["eptot"]
            tmp_volume = step_data[test_key][case][float(temp)]["volume"]

            # Extend the lists to concatenate data from all test_keys
            multi_run_data[float(temp)]["eptot"].extend(tmp_eptot)
            multi_run_data[float(temp)]["volume"].extend(tmp_volume)

    property_data = {"vital-cosmos-120_340": {"Cp": {}, "kT": {}, "Hvap": {}}}
    for temp in temps:
        # Perform the calculation: eptot_data + volume_data * 0.0144
        eptot_data = multi_run_data[float(temp)]["eptot"]
        volume_data = multi_run_data[float(temp)]["volume"]
        calculated_data = [e + v * 0.0144 for e, v in zip(eptot_data, volume_data)]

        qm_correction = calculate_quantum_vibrational_heat_capacity(temp)

        C_P_c = calculate_corrected_heat_capacity(temp, qm_correction, calculated_data)
        print(f"Corrected Heat Capacity: {C_P_c} cal/(mol·K)")

        kappa_T = calculate_isothermal_compressibility(volume_data, temp)
        print(
            f"{temp}: Isothermal compressibility (κ_T) per mole = {kappa_T} [10^-6 bar^-1]"
        )

        # New ΔHvap calculation
        avg_U_liquid = sum(eptot_data) / len(eptot_data)
        avg_V_liquid = sum(volume_data) / len(volume_data)
        # delta_H_vap = calculate_delta_H_vap(temp, avg_U_liquid, avg_V_liquid, pgm=True)
        with open("avg_data.pkl", "rb") as f:
            avg_data = pickle.load(f)

        print(avg_data)
        if temp == 240:
            alpha_p = avg_data["vital-cosmos-120_340"]["alpha_p"][float(temp) + 5]
        elif temp == 380:
            alpha_p = avg_data["vital-cosmos-120_340"]["alpha_p"][float(temp) - 5]
        else:
            alpha_p = avg_data["vital-cosmos-120_340"]["alpha_p"][float(temp)]
        delta_H_vap = calculate_delta_H_vap(
            temp, avg_U_liquid, avg_V_liquid, kappa_T, alpha_p, pgm=True
        )
        print(f"Heat of Vaporization (ΔHvap): {delta_H_vap:.4f} kcal/mol")

        property_data["vital-cosmos-120_340"]["Cp"][float(temp)] = C_P_c
        property_data["vital-cosmos-120_340"]["kT"][float(temp)] = kappa_T
        property_data["vital-cosmos-120_340"]["Hvap"][float(temp)] = delta_H_vap

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

    bench_folder_pattern = "MD_data_Benchamrk-pgm_mpi"
    bench_proj_path = f"/home8/yxwu/pGM_water_model/{bench_folder_pattern}"

    bench_multi_run_data = {}
    for bench_case in bench_cases:
        if bench_case not in property_data:
            property_data[bench_case] = {
                "Cp": {},
                "kT": {},
                "Hvap": {},
            }  # Initialize for each case

        # EPtot extraction for the main case
        bench_eptot_extractor = CpStepExtractor(
            bench_proj_path, bench_case, wat_num, eptot_file_pattern
        )
        bench_step_data = bench_eptot_extractor.run()
        for temp in temps:
            bench_multi_run_data[float(temp)] = {
                "eptot": [],
                "volume": [],
            }  # Initialize the lists

            for bench_test_key in bench_test_keys:
                # Get the data for the current test_key
                tmp_eptot = bench_step_data[bench_test_key][bench_case][float(temp)][
                    "eptot"
                ]
                tmp_volume = bench_step_data[bench_test_key][bench_case][float(temp)][
                    "volume"
                ]

                # Extend the lists to concatenate data from all test_keys
                bench_multi_run_data[float(temp)]["eptot"].extend(tmp_eptot)
                bench_multi_run_data[float(temp)]["volume"].extend(tmp_volume)

                # Perform the calculation: eptot_data + volume_data * 0.0144
                eptot_data = bench_multi_run_data[float(temp)]["eptot"]
                volume_data = bench_multi_run_data[float(temp)]["volume"]
                calculated_data = [
                    e + v * 0.0144 for e, v in zip(eptot_data, volume_data)
                ]

                qm_correction = calculate_quantum_vibrational_heat_capacity(temp)

                C_P_c = calculate_corrected_heat_capacity(
                    temp, qm_correction, calculated_data
                )
                print(f"Corrected Heat Capacity: {C_P_c} cal/(mol·K)")

                kappa_T = calculate_isothermal_compressibility(volume_data, temp)
                print(
                    f"{temp}: Isothermal compressibility (κ_T) per mole = {kappa_T} [10^-6 bar^-1]"
                )

                # New ΔHvap calculation
                avg_U_liquid = sum(eptot_data) / len(eptot_data)
                avg_V_liquid = sum(volume_data) / len(volume_data)
                # delta_H_vap = calculate_delta_H_vap(
                #     temp, avg_U_liquid, avg_V_liquid, pgm=False
                # )
                if temp == 240:
                    pass
                    alpha_p = avg_data["vital-cosmos-120_340"]["alpha_p"][
                        float(temp) + 5
                    ]
                elif temp == 380:
                    alpha_p = avg_data["vital-cosmos-120_340"]["alpha_p"][
                        float(temp) - 5
                    ]
                else:
                    alpha_p = avg_data["vital-cosmos-120_340"]["alpha_p"][float(temp)]
                delta_H_vap = calculate_delta_H_vap(
                    temp, avg_U_liquid, avg_V_liquid, kappa_T, alpha_p, pgm=False
                )
                print(f"Heat of Vaporization (ΔHvap): {delta_H_vap:.4f} kcal/mol")
                print(f"Heat of Vaporization (ΔHvap): {delta_H_vap:.4f} kcal/mol")

                property_data[bench_case]["Cp"][float(temp)] = C_P_c
                property_data[bench_case]["kT"][float(temp)] = kappa_T
                property_data[bench_case]["Hvap"][float(temp)] = delta_H_vap

    # print(property_data)
    # print(volume_data)

    # with open("avg_data.pkl", "rb") as f:
    #     avg_data = pickle.load(f)

    # property_data = calculate_isothermal_compressibility(avg_data, property_data)

    expt_path_with_names = {
        "Cp": "datasets/expt/expt_temp_heat_capacity-TIP4PEw.txt",
        "kT": "datasets/expt/expt_temp_isothermal_compressibility-TIP4PEw.txt",
        "Hvap": "datasets/expt/expt_temp_heat_vap-TIP4PEw.txt",
    }

    # read experiment data
    for property in expt_path_with_names.keys():

        file_path = expt_path_with_names[property]
        new_data = read_expt_data(file_path, name=property)

        if "Expt" not in property_data:
            property_data["Expt"] = {}

        property_data["Expt"].update(new_data)

    print(property_data)

    with open("property_data.pkl", "wb") as f:
        pickle.dump(property_data, f)

    with open("all_data.pkl", "rb") as f:
        all_data = pickle.load(f)

    # data reorganization
    alpha_p = {}
    for model in all_data.keys():
        alpha_p[model] = {}
        if model in bench_cases:
            test_keys = bench_test_keys
            for test_key in test_keys:
                alpha_p[model][test_key] = {}
                density_dict = all_data[model]["density"][test_key][model]
                coefficients = calculate_thermal_expansion_coefficients(density_dict)
                alpha_p[model][test_key]["alpha_p"] = coefficients

    print("*" * 180)
