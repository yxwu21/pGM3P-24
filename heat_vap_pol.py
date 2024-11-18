import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import rdf as mda_rdf
from scipy.integrate import simpson as simps
import glob


def extract_parameters(trajectories, topology, n_frames_per_traj=1000):
    """
    Extract parameters needed for Kirkwood-Buff heat of vaporization calculation from multiple trajectories.

    :param trajectories: List of paths to trajectory files
    :param topology: Path to topology file
    :param n_frames_per_traj: Number of frames to analyze per trajectory
    :return: rdf, r, density, U_config
    """
    all_volumes = []
    all_energies = []
    rdf_analyzers = []

    for trajectory in trajectories:
        u = mda.Universe(topology, trajectory)
        oxygen = u.select_atoms("name O")

        # RDF calculation
        rdf_analyzer = mda_rdf.InterRDF(oxygen, oxygen, nbins=500, range=(0, 15))
        rdf_analyzer.run()
        rdf_analyzers.append(rdf_analyzer)

        # Volume and energy collection
        for ts in u.trajectory[:n_frames_per_traj]:
            all_volumes.append(ts.volume)
            all_energies.append(ts.data.get("potential_energy", 0))

    # Combine RDF data
    combined_rdf = np.mean([analyzer.results.rdf for analyzer in rdf_analyzers], axis=0)
    r = rdf_analyzers[0].results.bins  # Assuming all have the same bins

    # Calculate average density
    avg_volume = np.mean(all_volumes)
    density = len(oxygen) / avg_volume  # molecules/Å³

    # Calculate average configurational energy per molecule
    U_config = np.mean(all_energies) / len(oxygen)  # kcal/mol

    return combined_rdf, r, density, U_config


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
    KB_integral = simps(y=4 * np.pi * r**2 * (rdf - 1), x=r)

    # Calculate isothermal compressibility
    kT = 1 / (density * R * T * (1 + density * KB_integral))

    # Calculate internal pressure
    P_int = density * R * T * (1 - density * KB_integral / 3)

    # Calculate heat of vaporization
    delta_H_vap = R * T - U_config + P_int / density

    return delta_H_vap


if __name__ == "__main__":

    # path = "/home8/yxwu/pGM_water_model/MD_data-pGM-temp-1/vital-cosmos-120_340/512/vital-cosmos-120_340-298/MD/3_Prod/case_1_pgm_512wat.npt.sander.*.nc"
    # topology = "/home8/yxwu/pGM_water_model/MD_data-pGM-temp-1/vital-cosmos-120_340/512/vital-cosmos-120_340-298/MD/prep_files/case_1_pgm_512wat.prmtop"

    # path = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi/TIP3P/512/TIP3P-298/MD/3_Prod/tip3p_512wat.npt.pmemd.*.nc"
    # topology = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi/TIP3P/512/TIP3P-298/MD/prep_files/tip3p_512wat.prmtop"

    # path = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/TIP4PEW/512/TIP4PEW-298/MD/3_Prod/tip4pew_512wat.npt.pmemd.*.nc"
    # topology = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/TIP4PEW/512/TIP4PEW-298/MD/prep_files/tip4pew_512wat.prmtop"

    # path = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/SPCE/512/SPCE-298/MD/3_Prod/spce_512wat.npt.pmemd.*.nc"
    # topology = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/SPCE/512/SPCE-298/MD/prep_files/spce_512wat.prmtop"

    # path = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/OPC/512/OPC-298/MD/3_Prod/opc_512wat.npt.pmemd.*.nc"
    # topology = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/OPC/512/OPC-298/MD/prep_files/opc_512wat.prmtop"

    # path = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/OPC3/512/OPC3-298/MD/3_Prod/opc3_512wat.npt.pmemd.*.nc"
    # topology = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/OPC3/512/OPC3-298/MD/prep_files/opc3_512wat.prmtop"

    path = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/TIP5P/512/TIP5P-298/MD/3_Prod/tip5p_512wat.npt.pmemd.*.nc"
    topology = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-1/TIP5P/512/TIP5P-298/MD/prep_files/tip5p_512wat.prmtop"

    trajectories = glob.glob(path)

    T = 298

    rdf, r, density, U_config = extract_parameters(trajectories, topology)

    print(f"Density: {density:.6f} molecules/Å³")
    print(f"Configurational Energy: {U_config:.4f} kcal/mol")

    delta_H_vap = calculate_delta_H_vap_KB(T, rdf, r, density, U_config)
    print(f"Heat of Vaporization: {delta_H_vap:.4f} kcal/mol")
