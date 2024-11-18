import MDAnalysis as mda
import numpy as np
import glob
import os
import matplotlib.pyplot as plt


def calculate_average_dipole_moment_per_molecule(universe, selection="resname WAT"):
    """
    Calculate the average dipole moment of individual water molecules from an MD simulation.

    Parameters:
        universe (MDAnalysis.Universe): The MDAnalysis Universe containing trajectory and topology.
        selection (str): Atom selection string for water molecules (default is "resname WAT").

    Returns:
        float: Average dipole moment (in Debye) of water molecules.
    """
    # Select water molecules
    water_molecules = universe.select_atoms(selection)

    # Check that the selection contains water molecules
    if water_molecules.n_atoms == 0:
        raise ValueError(
            "No atoms found in selection. Check the 'resname' or atom naming convention."
        )

    # Group water molecules by residues (each water molecule is a residue)
    water_residues = water_molecules.residues

    # Conversion factor from e·Å to Debye (1 e·Å = 4.80321 Debye)
    conversion_factor = 4.80321

    # Initialize list to store dipole moments for each frame
    dipole_moments = []

    # Loop through each frame of the trajectory
    for ts in universe.trajectory:
        frame_dipoles = []

        # Calculate dipole for each water molecule
        for residue in water_residues:
            charges = residue.atoms.charges
            positions = residue.atoms.positions
            water_com = residue.atoms.center_of_mass()

            # Calculate the dipole moment of this water molecule
            dipole = np.sum(charges[:, np.newaxis] * (positions - water_com), axis=0)

            # Convert to Debye and compute the magnitude of the dipole moment
            dipole_magnitude = np.linalg.norm(dipole) * conversion_factor
            frame_dipoles.append(dipole_magnitude)

        # Store the average dipole moment for this frame
        dipole_moments.append(np.mean(frame_dipoles))

    # Calculate the average dipole moment over all frames
    average_dipole_moment = np.mean(dipole_moments)
    return average_dipole_moment


def calculate_dipole_moment_per_frame(universe, selection="resname WAT"):
    """
    Calculate the dipole moment of water molecules for each frame from an MD simulation.

    Parameters:
        universe (MDAnalysis.Universe): The MDAnalysis Universe containing trajectory and topology.
        selection (str): Atom selection string for water molecules (default is "resname WAT").

    Returns:
        list: Dipole moment (in Debye) of water molecules for each frame.
    """
    water_molecules = universe.select_atoms(selection)

    if water_molecules.n_atoms == 0:
        raise ValueError(
            "No atoms found in selection. Check the 'resname' or atom naming convention."
        )

    water_residues = water_molecules.residues
    conversion_factor = 4.80321
    dipole_moments = []

    for ts in universe.trajectory:
        frame_dipole_sum = 0.0

        for residue in water_residues:
            charges = residue.atoms.charges
            positions = residue.atoms.positions
            water_com = residue.atoms.center_of_mass()
            dipole = np.sum(charges[:, np.newaxis] * (positions - water_com), axis=0)
            dipole_magnitude = np.linalg.norm(dipole) * conversion_factor
            frame_dipole_sum += dipole_magnitude

        dipole_moments.append(frame_dipole_sum / len(water_residues))

    return dipole_moments


def calculate_average_quadrupole_moment_per_molecule(universe, selection="resname WAT"):
    """
    Calculate the average quadrupole moment of individual water molecules from an MD simulation.

    Parameters:
        universe (MDAnalysis.Universe): The MDAnalysis Universe containing trajectory and topology.
        selection (str): Atom selection string for water molecules (default is "resname WAT").

    Returns:
        np.ndarray: Average quadrupole moment tensor (in Debye·Å^2) of water molecules.
    """
    # Select water molecules
    water_molecules = universe.select_atoms(selection)

    # Check that the selection contains water molecules
    if water_molecules.n_atoms == 0:
        raise ValueError(
            "No atoms found in selection. Check the 'resname' or atom naming convention."
        )

    # Group water molecules by residues (each water molecule is a residue)
    water_residues = water_molecules.residues

    # Conversion factor from e·Å^2 to Debye·Å^2 (1 e·Å^2 = 4.80321 Debye·Å)
    conversion_factor = 4.80321

    # Initialize list to store quadrupole moments for each frame
    quadrupole_moments = []

    # Loop through each frame of the trajectory
    for ts in universe.trajectory:
        frame_quadrupoles = []

        # Calculate quadrupole for each water molecule
        for residue in water_residues:
            charges = residue.atoms.charges
            positions = residue.atoms.positions
            water_com = residue.atoms.center_of_mass()

            # Calculate the quadrupole moment tensor of this water molecule
            q_tensor = np.zeros((3, 3))
            for charge, position in zip(charges, positions):
                r = position - water_com
                q_tensor += charge * (3 * np.outer(r, r) - np.eye(3) * np.dot(r, r))

            # Convert to Debye·Å^2
            q_tensor *= conversion_factor
            frame_quadrupoles.append(q_tensor)

        # Store the average quadrupole moment tensor for this frame
        quadrupole_moments.append(np.mean(frame_quadrupoles, axis=0))

    # Calculate the average quadrupole moment tensor over all frames
    average_quadrupole_moment = np.mean(quadrupole_moments, axis=0)
    return average_quadrupole_moment


if __name__ == "__main__":

    # Specify the paths to your topology and trajectory folder
    # topology_file = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-2/TIP4PEW/512/TIP4PEW-298/MD/prep_files/tip4pew_512wat.prmtop"
    # trajectory_folder = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-2/TIP4PEW/512/TIP4PEW-298/MD/3_Prod"
    # topology_file = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-2/TIP4PEW/512/TIP4PEW-245/MD/prep_files/tip4pew_512wat.prmtop"
    # trajectory_folder = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi-2/TIP4PEW/512/TIP4PEW-245/MD/3_Prod"
    topology_file = "/home8/yxwu/pGM_water_model/MD_data-pGM-temp-2/vital-cosmos-120_340/512/vital-cosmos-120_340-298/MD/prep_files/case_1_pgm_512wat.prmtop"
    trajectory_folder = "/home8/yxwu/pGM_water_model/MD_data-pGM-temp-2/vital-cosmos-120_340/512/vital-cosmos-120_340-298/MD/3_Prod"
    # topology_file = "/home8/yxwu/pGM_water_model/MD_data-pGM-temp-1/vital-cosmos-120_340/512/vital-cosmos-120_340-245/MD/prep_files/case_1_pgm_512wat.prmtop"
    # trajectory_folder = "/home8/yxwu/pGM_water_model/MD_data-pGM-temp-1/vital-cosmos-120_340/512/vital-cosmos-120_340-245/MD/3_Prod"
    # topology_file = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi/OPC/512/OPC-298/MD/prep_files/opc_512wat.prmtop"
    # trajectory_folder = "/home8/yxwu/pGM_water_model/MD_data_Benchamrk-pgm_mpi/OPC/512/OPC-298/MD/3_Prod"

    # Find all .nc trajectory files in the folder
    trajectory_files = sorted(glob.glob(os.path.join(trajectory_folder, "*.nc")))

    # Load the universe with the topology and all trajectory files
    universe = mda.Universe(topology_file, trajectory_files)

    # Calculate the average dipole moment for water in the liquid phase
    mu_liquid = calculate_average_dipole_moment_per_molecule(universe)
    print(f"Average dipole moment in the liquid phase: {mu_liquid:} Debye")

    # # Use this function to calculate and plot dipole moment per frame:
    # dipole_moments = calculate_dipole_moment_per_frame(universe)
    # plt.plot(dipole_moments)
    # plt.xlabel("Frame")
    # plt.ylabel("Dipole Moment (Debye)")
    # plt.title("Dipole Moment per Frame over Time")
    # plt.savefig("dipole_moment_plot.png")

    # # Calculate the quadrupole moment for water in the liquid phase
    # quad_universe = mda.Universe(topology_file, trajectory_files)
    # results = calculate_multipole_moments(quad_universe)
    # print(f"Average dipole_moment: {results['dipole_moment']} D")
    # print(f"Average Q0: {results['Q0']} e·Å²")
    # print(f"Average QT: {results['QT']} e·Å²")
    # print(f"Average omega0: {results['Omega0']} e·Å²")
    # print(f"Average omegaT: {results['OmegaT']} e·Å²")

    # Calculate the quadrupole moment for water in the liquid phase
    quad_universe = mda.Universe(topology_file, trajectory_files)
    results = calculate_average_quadrupole_moment_per_molecule(quad_universe)
    # print(f"Average Q0: {results['Q0']} e·Å²")
    # print(f"Average QT: {results['QT']} e·Å²")
    print(results)
    Q = np.array(results)
    # Calculate Q_0: isotropic part
    Q_0 = np.trace(Q) / 3

    # Calculate Q_T: traceless part
    Q_T = Q - Q_0 * np.eye(3)

    print(f"Q_0 (isotropic part): {Q_0}")
    print(f"Q_T (traceless part):\n{Q_T}")

    Q_T_norm = np.linalg.norm(Q_T)
    print(f"Norm of Q_T: {Q_T_norm}")
