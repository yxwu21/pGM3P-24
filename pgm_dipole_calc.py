import parmed as pmd
import numpy as np
import re
import MDAnalysis as mda


class QuadrupoleCalculator:
    def __init__(self, parmtop_file, dipole_file):
        """
        Initialize the QuadrupoleCalculator with parmtop and dipole data files.

        Args:
            parmtop_file (str): Path to the parmtop file.
            dipole_file (str): Path to the dipole data text file.
        """
        self.parmtop_file = parmtop_file
        self.dipole_file = dipole_file
        self.masses = None
        self.charges = None
        self.coordinates = None
        self.dipole_data = None
        self.quadrupole_data = {}

        # Load dipole data immediately upon initialization
        self.load_dipole_data()

    def load_parmtop_data(self, coordinate_file):
        """Load atomic data from the parmtop file and coordinates from a separate file."""
        parmtop = pmd.load_file(self.parmtop_file)
        self.masses = np.array([atom.mass for atom in parmtop.atoms])
        self.charges = np.array([atom.charge for atom in parmtop.atoms])

        # Load coordinates from the specified coordinate file (e.g., .rst7, .inpcrd)
        coordinates = pmd.load_file(coordinate_file)
        self.coordinates = coordinates.coordinates.reshape(
            -1, 3
        )  # Reshape to get (N_atoms, 3)
        print(f"Loaded coordinates from {coordinate_file}.")
        print("Loaded parmtop data successfully.")

    def load_dipole_data(self):
        """Load dipole moment data from the provided text file using the verified function."""
        self.dipole_data = self.extract_dipole_data(self.dipole_file)
        print(f"Loaded dipole data successfully: {self.dipole_data}")

    @staticmethod
    def extract_dipole_data(file_path):
        """Extracts dipole moments from a text file and computes the averages."""
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

    def calculate_center_of_mass(self, positions):
        """Calculate the center of mass."""
        total_mass = np.sum(self.masses)
        com = np.sum(positions * self.masses[:, np.newaxis], axis=0) / total_mass
        return com

    def calculate_moment_of_inertia(self, positions):
        """Calculate the moment of inertia tensor."""
        inertia_tensor = np.zeros((3, 3))
        for i, mass in enumerate(self.masses):
            r = positions[i]
            inertia_tensor[0, 0] += mass * (r[1] ** 2 + r[2] ** 2)
            inertia_tensor[1, 1] += mass * (r[0] ** 2 + r[2] ** 2)
            inertia_tensor[2, 2] += mass * (r[0] ** 2 + r[1] ** 2)
            inertia_tensor[0, 1] -= mass * r[0] * r[1]
            inertia_tensor[0, 2] -= mass * r[0] * r[2]
            inertia_tensor[1, 2] -= mass * r[1] * r[2]

        inertia_tensor[1, 0] = inertia_tensor[0, 1]
        inertia_tensor[2, 0] = inertia_tensor[0, 2]
        inertia_tensor[2, 1] = inertia_tensor[1, 2]
        return inertia_tensor

    def reorient_dipole(self, positions):
        """Reorient dipole moments along principal axes."""
        com = self.calculate_center_of_mass(positions)
        positions -= com
        inertia_tensor = self.calculate_moment_of_inertia(positions)
        eigenvalues, eigenvectors = np.linalg.eigh(inertia_tensor)

        # Use averaged dipole vector
        dipole_vector = np.array(
            [self.dipole_data["x"], self.dipole_data["y"], self.dipole_data["z"]]
        )
        dipole_reoriented = np.dot(eigenvectors.T, dipole_vector)
        return dipole_reoriented, eigenvectors

    def calculate_quadrupole(self, positions, dipole, eigenvectors):
        """Calculate the quadrupole tensor."""
        q = np.zeros((3, 3))
        for i in range(len(self.charges)):
            r = np.dot(eigenvectors.T, positions[i])
            q[0, 0] += self.charges[i] * (r[0] ** 2)
            q[1, 1] += self.charges[i] * (r[1] ** 2)
            q[2, 2] += self.charges[i] * (r[2] ** 2)
            q[0, 1] += self.charges[i] * r[0] * r[1]
            q[0, 2] += self.charges[i] * r[0] * r[2]
            q[1, 2] += self.charges[i] * r[1] * r[2]

        q[0, 0] += 2 * dipole[0] * positions[:, 0].sum()
        q[1, 1] += 2 * dipole[1] * positions[:, 1].sum()
        q[2, 2] += 2 * dipole[2] * positions[:, 2].sum()
        q[0, 1] += dipole[0] * positions[:, 1].sum() + dipole[1] * positions[:, 0].sum()
        q[0, 2] += dipole[0] * positions[:, 2].sum() + dipole[2] * positions[:, 0].sum()
        q[1, 2] += dipole[1] * positions[:, 2].sum() + dipole[2] * positions[:, 1].sum()

        q[1, 0] = q[0, 1]
        q[2, 0] = q[0, 2]
        q[2, 1] = q[1, 2]
        return q

    def compute_quadrupole(self):
        """Compute the quadrupole tensor using averaged dipole values."""
        dipole_reoriented, eigenvectors = self.reorient_dipole(self.coordinates.copy())
        quadrupole_tensor = self.calculate_quadrupole(
            self.coordinates.copy(), dipole_reoriented, eigenvectors
        )
        self.quadrupole_data = quadrupole_tensor
        print("Computed quadrupole tensor using averaged dipole data.")

    def get_quadrupole_data(self):
        """Returns the calculated quadrupole data."""
        return self.quadrupole_data

    def calculate_qt(quadrupole_tensor):
        """
        Calculate the trace (QT) of the quadrupole tensor with optional unit conversion.

        Args:
            quadrupole_tensor (np.array): A 3x3 quadrupole tensor.
            conversion_factor (float): Factor to convert units, e.g., from Bohr^2 to Å^2.

        Returns:
            float: The trace of the quadrupole tensor after unit conversion.
        """
        trace = np.trace(quadrupole_tensor)
        return trace


class TrajQuadrupoleCalculator(QuadrupoleCalculator):
    def __init__(self, parmtop_file, trajectory_file, dipole_file):
        """
        Initialize the TrajQuadrupoleCalculator with topology and trajectory data files.

        Args:
            parmtop_file (str): Path to the topology file (e.g., .prmtop).
            trajectory_file (str): Path to the trajectory file (e.g., .dcd, .nc).
        """
        self.trajectory_file = trajectory_file
        super().__init__(parmtop_file, dipole_file)

    def load_dipole_data(self):
        """Load dipole moment data from the provided text file using the verified function."""
        self.dipole_data = self.extract_dipole_data()
        print(f"Loaded dipole data successfully: {self.dipole_data}")

    def extract_dipole_data(self):
        """
        Calculate the dipole moment from a molecular dynamics trajectory.

        Returns:
            dict: Averaged dipole moment components (x, y, z, total).
        """
        u = mda.Universe(self.parmtop_file, self.trajectory_file)
        charges = u.atoms.charges  # Get charges from the topology
        dipole_moments = []

        for ts in u.trajectory:
            positions = u.atoms.positions  # Positions in Å
            dipole = np.sum(
                charges[:, np.newaxis] * positions, axis=0
            )  # Calculate dipole moment vector (e·Å)
            dipole_moments.append(dipole)

        # Convert list to numpy array for easier manipulation
        dipole_moments = np.array(dipole_moments)

        # Calculate averages
        dipole_x_avg = np.mean(dipole_moments[:, 0])
        dipole_y_avg = np.mean(dipole_moments[:, 1])
        dipole_z_avg = np.mean(dipole_moments[:, 2])
        dipole_total_avg = np.sqrt(dipole_x_avg**2 + dipole_y_avg**2 + dipole_z_avg**2)

        dipole_data = {
            "x": dipole_x_avg,
            "y": dipole_y_avg,
            "z": dipole_z_avg,
            "t": dipole_total_avg,
        }
        print("Calculated dipole moment from trajectory data.")
        return dipole_data

    def compute_quadrupole_from_trajectory(self):
        """
        Compute the quadrupole tensor using dipole data obtained from the trajectory.
        """
        self.calculate_dipole_from_trajectory()
        self.compute_quadrupole()


# Main script
if __name__ == "__main__":
    parmtop_file = "/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-298/nvt-npt-dipole/vital-cosmos-120_340/MD/prep_files/case_1_pgm_512wat.prmtop"
    dipole_file = "/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-298/nvt-npt-dipole/vital-cosmos-120_340/MD/dipole_1k-200ps_output_npt/fort.200"
    coordinate_file = "/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-298/nvt-npt-dipole/vital-cosmos-120_340/MD/dipole_1k-200ps_output_npt/case_1_pgm_512wat.npt.sander.rst"
    trajectory_file = "/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-298/nvt-npt-dipole/vital-cosmos-120_340/MD/dipole_1k-200ps_output_npt/case_1_pgm_512wat.npt.sander.nc"

    # Initialize the calculator
    # calculator = QuadrupoleCalculator(parmtop_file, dipole_file)
    calculator = TrajQuadrupoleCalculator(parmtop_file, trajectory_file, dipole_file)

    # Load parmtop data with the specified coordinate file
    calculator.load_parmtop_data(coordinate_file=coordinate_file)

    # Compute the quadrupole tensor
    calculator.compute_quadrupole()

    # Retrieve and print the quadrupole data
    quadrupole_data = calculator.get_quadrupole_data()
    print("Quadrupole Tensor (e·Å²):\n", quadrupole_data)

    qt_value = TrajQuadrupoleCalculator.calculate_qt(quadrupole_data)
    print(f"Q_T (converted to Å²): {qt_value}")
