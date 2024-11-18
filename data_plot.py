import pickle
import matplotlib.pyplot as plt
import numpy as np
import os

from temp_result_analysis import extract_dipole_data
from pathlib import Path
from scipy.optimize import curve_fit


class DataPlotter:
    def __init__(
        self, avg_data, save_dir="plots/with_tip5p", set_xlim=True, b2_zoom_in=False
    ):
        self.avg_data = avg_data
        self.save_dir = save_dir
        self.set_xlim = set_xlim
        self.b2_zoom_in = b2_zoom_in
        os.makedirs(self.save_dir, exist_ok=True)

        self.custom_ylabels = {}

        self.title_size = 16
        self.axis_label_size = 14
        self.tick_label_size = 12
        self.legend_size = 10

    def set_custom_ylabel(self, property_key, ylabel):
        self.custom_ylabels[property_key] = ylabel

    def set_font_sizes(self, title=16, axis_label=14, tick_label=12, legend=10):
        self.title_size = title
        self.axis_label_size = axis_label
        self.tick_label_size = tick_label
        self.legend_size = legend

    def plot_cases(self, property_key, labels_markers):
        if self.b2_zoom_in:
            plt.figure(figsize=(4, 3))
        else:
            plt.figure(figsize=(8, 6))

        # Plot each case in the sequence of labels_markers
        for case, (
            label,
            marker,
            color,
            linestyle,
            linewidth,
            markersize,
            hollow_marker,
        ) in labels_markers.items():
            if case in self.avg_data and property_key in self.avg_data[case]:
                data = self.avg_data[case][property_key]

                # Ensure temperatures are sorted
                temperatures = sorted(
                    [float(temp) for temp in data.keys() if float(temp) != 298]
                )
                values = [data[temp] for temp in temperatures]

                # Plot data with hollow marker if specified
                markerfacecolor = "none" if hollow_marker else color
                plt.plot(
                    temperatures,
                    values,
                    label=label,
                    marker=marker,
                    color=color,
                    linestyle=linestyle,
                    linewidth=linewidth,
                    markersize=markersize,
                    markerfacecolor=markerfacecolor,  # Set hollow marker
                )

        if self.b2_zoom_in:
            plt.xlim(250, 500)
            plt.ylim(-6000, 200)
            file_name = f"{property_key}_zoom_in.pdf"
        else:
            # Customize plot appearance
            plt.xlabel("Temperature [K]", fontsize=self.axis_label_size)
            if property_key in self.custom_ylabels:
                plt.ylabel(
                    self.custom_ylabels[property_key], fontsize=self.axis_label_size
                )
            else:
                plt.ylabel(
                    f"{property_key.capitalize()}", fontsize=self.axis_label_size
                )

            # plt.title(f"{property_key.capitalize()} vs Temperature")
            plt.legend()
            plt.tight_layout()
            if self.set_xlim:
                plt.xlim(230, 390)

            # Save the plot
            file_name = f"{property_key}_comparison.pdf"
        plt.savefig(os.path.join(self.save_dir, file_name))
        plt.close()


def polynomial_fit(T, *coeffs):
    # General polynomial function for fitting
    return sum(c * T**i for i, c in enumerate(coeffs))


def calculate_tmd_for_all_models(avg_data, order=5):

    # Extract temperature list and density values
    temperature_list = list(avg_data[model]["density"].keys())
    density_values = list(avg_data[model]["density"].values())

    # Convert to numpy arrays
    temperatures = np.array(temperature_list)
    densities = np.array(density_values)

    # Fit the data using the polynomial function of the specified order
    initial_guess = np.ones(order + 1)  # Initial guess for the polynomial coefficients
    params, _ = curve_fit(polynomial_fit, temperatures, densities, p0=initial_guess)

    # Find the derivative to calculate the TMD (solving for the root where slope = 0)
    derivative_coeffs = [i * params[i] for i in range(1, len(params))]
    roots = np.roots(derivative_coeffs[::-1])  # Solve the derivative for roots

    # Filter real roots and select the TMD within the temperature range
    TMD_candidates = roots[np.isreal(roots)].real
    TMD = TMD_candidates[
        (TMD_candidates >= min(temperatures)) & (TMD_candidates <= max(temperatures))
    ]

    if len(TMD) == 0:
        TMD = None  # No valid TMD found
    else:
        TMD = TMD[0]  # Take the first valid TMD root

    # Generate the fit line for the plot
    T_fit = np.linspace(min(temperatures), max(temperatures), 100)
    D_fit = polynomial_fit(T_fit, *params)

    # Plot the data and the polynomial fit for each model
    plt.figure(figsize=(8, 6))
    plt.plot(temperatures, densities, "o", label="Density Data", color="blue")
    plt.plot(T_fit, D_fit, "-", label=f"Polynomial Fit (Order {order})", color="red")
    if TMD is not None:
        plt.axvline(x=TMD, color="green", linestyle="--", label=f"TMD = {TMD:.3f} K")
    plt.xlabel("Temperature T [K]")
    plt.ylabel("Density [g/cm³]")
    plt.title(f"Temperature vs. Density for {model}")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    # Create directory and save the figure
    folder = Path("plots/TMD")
    folder.mkdir(parents=True, exist_ok=True)
    plt.savefig(f"{folder}/{model}_order_{order}.pdf")

    # Show plot (optional, remove if running non-interactively)
    plt.show()

    # Return the TMD for this model
    return TMD


# Example usage
if __name__ == "__main__":

    # Load the all_data and avg_data dictionaries
    with open("all_data.pkl", "rb") as f:
        all_data = pickle.load(f)

    # print(all_data["vital-cosmos-120_340"]["eptot"])
    with open("avg_data.pkl", "rb") as f:
        avg_data = pickle.load(f)

    with open("property_data.pkl", "rb") as f:
        property_data = pickle.load(f)

    with open("std_dev_data.pkl", "rb") as f:
        std_dev_data = pickle.load(f)

    labels_markers = {
        "TIP3P": ("TIP3P", "o", "green", "--", 1.5, 8, True),
        # "TIP4P": ("TIP4P", "^", "red", "--", 1.5, 8, True),
        "TIP5P": ("TIP5P", "s", "purple", "--", 1.5, 8, True),
        "OPC": ("OPC", "<", "y", "--", 1.5, 8, True),
        "OPC3": ("OPC3", ">", "darkorange", "--", 1.5, 8, True),
        "SPCE": ("SPC/E", "8", "deeppink", "--", 1.5, 8, True),  # Hollow marker
        "TIP4PEW": ("TIP4P-Ew", "^", "red", "--", 1.5, 8, True),
        "Expt": ("Expt.", None, "black", "-", 1.5, 12, True),
        "vital-cosmos-120_340": (
            "pGM3P-24",
            "D",
            "blue",
            "-",
            1.5,
            8,
            True,
        ),
    }

    # Initialize DataPlotter with avg_data
    plotter = DataPlotter(avg_data)
    # print(avg_data)

    # Plot all cases for a given property (e.g., "density")
    plotter.set_custom_ylabel("density", r"$\rho$ $[g/cm^3]$")
    plotter.set_custom_ylabel("diffusion", r"$D$ $[10^{-5} \, \text{cm}^2/\text{s}]$")
    plotter.set_custom_ylabel("dielectric", r"Static Dielectric $\varepsilon_0$")
    plotter.set_custom_ylabel(
        "alpha_p",
        r"$\alpha_p \ [10^{-4} \ \mathrm{K}^{-1}]$",
    )
    plotter.plot_cases("density", labels_markers)
    plotter.plot_cases("diffusion", labels_markers)
    plotter.plot_cases("dielectric", labels_markers)
    plotter.plot_cases("alpha_p", labels_markers)

    print(property_data)
    property_plotter = DataPlotter(property_data)
    property_plotter.set_custom_ylabel("Cp", r"$C_p$ [cal/(K·mol)]")
    property_plotter.set_custom_ylabel("kT", r"10$^5$$\kappa$$_T$ [bar$^{-1}$]")
    property_plotter.set_custom_ylabel(
        "Hvap", r"$\Delta H_{\mathrm{vap}} \ [kcal/mol]$"
    )
    property_plotter.plot_cases("Cp", labels_markers)
    property_plotter.plot_cases("kT", labels_markers)
    property_plotter.plot_cases("Hvap", labels_markers)

    # printing properties
    print("Cp", property_data["vital-cosmos-120_340"]["Cp"][298])
    print("kT", property_data["vital-cosmos-120_340"]["kT"][298])

    temperature_list = list(avg_data["vital-cosmos-120_340"]["density"].keys())
    density_values = list(avg_data["vital-cosmos-120_340"]["density"].values())
    max_index = np.argmax(density_values)
    TMD = temperature_list[max_index]
    print("TMD:", TMD)

    print(avg_data.keys())

    print("diffusion:")
    for model in avg_data.keys():
        if model != "Expt":
            print(model, avg_data[model]["diffusion"][298])
        else:
            print(model, avg_data[model]["diffusion"])

    print("-" * 180)
    print("density:")
    for model in avg_data.keys():
        if model != "Expt":
            print(model, avg_data[model]["density"][298])
        else:
            print(model, avg_data[model]["density"])

    print("-" * 180)
    print("dielectric:")
    for model in avg_data.keys():
        print(model, avg_data[model]["dielectric"])

    print("-" * 180)
    folders = ["1d", "1e", "1f", "1g", "1h", "1i", "1j", "1k"]
    avg_dipole = []
    for folder in folders:
        dipole_data = extract_dipole_data(
            f"/home8/yxwu/pGM_water_model/MD_data-0/langevin_Berendsen-298/nvt-npt-dipole/vital-cosmos-120_340/MD/dipole_{folder}-200ps_output_npt/fort.200"
        )
        print(folder, ":", dipole_data)
        avg_dipole.append(dipole_data["t"])

    print(avg_dipole)
    print("average dipole:", np.average(avg_dipole) * 4.8)
    std_dev_dipole = np.std(avg_dipole)
    print("Standard deviation of dipole values:", std_dev_dipole * 4.8)

    print("-" * 180)
    print("Cp:")
    for model in property_data.keys():
        if model != "Expt":
            print(model, property_data[model]["Cp"][298])
        else:
            print(model, property_data[model]["Cp"])

    print("-" * 180)
    print("kT:")
    for model in property_data.keys():
        if model != "Expt":
            print(model, property_data[model]["kT"][298])
        else:
            print(model, property_data[model]["kT"])

    print("-" * 180)
    print("Hvap:")
    for model in property_data.keys():
        if model != "Expt":
            print(model, property_data[model]["Hvap"][298])
        else:
            print(model, property_data[model]["Hvap"])

    print("-" * 180)
    print("alpha_p:")
    for model in avg_data.keys():
        if model != "Expt":
            print(
                model,
                "-",
                "average:",
                avg_data[model]["alpha_p"][298.0],
                "|",
                "std_dev:",
                std_dev_data[model]["alpha_p"][298.0],
            )

    print("-" * 180)
    print("dipole:")
    for model in avg_data.keys():
        if model != "Expt":
            print(
                model,
                "-",
                "average:",
                avg_data[model]["dipole"][298.0],
                "|",
                "std_dev:",
                std_dev_data[model]["dipole"][298.0],
            )

    print("-" * 180)
    print("TMD:")
    for model in avg_data.keys():
        TMD = calculate_tmd_for_all_models(avg_data)
        print(model, "-", "average:", TMD)

    print("*" * 180)
    print("TIP5P:")
    print(avg_data["TIP5P"]["density"])
    print("-" * 180)
    print("Expt:")
    print(avg_data["Expt"]["density"])
