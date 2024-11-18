from data_extractor import read_expt_data
from data_plot import DataPlotter


if __name__ == "__main__":
    b2_data_to_paths = {
        "Expt": "datasets/B2/expt_temp_b2.txt",
        "TIP3P": "datasets/B2/tip3p_temp_b2.txt",
        "POL3-LT": "datasets/B2/pol3-lt_temp_b2.txt",
        "OPC": "datasets/B2/opc_temp_b2.txt",
        "OPC3": "datasets/B2/opc3_temp_b2.txt",
        "SPCE": "datasets/B2/spce_temp_b2.txt",
        "TIP4PEW": "datasets/B2/tip4pew_temp_b2.txt",
        "vital-cosmos-120_340": "datasets/B2/vital_cosmos-120-340_temp_b2.txt",
    }

    b2_data = {}
    for model in b2_data_to_paths.keys():

        file_path = b2_data_to_paths[model]
        new_data = read_expt_data(file_path, name="b2")

        b2_data[model] = new_data

    labels_markers = {
        "TIP3P": ("TIP3P", "o", "green", "--", 1.5, 8, True),
        # "TIP4P": ("TIP4P", "^", "red", "--", 1.5, 8, True),
        # "TIP5P": ("TIP5P", "s", "purple", "--", 1.5, 8, True),
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

    plotter = DataPlotter(b2_data, save_dir="plots/without_tip5p/b2", set_xlim=False)
    zoom_plotter = DataPlotter(
        b2_data, save_dir="plots/without_tip5p/b2", set_xlim=False, b2_zoom_in=True
    )

    plotter.set_custom_ylabel("b2", r"$B_2 \, [\mathrm{cm}^3/\mathrm{mol}]$")
    plotter.plot_cases("b2", labels_markers)
    zoom_plotter.plot_cases("b2", labels_markers)
