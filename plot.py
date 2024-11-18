import os
import numpy as np

from src.eval import peak_metrics
from src.yxwu_lib import read_analysis_data
from src.result_analysis_func import (
    plot_rdf_single,
    analyze_diffusion_constants,
    find_matching_paths,
    calculate_averages,
)
from pathlib import Path
from matplotlib import pyplot as plt
from tqdm import tqdm


def eval_rdf_result(
    identifier: str,
    rdf_type: str,
    rdf_analysis_folder: str,
    experimental_rdf_analysis_folder: str,
):
    # calculate rdf score
    data_file_pattern = f"{rdf_analysis_folder}/*.out"
    md_rdf, _ = read_analysis_data(rdf_analysis_folder, data_file_pattern, identifier)
    exp_rdf, _ = read_analysis_data(
        experimental_rdf_analysis_folder, data_file_pattern, identifier
    )
    rdf_score = peak_metrics(np.array(md_rdf), np.array(exp_rdf))

    data = {"rdf": md_rdf, "exp_rdf": exp_rdf}

    if rdf_type == "rdf_goo":
        soper_expt_analysis_folder = str(
            Path(experimental_rdf_analysis_folder).parents[2]
            / "expt_Soper_2013"
            / "analysis"
            / "rdf_goo"
        )
        exp_rdf_soper, _ = read_analysis_data(
            soper_expt_analysis_folder, data_file_pattern, identifier
        )
        data["exp_rdf_soper"] = np.array(exp_rdf_soper)
    return rdf_score, data


def process_rdf_plots(analysis_origin, rdf_types, identifiers, expt_rdf_folders):
    plot_lists = os.listdir(analysis_origin)
    for rdf_type in tqdm(rdf_types):
        for rdf_name in plot_lists:
            rdf_analysis_folder = f"{analysis_origin}/{rdf_name}/analysis/{rdf_type}"
            goo_score, goo_data = eval_rdf_result(
                identifiers[rdf_type],
                rdf_type,
                rdf_analysis_folder,
                expt_rdf_folders[rdf_type],
            )

            fig_folder = Path(rdf_analysis_folder).parents[1] / "figs"
            fig_folder.mkdir(exist_ok=True)

            plot_rdf_single(goo_data, rdf_type, rdf_name, fig_folder)


def plot_rdf(data: dict, rdf_type: str, rdf_name: str, ax):
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

    ax.set_xlabel("#Distance (Ang)")
    ax.set_ylabel(ylabels.get(rdf_type, "Unknown Type"))
    if rdf_type in ylimits:
        ax.set_ylim(ylimits[rdf_type])
    ax.legend()


def plot_rdfs(
    analysis_origin, rdf_types, identifiers, expt_rdf_folders, results, average_data
):
    plot_lists = os.listdir(analysis_origin)
    for rdf_name in tqdm(plot_lists):
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(8, 18))
        for i, rdf_type in enumerate(rdf_types):
            rdf_analysis_folder = f"{analysis_origin}/{rdf_name}/analysis/{rdf_type}"
            goo_score, goo_data = eval_rdf_result(
                identifiers[rdf_type],
                rdf_type,
                rdf_analysis_folder,
                expt_rdf_folders[rdf_type],
            )
            plot_rdf(goo_data, rdf_type, rdf_name, axes[i])

        add_figure_text(fig, rdf_name, results, average_data)
        save_figure(fig, rdf_name, analysis_origin)


def add_figure_text(fig, rdf_name, results, average_data):
    if rdf_name in results and "Diffusion Constant" in results[rdf_name]:
        diffusion_constant = results[rdf_name]["Diffusion Constant"]["r"]
        avg_density = average_data.get(rdf_name, {}).get("Density", "N/A")
        avg_etot = average_data.get(rdf_name, {}).get("Etot", "N/A")
        label_text = f"Diffusion constant: {diffusion_constant:.3f} \nDensity: {avg_density:.3f} \nΔHvap (Etot): {avg_etot:.3f}"
        fig.text(
            0.85,
            0.83,
            label_text,
            ha="right",
            va="top",
            fontsize=12,
            transform=fig.transFigure,
        )


def save_figure(fig, rdf_name, analysis_origin):
    fig_folder = Path(analysis_origin) / rdf_name / "figs"
    fig_folder.mkdir(exist_ok=True)
    output_path = fig_folder / f"{rdf_name}_rdf_total.pdf"
    fig.savefig(output_path)
    plt.close(fig)


def main():
    identifiers = {
        "rdf_goo": "#Distance_(Ang)     @O_=>_@O",
        "rdf_goh": "#Distance_(Ang)    @O_=>_@H*",
        "rdf_ghh": "#Distance_(Ang)   @H*_=>_@H*",
    }
    expt_rdf_folders = {
        "rdf_goo": "datasets/pgm_512wat_main/expt_Skinner_2013/analysis/rdf_goo",
        "rdf_goh": "datasets/pgm_512wat_main/expt_Soper_2013/analysis/rdf_goh",
        "rdf_ghh": "datasets/pgm_512wat_main/expt_Soper_2013/analysis/rdf_ghh",
    }
    analysis_origin = "/home/yxwu/pGM_water_model/MD_data/langevin_Berendsen-pmemd_backup/nvt-npt-dipole"
    rdf_types = ["rdf_goo", "rdf_goh", "rdf_ghh"]

    analysis_folder = "diffusion_-11_-1"
    results = analyze_diffusion_constants(analysis_origin, analysis_folder)
    average_data = calculate_averages(
        find_matching_paths(analysis_origin, "1c-200ps_10iter_output_npt"),
        "1c-200ps_10iter_output_npt",
        ["Density", "Etot"],
    )

    plot_rdfs(
        analysis_origin, rdf_types, identifiers, expt_rdf_folders, results, average_data
    )


if __name__ == "__main__":
    main()
