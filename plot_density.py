import glob
import os

from pathlib import Path
from re import A
from matplotlib import pyplot as plt

from src.utils.file_parser import parse_step_info
from src.utils.keys_list import keys_list


def plot_info(ax, f_dict, x, y, title, step_num, kwargs):
    steps = list(range(1, step_num + 1))
    plot_kwargs = kwargs.get("plot", {})
    ax.plot(
        [f_dict[str(i)][x] for i in steps],
        [f_dict[str(i)][y] for i in steps],
        **plot_kwargs,
    )
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    # ax.set_xscale('symlog')
    ax.set_title(title)


def plot_file(fig, ax, f, f_type, x="TIME(PS)", y="Etot", title="", kwargs={}):
    f_dict = parse_step_info(f, keys_list[f_type])
    plot_info(ax, f_dict, x, y, title, len(f_dict), kwargs)
    return fig, ax


output_path = Path("/home/yxwu/pGM_water_model/outputs/fig")
bench_output_folder = output_path / "Benchmark"
bench_output_folder.mkdir(parents=True, exist_ok=True)

fig, ax = plt.subplots(figsize=(10, 4))

outfile_pattern = "tip3p_512wat.npt.pmemd.*.out"

# outfile_folder = Path(
#     "/home/yxwu/pGM_water_model/MD_data_Benchmark/TIP3P/512/TIP3P_512wat-240/MD/4_Prod"
# )

# wat_model = "TIP3P"
wat_model = "TIP5P"

# for wat_model in wat_models:
#     wat_model_path = f"/home/yxwu/pGM_water_model/MD_data_Benchmark/{wat_model}/512"
#     temp_paths = os.listdir(wat_model_path)
#     for temp_path in temp_paths:
#         outfile_folders = [
#             f"{wat_model_path}/{temp_path}/MD/3_Langevin_Prod",
#             f"{wat_model_path}/{temp_path}/MD/4_Prod",
#         ]

#         for outfile_folder in outfile_folders:
#             paths = glob.glob(f"{outfile_folder}/{outfile_pattern}")
#             for path in paths:
#                 fig, ax = plt.subplots(figsize=(10, 4))
#                 fig, ax = plot_file(
#                     fig,
#                     ax,
#                     path,
#                     "npt",
#                     x,
#                     y,
#                     # kwargs={"plot": {"label": "TIP3P_512wat-240"}},
#                 )

#             ax.set_title(f"{wat_model}_{y}")
#             ax.legend()
#             plt.show()

#             save_dir = Path(bench_output_folder / wat_model)
#             save_dir.mkdir(parents=True, exist_ok=True)
#             fig.savefig(str(save_dir / f"{temp_path}.pdf"))


temps = ["298"] + [str(t) for t in range(240, 381, 5)]

for temp in temps:
    fig, ax = plt.subplots(figsize=(10, 4))
    outfile_folders = [
        f"/home/yxwu/pGM_water_model/MD_data_Benchmark/{wat_model}/512/{wat_model}_512wat-{temp}/MD/3_Langevin_Prod",
        f"/home/yxwu/pGM_water_model/MD_data_Benchmark/{wat_model}/512/{wat_model}_512wat-{temp}/MD/4_Prod",
    ]

    x = "TIME(PS)"
    y = "Density"

    for outfile_folder in outfile_folders:
        paths = glob.glob(f"{outfile_folder}/{outfile_pattern}")
        for path in paths:
            fig, ax = plot_file(
                fig,
                ax,
                path,
                "npt",
                x,
                y,
                # kwargs={"plot": {"label": "TIP3P_512wat-240"}},
            )

    ax.set_title(f"{wat_model}_{y}")
    ax.legend()
    plt.show()
    save_dir = bench_output_folder / wat_model
    save_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(save_dir / f"{wat_model}-{temp}.pdf"))
