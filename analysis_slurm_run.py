import os

from pathlib import Path
from tqdm import tqdm
from mdtool import yxwu_lib
from rst_slurm_run import find_matching_paths, extract_subpath_after_key

# analysis_type = "rdf_ghh"

analysis_types = [
    # "rdf_ghh",
    # "rdf_goh",
    # "rdf_goo",
    "diffusion",
]

input_dir_name = "1c-200ps_10iter_output_npt"

traj_range = "-11 -1"

# analysis_origin = Path(
#     "/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-305/nvt-npt-dipole"
# )

ref_path = Path(
    "/home8/yxwu/pGM_water_model/scripts/simulate/old_analysis_local_run_template.sh"
)

if not ref_path.exists():
    raise FileNotFoundError(f"Reference path {ref_path} does not exist")

temperatures = set(range(240, 381, 5))

# temps = [298] + list(temperatures)
temps = list(temperatures)
for temp in temps:
    key_dir = "MD_data-3"
    analysis_origin = Path(
        f"/home8/yxwu/pGM_water_model/{key_dir}/langevin_Berendsen-{temp}/nvt-npt-dipole"
    )

    paths = find_matching_paths(analysis_origin, input_dir_name)

    subpath = extract_subpath_after_key(analysis_origin, key_dir)
    output_path = Path("/home8/yxwu/pGM_water_model/outputs/slurm_script")
    output_folder = output_path / subpath
    output_folder.mkdir(parents=True, exist_ok=True)

    for analysis_type in tqdm(analysis_types):
        for path in paths:
            wandb_folder_path = path.parents[1]
            wandb_name = wandb_folder_path.name
            output_file_path = output_folder / f"{analysis_type}_{wandb_name}.sh"

            replacements = [
                ("{FOLDER_PATH}", str(wandb_folder_path)),
                ("{ANALYSIS_TYPE}", analysis_type),
                ("{TRAJ_RANGE}", traj_range),
                ("{WANDB_NAME}", wandb_name),
            ]

            yxwu_lib.replace_contents_and_save_new(
                str(ref_path), str(output_file_path), replacements
            )

            os.system(f"sh {output_file_path}")
