import os
import glob

from pathlib import Path
from tqdm import tqdm
from mdtool import yxwu_lib


analysis_types = [
    "rdf_ghh",
    "rdf_goh",
    "rdf_goo",
    "diffusion",
]

slurm_run = True

# input_folder = "MD_data-pGM-temp-3"
# input_folder = "MD_data_Benchamrk-pgm_mpi-2"
input_folders = [
    "MD_data_Benchamrk-pgm_mpi-3",
    "MD_data_Benchamrk-long-pgm_mpi",
    "MD_data_Benchamrk-long-pgm_mpi-1",
    "MD_data_Benchamrk-long-pgm_mpi-2",
    "MD_data_Benchamrk-long-pgm_mpi-3",
]

last_n_trajectories = 10

input_dir_name = "3_Prod"


proj_path = Path("/home8/yxwu/pGM_water_model")

if slurm_run:
    mode = "slurm"
else:
    mode = "local"

ref_path = proj_path / "scripts" / "simulate" / f"analysis_{mode}_run_template.sh"

if not ref_path.exists():
    raise FileNotFoundError(f"Reference path {ref_path} does not exist")

for input_folder in input_folders:

    wandb_folder_paths = glob.glob(f"{proj_path}/{input_folder}/*/*/*")

    output_path = proj_path / "outputs" / "slurm_script"
    output_folder = output_path / input_folder
    output_folder.mkdir(parents=True, exist_ok=True)

    for analysis_type in tqdm(analysis_types):
        for wandb_folder_path in wandb_folder_paths:
            wandb_name = Path(wandb_folder_path).name
            output_file_path = output_folder / f"{analysis_type}_{wandb_name}.sh"

            replacements = [
                ("{FOLDER_PATH}", str(wandb_folder_path)),
                ("{ANALYSIS_TYPE}", analysis_type),
                ("{LAST_N_TRAJ}", str(last_n_trajectories)),
                ("{WANDB_NAME}", wandb_name),
            ]

            yxwu_lib.replace_contents_and_save_new(
                str(ref_path), str(output_file_path), replacements
            )

            os.system(f"sh {output_file_path}")
        #     break
        # break
