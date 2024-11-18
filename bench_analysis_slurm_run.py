import os

from pathlib import Path
from tqdm import tqdm
from mdtool import yxwu_lib
from rst_slurm_run import find_matching_paths, extract_subpath_after_key


analysis_types = [
    "rdf_ghh",
    "rdf_goh",
    "rdf_goo",
    "diffusion",
]

traj_range = "-10 21"

bench_folder = "MD_data_Benchmark_reorg-MC"
input_dir_name = "4_Prod"
wat_num = "512"

ref_path = Path(
    "/home8/yxwu/pGM_water_model/scripts/simulate/analysis_slurm_run_template.sh"
)

if not ref_path.exists():
    raise FileNotFoundError(f"Reference path {ref_path} does not exist")

temperatures = set(range(240, 381, 5))

temps = [298] + list(temperatures)
for temp in tqdm(temps):
    analysis_origin = Path(
        f"/home8/yxwu/pGM_water_model/{bench_folder}/Benchmark-{temp}/{wat_num}"
    )

    paths = find_matching_paths(analysis_origin, input_dir_name)
    print(paths)

    key_dir = bench_folder
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
    #         break
    #     break
    # break
