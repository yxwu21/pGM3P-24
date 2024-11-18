import os

from pathlib import Path
from mdtool import yxwu_lib
from rst_simulation import increment_id


def find_matching_paths(rst_origin, input_dir_name):
    rst_origin_path = Path(rst_origin)
    pattern = f"*/MD/{input_dir_name}"
    matching_paths = list(rst_origin_path.glob(pattern))
    if not matching_paths:
        raise ValueError(f"No matching paths found for the pattern {pattern}")
    return matching_paths


def check_final_performance_info_in_file(filepath):
    target_line = "Final Performance Info:"
    try:
        with open(filepath, "r") as file:
            for line in file:
                # if line.strip() == target_line:
                if target_line in line.strip():
                    return True
        return False
    except FileNotFoundError:
        print(f"File not found: {filepath}")
        return False
    except Exception as e:
        print(f"An error occurred: {e}")
        return False


def find_rst_paths(rst_origin, input_dir_name, pattern):
    rst_origin_path = Path(rst_origin)
    # pattern = f"*/MD/{input_dir_name}/*.10.out"
    matching_paths = list(rst_origin_path.glob(pattern))

    if not matching_paths:
        raise ValueError(f"No matching paths found for the pattern {pattern}")

    rst_paths = []
    not_rst_paths = []

    for matching_path in matching_paths:
        md_path = matching_path.parents[1]

        # Check if there are any directories or files matching the new pattern
        if list(Path(md_path).glob("*" + increment_id(input_dir_name) + "*")):
            not_rst_paths.append(matching_path)
        else:
            if check_final_performance_info_in_file(str(matching_path)):
                rst_paths.append(matching_path.parents[0])
            else:
                not_rst_paths.append(matching_path.parents[0])

    return rst_paths, not_rst_paths


def extract_subpath_after_key(full_path: Path, key_directory: str) -> str:
    path_parts = full_path.parts
    try:
        index = path_parts.index(key_directory)
    except ValueError:
        raise ValueError(f"Key directory '{key_directory}' not found in path")
    relevant_path_parts = path_parts[index + 1 :]
    relevant_path = "/".join(relevant_path_parts)
    return relevant_path


if __name__ == "__main__":
    md_length = "200ps"

    DIPOLE_PRINT_INTERVAL = 10000
    LAST_EXECUTABLE = "sander"
    LAST_ITERATE = "False"

    rst_wandb_list = [
        # "vital-cosmos-120_320",
        "vital-cosmos-120_340",
        # "vital-cosmos-120_578",
    ]

    temperatures = ["298"] + [str(t) for t in range(240, 381, 5)]

    key_dir = "MD_data"
    for temp in temperatures:
        TEMP0 = temp
        print(temp, ":")

        rst_origin = Path(
            f"/home8/yxwu/pGM_water_model/{key_dir}/langevin_Berendsen-{temp}/nvt-npt-dipole"
        )
        # input_dir_name = "1c-200ps_10iter_output_npt"
        # input_dir_name = "dipole_1h-200ps_output_npt"
        # input_dir_name = "1b-200ps_10iter_output_npt"
        input_dir_name = "dipole_1j-200ps_output_npt"

        id = increment_id(input_dir_name)

        ref_path = Path(
            "/home8/yxwu/pGM_water_model/scripts/simulate/rst_slurm_run_from_wandb_template_temp.sh"
        )

        # ref_path = Path(
        #     "/home8/yxwu/pGM_water_model/scripts/simulate/rst_slurm_run_from_wandb_template_temp_dipole.sh"
        # )

        if not ref_path.exists():
            raise FileNotFoundError(f"Reference path {ref_path} does not exist")

        paths = find_matching_paths(rst_origin, input_dir_name)
        pattern = f"*/MD/{input_dir_name}/*.out"
        rst_paths, not_rst_paths = find_rst_paths(rst_origin, input_dir_name, pattern)
        # print(rst_paths)

        subpath = extract_subpath_after_key(rst_origin, key_dir)
        output_path = Path("/home8/yxwu/pGM_water_model/outputs/slurm_script")
        output_folder = output_path / subpath
        output_folder.mkdir(parents=True, exist_ok=True)

        for rst_path in rst_paths:
            wandb_folder_path = rst_path.parents[1]
            wandb_name = wandb_folder_path.name
            if wandb_name in rst_wandb_list:
                output_file_path = output_folder / f"{id}-rst_{wandb_name}.sh"

                replacements = [
                    ("{FOLDER_PATH}", str(wandb_folder_path)),
                    ("{INPUT_DIR_NAME}", input_dir_name),
                    ("{MD_LENGTH}", md_length),
                    ("{WANDB_NAME}", wandb_name),
                    ("{TEMP0}", str(TEMP0)),
                    ("{DIPOLE_PRINT_INTERVAL}", str(DIPOLE_PRINT_INTERVAL)),
                    ("{LAST_EXECUTABLE}", str(LAST_EXECUTABLE)),
                    ("{LAST_ITERATE}", str(LAST_ITERATE)),
                ]

                yxwu_lib.replace_contents_and_save_new(
                    str(ref_path), str(output_file_path), replacements
                )

                # os.system(f"sh {output_file_path}")
                print(output_file_path)

        #     break
        # break
