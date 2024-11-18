import os
import glob
from pathlib import Path


def check_section_in_file(file_path, section_name):
    """Check if a specific section exists in a file."""
    try:
        with open(file_path, "r") as file:
            for line in file:
                if section_name in line:
                    return True
        return False
    except FileNotFoundError:
        print(f"File not found: {file_path}")
        return False
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return False


def process_files(folder_path, section_name="5.  TIMINGS"):
    """Process all .out files in a directory, checking for a specific section."""
    finish_list = []
    not_finish_list = []

    out_file_paths = Path(folder_path).glob("*/MD/3*/*.out")
    for out_file_path in out_file_paths:
        case = out_file_path.parents[2].name
        if check_section_in_file(out_file_path, section_name):
            finish_list.append(case)
        else:
            not_finish_list.append(case)
            print(f"Section not found in: {out_file_path}")

    return finish_list, not_finish_list


if __name__ == "__main__":
    dataset_path = Path("/home8/yxwu/pGM_water_model")
    folder_path = dataset_path / "MD_data-pGM-temp-3" / "vital-cosmos-120_340" / "512"

    finish_list, not_finish_list = process_files(folder_path)

    print("Cases with finished section:")
    print(finish_list)
    print(f"Number of finished cases found: {len(finish_list)}")

    print("Cases without finished section:")
    print(not_finish_list)
    print(f"Number of not finished cases found: {len(not_finish_list)}")
