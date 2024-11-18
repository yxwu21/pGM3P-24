import os
import glob
import numpy as np
from pathlib import Path


class DielectricExtractor:

    def __init__(self, proj_path, case, wat_num):
        self.proj_path = Path(proj_path)
        self.case = case
        self.wat_num = wat_num
        self.temp_dielectric = {}

    def extract_last_number(self, file_path):
        try:
            with open(file_path, "r") as file:
                lines = file.readlines()
                if lines:
                    last_line = lines[-1].strip()
                    last_number_str = last_line.split()[-1]
                    try:
                        return float(last_number_str)
                    except ValueError:
                        print(
                            f"The value '{last_number_str}' is not a valid float in file: {file_path}"
                        )
                        return None
                else:
                    print(f"File {file_path} is empty.")
                    return None
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
            return None

    def get_file_paths(self, folder):
        temp_folders = glob.glob(f"{folder}/{self.case}/{self.wat_num}/*")
        return temp_folders

    def process_temp_folder(self, temp_folder):
        temp_folder_path = Path(temp_folder)
        temp = temp_folder_path.name.split("-")[-1]
        prod_folder = temp_folder_path / "MD" / "3_Prod"
        file_paths = glob.glob(f"{prod_folder}/fort.100_*")

        dielectric_values = [
            self.extract_last_number(file_path)
            for file_path in file_paths
            if self.extract_last_number(file_path) is not None
        ]

        return temp, dielectric_values

    def process_folder(self, folder):
        folder_name = Path(folder).name
        self.initialize_nested_dict(folder_name)

        temp_folders = self.get_file_paths(folder)
        for temp_folder in temp_folders:
            temp, dielectric_values = self.process_temp_folder(temp_folder)
            if dielectric_values:
                avg_temp_dielectric = np.average(dielectric_values)
                self.temp_dielectric[folder_name][self.case][temp] = avg_temp_dielectric
            else:
                print(f"No valid dielectric values found in {temp_folder}.")

    def initialize_nested_dict(self, folder_name):
        if folder_name not in self.temp_dielectric:
            self.temp_dielectric[folder_name] = {}
        if self.case not in self.temp_dielectric[folder_name]:
            self.temp_dielectric[folder_name][self.case] = {}

    def run(self):
        folder_pattern = self.proj_path / "MD_data-pGM-temp*"
        folders = glob.glob(str(folder_pattern))

        for folder in folders:
            self.process_folder(folder)

        return self.temp_dielectric


if __name__ == "__main__":
    proj_path = "/home8/yxwu/pGM_water_model"
    case = "vital-cosmos-120_340"
    wat_num = "512"

    extractor = DielectricExtractor(proj_path, case, wat_num)
    result = extractor.run()

    print(result)
