import os

from collections import defaultdict
from tqdm import tqdm

import sys
sys.path.append('/home/yxwu/local/script/')
import yxwu_lib


def error_check(path, hit_word):
    print("*"*160)
    print("Folders under error checking...")
    cases = []
    folders = os.listdir(path)
    sorted_folders = sorted(folders)
    for sorted_folder in sorted_folders:
        if f'{hit_word}' in sorted_folder:
            print(sorted_folder)
            cases.append(sorted_folder)
    sorted_cases = sorted(cases)
    print("*"*160)
    print("Simulation error output checking...")
    total_num = 0
    error_num = 0
    for sorted_case in sorted_cases:
        total_num += 1
        MD_output_dirs = os.listdir(f"{path}/{sorted_case}/MD")
        for MD_output_dir in MD_output_dirs:
            if '.run.e' in MD_output_dir:
                error_file_path = f"{path}/{sorted_case}/MD/{MD_output_dir}"
                lines = yxwu_lib.read_file(error_file_path)
                error_lines = []
                for line in lines:
                    if line == '':
                        continue
                    else:
                        error_lines.append(line)
                        # print(line)
                if len(error_lines) == 0:
                    print(f"Error message in {sorted_case}: NO ERROR")
                else:
                    error_num += 1
                    print(f"Error message in {sorted_case}:")
                    for i in error_lines:
                        print(i)
        print("-"*160)
    print(f"Total simulations: {total_num}")
    print(f"Error simulations: {error_num}")


if __name__ == '__main__':

    path = "/home/yxwu/water_model"
    # time_lens = [
    #     "10fs",
    #     '100ps',
    #     '1ns'
    # ]
    # simulation_data = defaultdict(list)
    # for time_len in time_lens:
    #     folders = os.listdir(path)
    #     for folder in folders:
    #         if time_len in folder:
    #             simulation_data[time_len].append(folder)

    # # print each time len
    # n = 0
    # for time_len in time_lens:
    #     print("*" * 80)
    #     print(f'{time_len} simulations: total_num = {len(simulation_data[time_len])}')
    #     n += len(simulation_data[time_len])
    #     # print sorted dirs
    #     for case in sorted(simulation_data[time_len]):
    #         print(case)
    # print("-" * 80)
    # print(f"Total simulations: {n}")

    hit_word = 'case'
    error_check(path, hit_word)

    # print('Error checking NVE simulations:')
    # hit_word_1 = 'nve'
    # error_check(path, hit_word_1)

    # print('Error checking NVT simulations:')
    # hit_word_2 = 'nvt'
    # error_check(path, hit_word_2)

    # print('Error checking NPT simulations:')
    # hit_word_3 = 'npt'
    # error_check(path, hit_word_3)

    # print('Error checking sander simulations:')
    # hit_word_4 = 'sander'
    # error_check(path, hit_word_4)

    # print('Error checking sander simulations:')
    # hit_word_5 = 'pmemd'
    # error_check(path, hit_word_5)