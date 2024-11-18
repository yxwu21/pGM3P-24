import os
import glob

from pathlib import Path
from src.param_dict import step_dict

# # folders = os.listdir("/tmp/yongxian_pgm_wat_cadidates")
# # print(folders)

# # for folder in folders:
# #     os.system(
# #         f"cp /home/yxwu/pGM_water_model/MD_data/langevin_Berendsen-pmemd_backup/nvt-npt-dipole/{folder}/model_construct/modi_qin.qout /tmp/yongxian_pgm_wat_cadidates/{folder}"
# #     )


# # for key in step_dict:
# #     for step in step_dict[key]:
# #         log_step = f"{key}_{step}"
# #         full_path = Path(
# #             "/home/yxwu/pGM_water_model/MD_data/langevin_Berendsen_full-298/nvt-npt-dipole"
# #         )
# #         file_path = full_path / log_step
# #         folder_path = Path(
# #             "/home/yxwu/pGM_water_model/MD_data/langevin_Berendsen-298/nvt-npt-dipole"
# #         )
# #         folder_path.mkdir(exist_ok=True)
# #         os.system(f"cp -r {file_path} {folder_path}")

# temps = [
#     # "280",
#     # "298",
#     # "300",
# ]
# cases = [
#     # "vital-cosmos-120_320",
#     # "vital-cosmos-120_338",
#     # "vital-cosmos-120_387",
#     # "vital-cosmos-120_360",
#     # "vital-cosmos-120_383",
#     # "vital-cosmos-120_343",
#     # "vital-cosmos-120_521",
#     # "vital-cosmos-120_578",
#     # "vital-cosmos-120_659",
#     # "vital-cosmos-120_687",
#     # "vital-cosmos-120_717",
#     # "vital-cosmos-120_730",
#     # "vital-cosmos-120_742",
#     # "vital-cosmos-120_758",
#     # "vague-sun-127_115",
#     # "atomic-puddle-128_177",
#     # "atomic-puddle-128_471",
#     # "atomic-puddle-128_226",
#     # "vital-cosmos-120_1455",
#     # "atomic-puddle-128_755",
#     # "desert-star-137_126",
#     # "atomic-puddle-128_1051",
#     # "eternal-terrain-135_109",
#     "vital-cosmos-120_325",
#     "vital-cosmos-120_340",
#     "vital-cosmos-120_386",
#     "vital-cosmos-120_611",
#     "vital-cosmos-120_707",
#     "vital-cosmos-120_737",
#     "vital-cosmos-120_993",
#     "vague-sun-127_42",
#     "atomic-puddle-128_51",
#     "atomic-puddle-128_256",
#     "atomic-puddle-128_389",
#     "atomic-puddle-128_672",
#     "atomic-puddle-128_752",
#     "atomic-puddle-128_977",
#     "atomic-puddle-128_1182",
#     "eager-meadow-129_374",
#     "eager-meadow-129_485",
#     "desert-star-137_111",
# ]
# for temp in temps:
#     md_path = Path(f"MD_data/langevin_Berendsen-{temp}")
#     print(md_path.parts[-1])
#     target_path = Path("temp_analysis_2")
#     for case in cases:
#         target_case_path = target_path / md_path.parts[-1] / case / "MD" / "prep_files"
#         target_case_path.mkdir(exist_ok=True, parents=True)

#         os.system(
#             f"cp {md_path}/nvt-npt-dipole/{case}/MD/prep_files/case_1_pgm_512wat.prmtop {target_case_path}/case_1_pgm_512wat.prmtop"
#         )

# src_bench_path = Path("MD_data_Benchmark-MC-Langevin")
# bench_path = Path("MD_data_Benchmark_reorg-MC-Langevin")
# bench_path.mkdir(exist_ok=True, parents=True)

# wat_models = [
#     # "TIP3P",
#     # "TIP4P",
#     # "TIP5P",
#     "OPC",
#     "OPC3",
#     "SPCE",
#     "TIP4PEW",
# ]
# wat_num = "512"
# temps = ["298"] + [str(t) for t in range(240, 381, 5)]
# for temp in temps:
#     target_temp_path = bench_path / f"Benchmark-{temp}" / wat_num
#     target_temp_path.mkdir(parents=True, exist_ok=True)
#     for wat_model in wat_models:
#         src_wat_model_path = (
#             src_bench_path / wat_model / wat_num / f"{wat_model}_{wat_num}wat-{temp}"
#         )
#         target_wat_model_path = target_temp_path / wat_model
#         os.system(f"cp -r {src_wat_model_path} {target_wat_model_path}")


# pattern = "/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-250/nvt-npt-dipole/vital-cosmos-120_*/MD/prep_files"
# cases = ["vague-sun-127_115", "atomic-puddle-128_1182", "atomic-puddle-128_672"]
# cases = ["vital-cosmos-120_320", "vital-cosmos-120_340", "vital-cosmos-120_578"]
cases = ["vital-cosmos-120_340"]
# for case in cases:
#     pattern = f"/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-298/nvt-npt-dipole/{case}"
#     target_path = Path("/home8/yxwu/pGM_water_model/MD_data-0")
#     paths = glob.glob(pattern)
#     for path in paths:
#         temp_folder = Path(path).parts[5]
#         case_folder = Path(path).parts[7]
#         target_path_folder = target_path / temp_folder / "nvt-npt-dipole"
#         target_path_folder.mkdir(exist_ok=True, parents=True)
#         os.system(f"cp -r {path} {target_path_folder}")
#         # print(path)


path = "/home8/yxwu/pGM_water_model/MD_data-0"
folders = os.listdir(path)
for folder in folders:
    temp = folder.split("-")[-1]
    print(temp)
    # for case in cases:
    #     pattern = f"/home8/yxwu/pGM_water_model/MD_data-0/langevin_Berendsen-240/nvt-npt-dipole/{case}"
