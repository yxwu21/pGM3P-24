import os

from pathlib import Path

temp0 = 298

proj_path = Path("/home/yxwu/pGM_water_model")
main_path = proj_path / "datasets" / "benchmark"
md_data_path = proj_path / "MD_data"
ref_path = proj_path / "scripts" / "md" / "benchmark"

wat_model = "tip3p"
wat_num = "512"

print(wat_model.upper())
file_name = f"{wat_model}_{wat_num}wat"


case_main_path = md_data_path / wat_model.upper() / wat_num

out_path = md_data_path / "benchmark" / f"{file_name}-{temp0}"
prep_path = out_path / "prep_files"
prep_path.mkdir(exist_ok=True, parents=True)
md_path = out_path / "MD"
md_path.mkdir(exist_ok=True, parents=True)

os.system(f"cp {case_main_path}/{file_name}.prmtop {prep_path}")
os.system(f"cp {case_main_path}/{file_name}.inpcrd {prep_path}")
