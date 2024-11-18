import os

from pathlib import Path
from src.param_dict import step_dict
from src.pipeline_settings_dict import (
    pipeline_dict,
    THERMOSTAT,
    BAROSTAT,
    PIPELINE_NAME,
)
from mdtool import yxwu_lib

# THERMOSTAT = "langevin"
# BAROSTAT = "Berendsen"
# PIPELINE_NAME = "nvt-npt-dipole"
# MD_LENGTHS = "10step 10step 1step"
# GAMMA_LNS = "100 100 100"

thermo_baro = f"{THERMOSTAT}_{BAROSTAT}"
MD_LENGTHS = pipeline_dict[thermo_baro][PIPELINE_NAME]["MD_LENGTHS"]
GAMMA_LNS = pipeline_dict[thermo_baro][PIPELINE_NAME]["GAMMA_LNS"]

# temperatures1 = set(range(240, 381, 20))
# temperatures2 = set(range(240, 381, 10))
temperatures3 = set(range(240, 381, 5))

# non_overlapping_temperatures = temperatures3 - temperatures2 - temperatures1

# non_overlapping_temperatures = list(non_overlapping_temperatures)
# print(non_overlapping_temperatures)

WAT_MODEL = "TIP3P"
WAT_NUM = "512"

for temp in temperatures1:

    TEMP0 = temp
    print(temp, ":")

    ref_path = Path(
        "/home8/yxwu/pGM_water_model/scripts/simulate/slurm_run_template_bench.sh"
    )

    output_path = Path("/home8/yxwu/pGM_water_model/outputs/slurm_script")
    md_output_folder = Path(
        f"/home8/yxwu/pGM_water_model/MD_data/{thermo_baro}-{str(TEMP0)}/{PIPELINE_NAME}"
    )
    md_output_folder.mkdir(parents=True, exist_ok=True)

    for key in step_dict:
        for step in step_dict[key]:

            output_folder = (
                output_path
                / "Benchmark"
                / f"wat-{WAT_NUM}"
                / f"{thermo_baro}-{str(TEMP0)}"
                / PIPELINE_NAME
            )
            output_folder.mkdir(parents=True, exist_ok=True)
            output_file_path = output_folder / f"{WAT_MODEL}.sh"

            replacements = [
                ("{THERMOSTAT}", THERMOSTAT),
                ("{BAROSTAT}", BAROSTAT),
                ("{PIPELINE_NAME}", PIPELINE_NAME),
                ("{MD_LENGTHS}", MD_LENGTHS),
                ("{GAMMA_LNS}", GAMMA_LNS),
                ("{TEMP0}", str(TEMP0)),
                ("{WAT_MODEL}", WAT_MODEL),
                ("{WAT_NUM}", WAT_NUM),
            ]

            yxwu_lib.replace_contents_and_save_new(
                str(ref_path), str(output_file_path), replacements
            )

            os.system(f"sh {output_file_path}")
            break
        break
    break
    print("-" * 180)
