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

ref_path = Path(
    "/home/yxwu/pGM_water_model/scripts/simulate/slurm_run_from_wandb_template.sh"
)

output_path = Path("/home/yxwu/pGM_water_model/outputs/slurm_script")
md_output_folder = Path(
    f"/home/yxwu/pGM_water_model/MD_data/{thermo_baro}/{PIPELINE_NAME}"
)
md_output_folder.mkdir(parents=True, exist_ok=True)

for key in step_dict:
    for step in step_dict[key]:
        log_step = f"{key}_{step}"

        output_folder = output_path / f"{thermo_baro}" / PIPELINE_NAME
        output_folder.mkdir(parents=True, exist_ok=True)
        output_file_path = output_folder / f"{log_step}.sh"

        replacements = [
            ("{WANDB_LOG_NAME}", str(key)),
            ("{LOG_STEP}", str(step)),
            ("{THERMOSTAT}", THERMOSTAT),
            ("{BAROSTAT}", BAROSTAT),
            ("{PIPELINE_NAME}", PIPELINE_NAME),
            ("{MD_LENGTHS}", MD_LENGTHS),
            ("{GAMMA_LNS}", GAMMA_LNS),
        ]

        yxwu_lib.replace_contents_and_save_new(
            str(ref_path), str(output_file_path), replacements
        )

        # os.system(f"sh {output_file_path}")
    #     break
    # break
