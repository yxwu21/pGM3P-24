from pathlib import Path
from dataclasses import dataclass
from mdtool.slurm import SlurmArgs, slurm_launcher
from src.temp_simulation import TempSimulationManager


@dataclass
class Args:
    slurm: SlurmArgs

    md_path: str = "/home/yxwu/pGM_water_model/MD_data/langevin/vital-cosmos-120_144/MD"
    case: str = "case_1"
    file_name: str = "pgm_512wat"
    input_dir_name: str = "1c-200ps_10iter_output_npt"
    md_length: str = "200ps"
    executable: str = "sander"
    temp0: int = 298


def increment_id(input_id: str) -> str:
    # Splitting the ID into parts: numeric, letter, and the rest
    parts = input_id.split("-")
    numeric, letter = parts[0][:-1], parts[0][-1]

    # Increment the letter
    if letter != "z":
        # Increment the letter if it's not 'z'
        next_letter = chr(ord(letter) + 1)
    else:
        # If the letter is 'z', roll over to 'a' and increment the number
        next_letter = "a"
        numeric = str(int(numeric) + 1)

    # Reconstruct the ID
    next_id = f"{numeric}{next_letter}"
    return next_id


@slurm_launcher(Args)
def main(args: Args):
    simulation = TempSimulationManager()

    simulation.run_or_restart_simulation(
        Path(args.md_path),
        f"{args.case}_{args.file_name}",
        ensemble_type="npt",
        executable="sander",
        md_length=args.md_length,
        id=increment_id(args.input_dir_name),
        ntt=3,
        gamma_ln=2,
        skinnb=1,
        barostat=1,
        dipole_print=True,
        iterate=True,
        total_iterations=10,
        last_iterate=True,
        last_ensemble="npt",
        last_executable="pmemd-pgm",
        last_iter_num=10,
        restart=True,
        input_dir_name=args.input_dir_name,
        temp0=args.temp0,
    )


if __name__ == "__main__":
    # args = tyro.cli(Args)
    main()
