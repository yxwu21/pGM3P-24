import tyro

from mdtool.slurm import SlurmArgs, slurm_launcher
from dataclasses import dataclass
from src.md import (
    MDRunner,
    MDModelParameters,
    MDConfigurations,
)
from src.simulation import SimulationManager
from src.get_wandb_data import process_all_csv_files
from src.md_pipeline import PipelineMDConfigurations, RDFExtAnalysisMDRunner
from src.barostat_simulation import BarostatSimulationManager
from src.temp_simulation import TempSimulationManager
from src.temp_md_pipeline import TempMDConfigurations, TempRDFExtAnalysisMDRunner


@dataclass
class Args:
    md_config: TempMDConfigurations

    slurm: SlurmArgs

    md_param_key: str

    md_param_step: int

    directory_path: str


def md_simulation(args: Args, wandb_data: dict):
    # md_config = MDConfigurations(
    #     case=args.case, output_dir=args.output_dir, md_lengths=args.md_lengths
    # )
    md_config = args.md_config
    # runner = MDRunner(md_config)
    # runner = RDFExtAnalysisMDRunner(md_config)
    runner = TempRDFExtAnalysisMDRunner(md_config)

    step_wise_data = wandb_data.get(args.md_param_key, {})
    model_params = MDModelParameters(
        vdw_a=step_wise_data.get("config.vdw_a").get(args.md_param_step),
        vdw_b=step_wise_data.get("config.vdw_b").get(args.md_param_step),
        scale_q=step_wise_data.get("config.scale_q").get(args.md_param_step),
        scale_p=step_wise_data.get("config.scale_p").get(args.md_param_step),
        scale_pol=step_wise_data.get("config.scale_pol").get(args.md_param_step),
        scale_rad=step_wise_data.get("config.scale_rad").get(args.md_param_step),
    )
    # sim_manager = SimulationManager()
    # sim_manager = BarostatSimulationManager()
    sim_manager = TempSimulationManager()
    md_result = runner.run(model_params=model_params, sim_manager=sim_manager)
    return md_result


def main(args: Args):
    wandb_data = process_all_csv_files(args.directory_path)
    md_result = md_simulation(args=args, wandb_data=wandb_data)
    return md_result


@slurm_launcher(Args)
def slurm_main(args: Args):
    wandb_data = process_all_csv_files(args.directory_path)
    md_result = md_simulation(args=args, wandb_data=wandb_data)
    return md_result


if __name__ == "__main__":
    # args = tyro.cli(Args)
    # main()

    slurm_main()
