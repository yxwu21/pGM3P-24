import tyro

from mdtool.slurm import SlurmArgs, slurm_launcher
from dataclasses import dataclass

from src.bench_md_pipeline import BenchMDConfigurations, BenchRDFExtAnalysisMDRunner
from src.bench_simulation import BenchSimulationManager


@dataclass
class Args:
    md_config: BenchMDConfigurations

    slurm: SlurmArgs


def md_simulation(args: Args):

    md_config = args.md_config

    runner = BenchRDFExtAnalysisMDRunner(md_config)

    sim_manager = BenchSimulationManager()
    md_result = runner.run(sim_manager=sim_manager)
    return md_result


def main(args: Args):
    md_result = md_simulation(args=args)
    return md_result


@slurm_launcher(Args)
def slurm_main(args: Args):
    md_result = md_simulation(args=args)
    return md_result


if __name__ == "__main__":
    # args = tyro.cli(Args)
    # main()

    slurm_main()
