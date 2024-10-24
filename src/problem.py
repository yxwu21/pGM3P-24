from pathlib import Path
from ConfigSpace import ConfigurationSpace
from smac import Scenario
from dataclasses import dataclass


@dataclass
class OptimizationProblemParameters:
    name: str = "water_model"

    output_directory: str = "outputs"

    deterministic: bool = True

    termination_cost_threshold: float = float("inf")

    crash_cost: float = 100.0

    walltime_limit: float = float("inf")

    # trial time limit in hours
    trial_walltime_limit: float = 12

    # max number of trials
    n_trials: int = 1000

    # min trial time limit in hours
    min_budget: float = 8

    # max trial time limit in hours
    max_budget: float = 12

    # random seed
    seed: int = 2023

    # number of parallel workers
    num_workers: int = 1


class OptimizationProblem:
    def __init__(self, params: OptimizationProblemParameters) -> None:
        self.params = params

    def define(self, configspace: ConfigurationSpace) -> Scenario:
        output_dir = Path(self.params.output_directory)
        output_dir.mkdir(exist_ok=True)

        prob = Scenario(
            configspace,
            name=self.params.name,
            output_directory=output_dir,
            deterministic=self.params.deterministic,
            crash_cost=self.params.crash_cost,
            termination_cost_threshold=self.params.termination_cost_threshold,
            walltime_limit=self.params.walltime_limit,
            # trial_memory_limit=self.params.trial_walltime_limit,
            n_trials=self.params.n_trials,
            min_budget=self.params.min_budget,
            max_budget=self.params.max_budget,
            seed=self.params.seed,
            n_workers=self.params.num_workers,
        )
        return prob
