from dataclasses import dataclass
from smac import Scenario
from smac import HyperparameterOptimizationFacade as HPOFacade


@dataclass
class OptimizerParameters:
    optimizer: str = "HyperbandFacade"

    n_configs_per_hyperparamter: int = 100

    max_tries_for_challenging: int = 512

    # wait for submitted trials
    wait_at_submitted_num: int = -1


def get_setting_for_hyperbandfacade(problem: Scenario, arg: OptimizerParameters):
    init_design = HPOFacade.get_initial_design(
        problem, n_configs_per_hyperparamter=arg.n_configs_per_hyperparamter
    )
    intensifier = HPOFacade.get_intensifier(problem)
    intensifier._retries = arg.max_tries_for_challenging
    return init_design, intensifier
