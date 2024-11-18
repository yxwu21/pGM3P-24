import tyro

from dataclasses import dataclass, field, asdict
from smac import HyperparameterOptimizationFacade as HPOFacade
from src.md import MDConfigurations
from src.model import WarterModelParameters, WaterModel
from src.problem import OptimizationProblemParameters, OptimizationProblem
from src.optim import OptimizerParameters, get_setting_for_hyperbandfacade
from src.cluster import ClusterConfigurations, Cluster
from src.wandb import WandBConfig, WandBCallback
from src.visualization import visualize_rdf
from src.callback import WaitRunningTrialsAfterSubmmision


@dataclass
class Experiments:
    # use cluster
    use_cluster: bool = False

    # water model settings
    model: WarterModelParameters = field(default_factory=WarterModelParameters)

    # optimization problem settings
    problem: OptimizationProblemParameters = field(
        default_factory=OptimizationProblemParameters
    )

    # optimizer settings
    optimizer: OptimizerParameters = field(default_factory=OptimizerParameters)

    # configurations of cluster system
    cluster: ClusterConfigurations = field(default_factory=ClusterConfigurations)

    # configurations of md run
    md: MDConfigurations = field(default_factory=MDConfigurations)

    # logger
    wandb: WandBConfig = field(default_factory=WandBConfig)


if __name__ == "__main__":
    arg = tyro.cli(Experiments)

    # define water model and problem
    water_model = WaterModel(arg.model, arg.md)
    problem = OptimizationProblem(arg.problem).define(water_model.configspace())

    # assign task client
    client = None
    if arg.use_cluster:
        client = Cluster(arg.cluster).get_client()

    # Now we use SMAC to find the best hyperparameters
    wandb_callback = WandBCallback(
        arg.wandb,
        config_to_save=asdict(arg),
        wandb_function=visualize_rdf,
        wandb_function_step=1,
    )
    wait_callback = WaitRunningTrialsAfterSubmmision(
        arg.optimizer.wait_at_submitted_num
    )
    init_design, intensifier = get_setting_for_hyperbandfacade(problem, arg.optimizer)
    optimizer = HPOFacade(
        problem,
        water_model.evaluate,  # We pass the target function here
        initial_design=init_design,
        intensifier=intensifier,
        overwrite=True,  # Overrides any previous results that are found that are inconsistent with the meta-data
        dask_client=client,
        callbacks=[wandb_callback, wait_callback],
        # logging_level=1,
    )

    # Start optimization
    incumbent = optimizer.optimize()
