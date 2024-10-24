import os
import shutil

from dataclasses import dataclass
from typing import Tuple
from ConfigSpace import Configuration, ConfigurationSpace, Float
from .md import MDModelParameters, MDConfigurations, MDRunner
from .eval import eval_md_results
from .simulation import SimulationManager


@dataclass
class WarterModelParameters:
    vdw_a_range: Tuple[float, float] = (6e5, 8e5)
    vdw_a_default: float = 6.24e5

    vdw_b_range: Tuple[float, float] = (6e2, 8e2)
    vdw_b_default: float = 6.72e2

    scale_q_range: Tuple[float, float] = (0.5, 1.5)
    scale_q_default: float = 1.18
    scale_p_range: Tuple[float, float] = (0.5, 1.5)
    scale_p_default: float = 0.90
    scale_pol_range: Tuple[float, float] = (0.5, 1.5)
    scale_pol_default: float = 0.97
    scale_rad_range: Tuple[float, float] = (0.5, 1.5)
    scale_rad_default: float = 0.82

    rdf_score_weight: float = 1.0
    density_score_weight: float = 1.0
    dipole_score_weight: float = 1.0


class WaterModel:
    def __init__(
        self,
        params: WarterModelParameters,
        md_config: MDConfigurations,
    ) -> None:
        self.params = params
        self.md_config = md_config

    def configspace(
        self,
    ) -> ConfigurationSpace:
        cs = ConfigurationSpace(seed=0)
        vdw_a = Float(
            "vdw_a", self.params.vdw_a_range, default=self.params.vdw_a_default
        )
        vdw_b = Float(
            "vdw_b", self.params.vdw_b_range, default=self.params.vdw_b_default
        )
        scale_q = Float(
            "scale_q", self.params.scale_q_range, default=self.params.scale_q_default
        )
        scale_p = Float(
            "scale_p", self.params.scale_p_range, default=self.params.scale_p_default
        )
        scale_pol = Float(
            "scale_pol",
            self.params.scale_pol_range,
            default=self.params.scale_pol_default,
        )
        scale_rad = Float(
            "scale_rad",
            self.params.scale_rad_range,
            default=self.params.scale_rad_default,
        )

        cs.add_hyperparameters([vdw_a, vdw_b, scale_q, scale_p, scale_pol, scale_rad])
        return cs

    def delete_evaluate_cache(self, cache_path: str):
        if os.path.exists(cache_path):
            shutil.rmtree(cache_path)

    def evaluate(self, config: Configuration, seed: int = 0) -> float:
        md_model_params = MDModelParameters(
            vdw_a=config["vdw_a"],
            vdw_b=config["vdw_b"],
            scale_q=config["scale_q"],
            scale_p=config["scale_p"],
            scale_pol=config["scale_pol"],
            scale_rad=config["scale_rad"],
        )

        # create md runner
        md_runner = MDRunner(self.md_config)
        sim_manager = SimulationManager()

        # run md
        md_results, running_folder_path = md_runner.run(
            md_model_params, sim_manager, return_folder_path=True
        )

        # eval md results
        if md_results.success:
            try:
                metrics, data = eval_md_results(md_results)

                # return metrics
                eval_metrics = (
                    metrics["rdf_score"] * self.params.rdf_score_weight
                    + metrics["density_score"] * self.params.density_score_weight
                    + metrics["dipole_score"] * self.params.dipole_score_weight
                )
                metrics["weighted_total_score"] = eval_metrics
                metrics["wandb"] = data
            except:
                raise ValueError("MD eval failed")
            finally:
                # remove evaluation cache
                self.delete_evaluate_cache(running_folder_path)
        else:
            # remove evaluation cache
            self.delete_evaluate_cache(running_folder_path)
            raise ValueError("MD run failed")
        return eval_metrics, metrics
