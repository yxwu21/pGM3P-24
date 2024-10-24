import os

from pathlib import Path
from dataclasses import dataclass, field
from src.md import (
    MDConfigurations,
    MDRunner,
    MDResluts,
    MDModelParameters,
    SimulationManager,
)
from src.barostat_simulation import BarostatSimulationManager
from src.md_pipeline import PipelineMDConfigurations, RDFExtAnalysisMDRunner
from src.temp_simulation import TempSimulationManager


@dataclass
class TempMDConfigurations(PipelineMDConfigurations):
    temp0: int = 298


class SimulationStep:
    def __init__(self, **kwargs):
        self.params = kwargs

    def step_run(self, sim_manager, md_path, comb_filename, previous_output_dir=None):
        self.params.setdefault("restart", False)
        self.params.setdefault("timeout_limit", None)
        if previous_output_dir and self.params.get("restart"):
            self.params["input_dir_name"] = previous_output_dir
        return sim_manager.run_or_restart_simulation(
            md_path, comb_filename, **self.params
        )


class TempRDFExtAnalysisMDRunner(RDFExtAnalysisMDRunner):
    def __init__(self, md_config: TempMDConfigurations) -> None:
        super().__init__(md_config)
        self.pipeline_builder = TempPipelineBuilder(md_config)
        self.pipeline = self.pipeline_builder.pipeline_config()

    def run(
        self,
        model_params: MDModelParameters,
        sim_manager: TempSimulationManager,
        return_folder_path: bool = False,
    ) -> MDResluts:
        """Execute a MD run.

        :param model_params: parameters for MD model definition
        :return: a MDResluts contains all states, e.g. success: True if MD run is
        successful, False otherwise
        """
        current_path = Path.cwd()
        last_output_dir = None  # Initialize last_output_dir before the loop
        try:
            # model setup
            md_path, case_folder_name = self.model_setup(model_params)

            # save current model params
            self.save_current_model_params(
                self.output_path / case_folder_name, model_params
            )

            # Run each step in the pipeline
            for step in self.pipeline:
                last_output_dir = step.step_run(
                    sim_manager, md_path, self.comb_filename, last_output_dir
                )

            # rdf analysis
            rdf_folder = self.rdf_analysis(md_path, "rdf_goo", case_folder_name)
            self.rdf_analysis(md_path, "rdf_goh", case_folder_name)
            self.rdf_analysis(md_path, "rdf_ghh", case_folder_name)
        except:
            md_results = MDResluts(
                success=False,
                output_dir=None,
                dipole_folder=None,
                rdf_analysis_folder=None,
                experimental_rdf_analysis_folder=None,
            )
            os.chdir(current_path)
            return md_results

        md_results = MDResluts(
            success=True,
            output_dir=md_path,
            dipole_folder=md_path / last_output_dir,
            rdf_analysis_folder=Path(rdf_folder),
            experimental_rdf_analysis_folder=Path(
                self.config.experimental_rdf_analysis_folder
            ),
        )
        os.chdir(current_path)

        if return_folder_path:
            return md_results, self.output_path / case_folder_name
        else:
            return md_results


@dataclass
class TempPipelineBuilder:
    md_config: TempMDConfigurations

    def pipeline_config(self):
        thermo_baro = f"{self.md_config.thermostat}_{self.md_config.barostat}"
        condition = [thermo_baro, self.md_config.pipeline_name, self.md_config.exe]

        pipeline = None
        if condition == ["langevin_Berendsen", "nvt-npt-dipole", "pmemd"]:
            pipeline = [
                SimulationStep(
                    ensemble_type="nvt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[0],
                    id="1a",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[0],
                    skinnb=1,
                    barostat=1,
                    dipole_print=False,
                    iterate=False,
                    total_iterations=1,
                    last_iterate=False,
                    last_ensemble=None,
                    last_executable=None,
                    last_iter_num=None,
                    restart=False,
                    temp0=self.md_config.temp0,
                ),  # nvt run for 10ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[1],
                    id="1b",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[1],
                    skinnb=1,
                    barostat=1,
                    dipole_print=False,
                    iterate=True,
                    total_iterations=10,
                    last_iterate=False,
                    last_ensemble="nvt",
                    last_executable="pmemd-pgm",
                    last_iter_num=None,
                    restart=True,
                    temp0=self.md_config.temp0,
                ),  # npt run for 100ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[1],
                    id="1c",
                    ntt=1,
                    gamma_ln=0,
                    skinnb=1,
                    barostat=1,
                    dipole_print=False,
                    iterate=True,
                    total_iterations=10,
                    last_iterate=True,
                    last_ensemble="npt",
                    last_executable="pmemd-pgm",
                    last_iter_num=10,
                    restart=True,
                    temp0=self.md_config.temp0,
                ),  # npt run for 100ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="sander",
                    md_length=self.md_config.md_lengths[2],
                    id="1c",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[2],
                    skinnb=1,
                    barostat=1,
                    dipole_print=True,
                    iterate=False,
                    total_iterations=1,
                    last_iterate=True,
                    last_ensemble="npt",
                    last_executable="pmemd-pgm",
                    last_iter_num=10,
                    restart=True,
                    temp0=self.md_config.temp0,
                ),  # dipole_print run for 1step
                SimulationStep(
                    ensemble_type="npt",
                    executable="sander",
                    md_length=self.md_config.md_lengths[2],
                    id="1d",
                    ntt=1,
                    gamma_ln=0,
                    skinnb=1,
                    barostat=1,
                    dipole_print=True,
                    iterate=False,
                    total_iterations=1,
                    last_iterate=True,
                    last_ensemble="npt",
                    last_executable="pmemd-pgm",
                    last_iter_num=10,
                    restart=True,
                    temp0=self.md_config.temp0,
                ),  # dipole_print run for 1step
            ]
        elif condition == ["langevin_constTemp", "nvt-npt-npt-dipole"]:
            pipeline = [
                SimulationStep(
                    ensemble_type="nvt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[0],
                    id="1a",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[0],
                    skinnb=1,
                    barostat=1,
                    dipole_print=False,
                    iterate=False,
                    total_iterations=1,
                    last_iterate=False,
                    last_ensemble=None,
                    last_executable=None,
                    last_iter_num=None,
                    restart=False,
                ),  # nvt run for 10ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[1],
                    id="1b",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[1],
                    skinnb=1,
                    barostat=1,
                    dipole_print=False,
                    iterate=True,
                    total_iterations=10,
                    last_iterate=False,
                    last_ensemble="nvt",
                    last_executable="pmemd-pgm",
                    last_iter_num=None,
                    restart=True,
                ),  # npt run for 100ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[2],
                    id="1c",
                    ntt=1,
                    gamma_ln=self.md_config.gamma_lns[2],
                    skinnb=1,
                    barostat=1,
                    dipole_print=False,
                    iterate=True,
                    total_iterations=10,
                    last_iterate=True,
                    last_ensemble="npt",
                    last_executable="pmemd-pgm",
                    last_iter_num=None,
                    restart=True,
                ),  # npt run for 100ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="sander",
                    md_length=self.md_config.md_lengths[3],
                    id="1d",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[3],
                    skinnb=1,
                    barostat=1,
                    dipole_print=True,
                    iterate=False,
                    total_iterations=1,
                    last_iterate=True,
                    last_ensemble="npt",
                    last_executable="pmemd-pgm",
                    last_iter_num=10,
                    restart=True,
                ),  # dipole_print run for 1step
            ]
        elif condition == ["langevin_MC", "nvt-npt-dipole", "pmemd"]:
            pipeline = [
                SimulationStep(
                    ensemble_type="nvt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[0],
                    id="1a",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[0],
                    skinnb=1,
                    barostat=2,
                    dipole_print=False,
                    iterate=False,
                    total_iterations=1,
                    last_iterate=False,
                    last_ensemble=None,
                    last_executable=None,
                    last_iter_num=None,
                    restart=False,
                    temp0=self.md_config.temp0,
                ),  # nvt run for 10ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[1],
                    id="1b",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[1],
                    skinnb=1,
                    barostat=2,
                    dipole_print=False,
                    iterate=True,
                    total_iterations=10,
                    last_iterate=False,
                    last_ensemble="nvt",
                    last_executable="pmemd-pgm",
                    last_iter_num=None,
                    restart=True,
                    temp0=self.md_config.temp0,
                ),  # npt run for 100ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="pmemd-pgm",
                    md_length=self.md_config.md_lengths[1],
                    id="1c",
                    ntt=1,
                    gamma_ln=0,
                    skinnb=1,
                    barostat=2,
                    dipole_print=False,
                    iterate=True,
                    total_iterations=10,
                    last_iterate=True,
                    last_ensemble="npt",
                    last_executable="pmemd-pgm",
                    last_iter_num=10,
                    restart=True,
                    temp0=self.md_config.temp0,
                ),  # npt run for 100ps
                SimulationStep(
                    ensemble_type="npt",
                    executable="sander",
                    md_length=self.md_config.md_lengths[2],
                    id="1c",
                    ntt=3,
                    gamma_ln=self.md_config.gamma_lns[2],
                    skinnb=1,
                    barostat=2,
                    dipole_print=True,
                    iterate=False,
                    total_iterations=1,
                    last_iterate=True,
                    last_ensemble="npt",
                    last_executable="pmemd-pgm",
                    last_iter_num=10,
                    restart=True,
                    temp0=self.md_config.temp0,
                ),  # dipole_print run for 1step
                SimulationStep(
                    ensemble_type="npt",
                    executable="sander",
                    md_length=self.md_config.md_lengths[2],
                    id="1d",
                    ntt=1,
                    gamma_ln=0,
                    skinnb=1,
                    barostat=2,
                    dipole_print=True,
                    iterate=False,
                    total_iterations=1,
                    last_iterate=True,
                    last_ensemble="npt",
                    last_executable="pmemd-pgm",
                    last_iter_num=10,
                    restart=True,
                    temp0=self.md_config.temp0,
                ),  # dipole_print run for 1step
            ]

        return pipeline