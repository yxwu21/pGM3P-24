import os
import json
import datetime

from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict
from src.yxwu_lib import (
    edit_bond_length,
    modify_and_write_qin,
    edit_vdw_parms,
    list_folders,
    get_traj_info,
)
from src.constant_dict import bond_val_dict
from src.simulation import SimulationManager


@dataclass
class MDConfigurations:
    """General configurations for MD run.

    All general configurations needed for MD run should be defined here, such as the
    input data, output folder, and execution program, etc.
    """

    case: str = "case_1"

    file_name: str = "pgm_512wat"

    qout_type: str = "wat-ms"

    gamma_ln: int = 100

    skinnb: int = 1

    ntt: int = 3

    output_dir: str = "outputs/md_data"

    analysis_script: Dict[str, str] = field(
        default_factory=lambda: {
            "rdf_goo": "cpptraj_rdf_goo_MPI.sh",
            "rdf_goh": "cpptraj_rdf_goh_MPI.sh",
            "rdf_ghh": "cpptraj_rdf_ghh_MPI.sh",
        }
    )

    analysis_ref: str = "/home/yxwu/pGM_water_model/scripts/analysis"

    md_ref: str = "scripts/md"

    experimental_rdf_analysis_folder: str = (
        "datasets/pgm_512wat_main/expt_Skinner_2013/analysis/rdf_goo"
    )

    md_lengths: list[str] = field(default_factory=lambda: ["10ps", "100ps", "1step"])

    use_custom_case_folder_name: bool = False

    custom_case_folder_name: str = ""


@dataclass
class MDModelParameters:
    """Parameters for MD model definition.

    You should define all parameters needed for MD model definition here, such as hyper
    parameters, model parameters.
    """

    vdw_a: float
    vdw_b: float
    scale_q: float
    scale_p: float
    scale_pol: float
    scale_rad: float


@dataclass
class MDResluts:
    """Results of a MD run."""

    success: bool
    output_dir: str
    dipole_folder: Path
    rdf_analysis_folder: Path
    experimental_rdf_analysis_folder: Path


class MDRunner:
    def __init__(self, md_config: MDConfigurations) -> None:
        self.config = md_config

        self.case = md_config.case
        self.file_name = md_config.file_name
        self.qout_type = md_config.qout_type

        self.comb_filename = f"{self.case}_{self.file_name}"
        self.input_dir: str = f"datasets/{self.file_name}_main"
        self.ref_folder_path = f"datasets/{self.file_name}_ref"
        self.prmtop_file: str = f"{self.file_name}.prmtop"
        self.rst_file: str = f"{self.file_name}.restrt0"
        self.prmtop_file_path: str = f"{self.input_dir}/{self.prmtop_file}"
        self.rst_file_path: str = f"{self.input_dir}/{self.rst_file}"
        self.bond_val = bond_val_dict[self.case]

        if "wat-ms" in self.qout_type:
            self.pgm_qout = f"{self.input_dir}/wat_maqz_opt_maqz_wat_ms_esp.wat-ms.qout"
        elif "v08.e.test.q.chg" in self.qout_type:
            self.pgm_qout = f"{self.input_dir}/v08.e.test.q.chg"
        else:
            self.pgm_qout = f"{self.input_dir}/wat10x10_m6311_wat_opt_matz_wat_sol_pol.{self.qout_type}.qout"

        # create output folder
        self.output_path = Path(self.config.output_dir)
        self.output_path.mkdir(exist_ok=True)

        self.use_custom_case_folder_name = md_config.use_custom_case_folder_name
        self.custom_case_folder_name = md_config.custom_case_folder_name

    def save_current_model_params(
        self, save_path: str, model_params: MDModelParameters
    ) -> None:
        data = asdict(model_params)
        sigma, epsilon, r_min = self.cal_lennard_jones_info(model_params)
        data["sigma"] = sigma
        data["epsilon"] = epsilon
        data["r_min"] = r_min
        with open(save_path / "model_params.json", "w") as f:
            json.dump(data, f, indent=4)

    def run(
        self,
        model_params: MDModelParameters,
        sim_manager: SimulationManager,
        return_folder_path: bool = False,
    ) -> MDResluts:
        """Execute a MD run.

        :param model_params: parameters for MD model definition
        :return: a MDResluts contains all states, e.g. success: True if MD run is
        successful, False otherwise
        """
        current_path = Path.cwd()
        try:
            # model setup
            md_path, case_folder_name = self.model_setup(model_params)

            # save current model params
            self.save_current_model_params(
                self.output_path / case_folder_name, model_params
            )

            # nvt run for 10ps
            nvt_output_dir_1a = sim_manager.run_or_restart_simulation(
                md_path,
                self.comb_filename,
                ensemble_type="nvt",
                executable="pmemd-pgm",
                md_length=self.config.md_lengths[0],
                id="1a",
                ntt=self.config.ntt,
                gamma_ln=self.config.gamma_ln,
                skinnb=self.config.skinnb,
                dipole_print=False,
                iterate=False,
                total_iterations=1,
                last_iterate=False,
                last_ensemble=None,
                last_executable=None,
                last_iter_num=None,
                restart=False,
                input_dir_name=None,
                timeout_limit=None,
            )

            # npt run for 100ps
            npt_output_dir_1b = sim_manager.run_or_restart_simulation(
                md_path,
                self.comb_filename,
                ensemble_type="npt",
                executable="pmemd-pgm",
                md_length=self.config.md_lengths[1],
                id="1b",
                ntt=self.config.ntt,
                gamma_ln=self.config.gamma_ln,
                skinnb=self.config.skinnb,
                dipole_print=False,
                iterate=True,
                total_iterations=10,
                last_iterate=False,
                last_ensemble="nvt",
                last_executable="pmemd-pgm",
                last_iter_num=None,
                restart=True,
                input_dir_name=nvt_output_dir_1a,
                timeout_limit=None,
            )

            # dipole_print run for 1step
            dipole_output_dir_1c = sim_manager.run_or_restart_simulation(
                md_path,
                self.comb_filename,
                ensemble_type="npt",
                executable="sander",
                md_length=self.config.md_lengths[2],
                id="1c",
                ntt=self.config.ntt,
                gamma_ln=self.config.gamma_ln,
                skinnb=self.config.skinnb,
                dipole_print=True,
                iterate=False,
                total_iterations=1,
                last_iterate=True,
                last_ensemble="npt",
                last_executable="pmemd-pgm",
                last_iter_num=10,
                restart=True,
                input_dir_name=npt_output_dir_1b,
                timeout_limit=None,
            )

            # rdf analysis
            rdf_folder = self.rdf_analysis(md_path, "rdf_goo", case_folder_name)
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
            dipole_folder=md_path / dipole_output_dir_1c,
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

    def cal_lennard_jones_info(self, model_params: MDModelParameters):
        sigma: float = (model_params.vdw_a / model_params.vdw_b) ** (1 / 6)
        epsilon: float = (model_params.vdw_b**2) / (4 * model_params.vdw_a)
        r_min = (2 * model_params.vdw_a / model_params.vdw_b) ** (1 / 6)
        return sigma, epsilon, r_min

    def model_setup(self, model_params: MDModelParameters) -> tuple[Path, str]:
        """Build promtop and rst files for MD run.

        file structure:
        e.g.
            output_dir/case_folder_name/MD/*

        :param model_params: model parameters
        :return: md result folder Path(output_dir/case_folder_name/MD/), case_folder_name
        """
        # Getting the current time
        current_time = datetime.datetime.now()
        formatted_time = current_time.strftime("%Y%m%d_%H%M%S")

        if self.use_custom_case_folder_name:
            case_folder_name = self.custom_case_folder_name
        else:
            # case_folder_name: str = f"a{model_params.vdw_a}-b{model_params.vdw_b}-q{model_params.scale_q}-p{model_params.scale_p}-x{model_params.scale_pol}-r{model_params.scale_rad}-{formatted_time}"
            case_folder_name: str = f"a{model_params.vdw_a:.4f}-b{model_params.vdw_b:.4f}-q{model_params.scale_q:.4f}-p{model_params.scale_p:.4f}-x{model_params.scale_pol:.4f}-r{model_params.scale_rad:.4f}-{formatted_time}"

        # create folders
        print(type(case_folder_name))
        case_folder: Path = Path(self.config.output_dir) / case_folder_name
        case_folder.mkdir(exist_ok=True)
        model_folder: Path = case_folder / "model_construct"
        model_folder.mkdir(exist_ok=True)

        a_factor = [model_params.vdw_a, 0, 0]
        b_factor = [model_params.vdw_b, 0, 0]

        # edit bond length
        ref_prmtop_1: str = self.prmtop_file_path
        tar_prmtop_name_1 = "TIP3Pmod.prmtop"
        tar_prmtop_path_1: Path = model_folder / tar_prmtop_name_1
        edit_log_bond = edit_bond_length(
            str(tar_prmtop_path_1), ref_prmtop_1, self.bond_val
        )
        # print("%FLAG BOND_EQUIL_VALUE  :", edit_log_bond)

        # modify qin
        input_qin_path = self.pgm_qout
        output_qin_path = model_folder / "modi_qin.qout"
        modify_and_write_qin(
            input_qin_path,
            str(output_qin_path),
            model_params.scale_q,
            model_params.scale_p,
            model_params.scale_pol,
            model_params.scale_rad,
        )

        # amend TIP3Pmod.prmtop with charges
        prmtop = tar_prmtop_path_1
        prmtop_amend = "amend.TIP3Pmod.prmtop"
        prmtop_out = model_folder / prmtop_amend
        amend_cmd = f"python src/new_amend_pgm_wat.py -q {output_qin_path} -p {prmtop} -o {prmtop_out}"
        os.system(amend_cmd)

        # edit VDW parameters
        factor_type_a = "a"
        ref_prmtop_2 = prmtop_out
        tar_prmtop_name_2 = "amend_vdw_a.TIP3Pmod.prmtop"
        tar_prmtop_path_2 = model_folder / tar_prmtop_name_2
        edit_log_a = edit_vdw_parms(
            str(tar_prmtop_path_2), ref_prmtop_2, a_factor, factor_type_a
        )
        # print("%FLAG LENNARD_JONES_ACOEF:", edit_log_a)

        factor_type_b = "b"
        ref_prmtop_3 = tar_prmtop_path_2
        tar_prmtop_name_3 = "amend_vdw.TIP3Pmod.prmtop"
        tar_prmtop_path_3 = model_folder / tar_prmtop_name_3
        edit_log_b = edit_vdw_parms(
            str(tar_prmtop_path_3), ref_prmtop_3, b_factor, factor_type_b
        )
        # print("%FLAG LENNARD_JONES_BCOEF:", edit_log_b)
        os.system(f"rm {tar_prmtop_path_2}")

        # prepare files for md
        md_path: Path = case_folder / "MD"
        md_path.mkdir(exist_ok=True)
        prep_files_path = md_path / "prep_files"
        prep_files_path.mkdir(exist_ok=True)

        os.system(
            f"cp {self.rst_file_path} {prep_files_path}/{self.comb_filename}.restrt0"
        )
        os.system(
            f"cp {tar_prmtop_path_3} {prep_files_path}/{self.comb_filename}.prmtop"
        )
        return md_path, case_folder_name

    def rdf_analysis(self, md_path: Path, analysis_type, case_folder_name):
        trajectories = {}
        trajin_contents = {}
        parmtops = {}
        run_script = self.config.analysis_script[analysis_type]
        exclude_folders = ["log", "prep_files"]

        traj_folders = list_folders(md_path, exclude_folders)
        sorted_traj_folders = sorted(traj_folders)
        tmp_trajectory, tmp_trajin_content, tmp_parmtop = get_traj_info(
            sorted_traj_folders, md_path
        )

        # append to dict
        trajectories[case_folder_name] = tmp_trajectory
        trajin_contents[case_folder_name] = tmp_trajin_content
        parmtops[case_folder_name] = tmp_parmtop

        analysis_path = f"{self.config.output_dir}/{case_folder_name}/analysis"
        analysis_folder = f"{analysis_path}/{analysis_type}"
        os.makedirs(analysis_folder, exist_ok=True)
        os.chdir(analysis_folder)

        exp_name = f"{analysis_type}_{case_folder_name}"
        cwrk = analysis_path
        output_dir_name = analysis_folder
        parmfile = parmtops[case_folder_name]
        trajin_content = trajin_contents[case_folder_name]
        prepared_trajin_content = trajin_content.replace("\n", "##NEWLINE##")

        analysis_type_ref = f"{self.config.analysis_ref}/{analysis_type}"

        cmd = f'{analysis_type_ref}/{run_script} {exp_name} {cwrk} {output_dir_name} {parmfile} "{prepared_trajin_content}"'

        cwrk_dir_path = Path.cwd()
        os.system(cmd)
        os.chdir(cwrk_dir_path)
        return analysis_folder
