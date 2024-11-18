import os
import glob
import tyro

from pathlib import Path
from dataclasses import dataclass, field
from typing import Dict
from mdtool.slurm import SlurmArgs, slurm_launcher


@dataclass
class MDRunConfig:
    analysis_script: Dict[str, str] = field(
        default_factory=lambda: {
            "rdf_goo": "cpptraj_rdf_goo_MPI.sh",
            "rdf_goh": "cpptraj_rdf_goh_MPI.sh",
            "rdf_ghh": "cpptraj_rdf_ghh_MPI.sh",
            "diffusion": "cpptraj_diffusion_MPI.sh",
        }
    )

    analysis_ref: str = "/home8/yxwu/pGM_water_model/scripts/analysis"

    experimental_rdf_analysis_folder: str = (
        "datasets/pgm_512wat_main/expt_Skinner_2013/analysis/rdf_goo"
    )


class TrajectoryAnalysis:
    def __init__(self, config: MDRunConfig):
        self.config = config

    def list_folders(self, path, exclude_folders):
        try:
            all_entries = os.listdir(path)
            folders = [
                entry
                for entry in all_entries
                if os.path.isdir(os.path.join(path, entry))
                and entry not in exclude_folders
            ]
            return folders

        except FileNotFoundError:
            print(f"The system cannot find the specified path: {path}")
        except PermissionError:
            print(f"Permission denied for the specified path: {path}")
        except Exception as e:
            print(f"An error occurred: {e}")
        return []

    def get_trajectory_info_custom_range(
        self, traj_folders, md_folder_path, last_n_trajectories
    ):
        trajectory = []
        parm_pattern = f"{md_folder_path}/prep_files/*.prmtop"
        parm_path = glob.glob(parm_pattern)[0]
        parmtop = parm_path.split("/")[-1]
        parmtop = f"../../MD/prep_files/{parmtop}"

        for traj_folder in traj_folders:
            traj_path_pattern = f"{md_folder_path}/{traj_folder}/*.nc"
            traj_paths = sorted(glob.glob(traj_path_pattern))
            if traj_paths:
                for traj_path in traj_paths:
                    trajectory.append(traj_path)

        def safe_extract_number(filename):
            parts = filename.split(".")
            if parts[-2].isdigit():
                return int(parts[-2])
            return 0

        trajin_lists = []
        for i in sorted(
            trajectory, key=lambda x: (x.split("/")[-2], safe_extract_number(x))
        ):
            traj_case = i.split("/")[-2]
            traj_nc = i.split("/")[-1]
            relative_traj_nc = f"trajin ../../MD/{traj_case}/{traj_nc}"
            trajin_lists.append(relative_traj_nc)

        start = len(trajin_lists) - last_n_trajectories
        end = len(trajin_lists)
        trajin_content = "\n".join(trajin_lists[start:end])
        print(trajin_lists)
        print(len(trajin_lists))
        return trajectory, trajin_content, parmtop

    def traj_analysis(
        self, md_path: Path, analysis_type: str, last_n_trajectories: int = 10
    ):
        trajectories = {}
        trajin_contents = {}
        parmtops = {}
        run_script = self.config.analysis_script[analysis_type]
        exclude_folders = ["log", "prep_files"]

        traj_folders = self.list_folders(md_path, exclude_folders)
        sorted_traj_folders = sorted(traj_folders)
        (
            tmp_trajectory,
            tmp_trajin_content,
            tmp_parmtop,
        ) = self.get_trajectory_info_custom_range(
            sorted_traj_folders, md_path, last_n_trajectories
        )

        md_path = Path(md_path)
        case_folder_name = md_path.parent.name
        trajectories[case_folder_name] = tmp_trajectory
        trajin_contents[case_folder_name] = tmp_trajin_content
        parmtops[case_folder_name] = tmp_parmtop

        analysis_path = md_path.parent / "analysis"

        analysis_folder = f"{analysis_path}/{analysis_type}-last_{last_n_trajectories}"
        # if analysis_type == "diffusion":
        #     analysis_folder = (
        #         f"{analysis_path}/{analysis_type}-last_{last_n_trajectories}"
        #     )
        # else:
        #     analysis_folder = f"{analysis_path}/{analysis_type}"
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


@dataclass
class Args:
    slurm: SlurmArgs

    md_path: str = (
        "/home8/yxwu/pGM_water_model/MD_data/langevin/vital-cosmos-120_144/MD"
    )
    analysis_type: str = "diffusion"
    last_n_trajectories: int = 10


@slurm_launcher(Args)
def main(args: Args):
    config = MDRunConfig()
    analysis = TrajectoryAnalysis(config)

    analysis.traj_analysis(args.md_path, args.analysis_type, args.last_n_trajectories)


if __name__ == "__main__":
    args = tyro.cli(Args)
    main()
