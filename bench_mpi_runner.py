import os
from pathlib import Path
import shutil
import logging


class pGMRunner:
    def __init__(
        self,
        wat_model,
        wat_num,
        run_script,
        temp,
        out_folder_path,
        run_md_data,
        md_data_ref,
        case,
        rst_folder,
        slurm_run,
    ):
        self.wat_model = wat_model
        self.wat_num = wat_num
        self.run_script = run_script
        self.temp = temp
        self.out_folder_path = Path(out_folder_path)
        self.run_md_data = run_md_data
        self.md_data_ref = md_data_ref
        self.case = case
        self.rst_folder = rst_folder
        self.slurm_run = slurm_run
        logging.basicConfig(level=logging.INFO)

    def prepare_directories(self, folder_name):
        file_name = f"{self.wat_model.lower()}_{self.wat_num}wat"
        out_folder = self.out_folder_path / self.case / str(self.wat_num) / folder_name
        md_path = out_folder / "MD"
        md_prep_path = md_path / "prep_files"
        md_prep_path.mkdir(parents=True, exist_ok=True)

        return file_name, out_folder, md_path, md_prep_path

    def submit_job(self, md_path, out_folder, output_dir_name, cmd):
        output_dir = md_path / output_dir_name
        output_dir.mkdir(parents=True, exist_ok=True)

        slurm_output = out_folder / "slurm"
        slurm_output.mkdir(parents=True, exist_ok=True)

        os.chdir(slurm_output)
        os.system(cmd)

    def run_equilibrium(self):
        run_script = self.run_script / "equi1.run"
        output_dir_name = "1_Equi"
        folder_name = f"{self.case}-{self.temp}"

        exp_name = f"{self.case}-{self.temp}"
        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        case_name = f"{self.case}_{self.wat_num}wat-{self.temp}"

        init_prmtop = (
            self.md_data_ref
            / case
            / str(self.wat_num)
            / case_name
            / "MD"
            / "prep_files"
            / f"{file_name}.prmtop"
        )

        init_inpcrd = (
            self.md_data_ref
            / case
            / str(self.wat_num)
            / case_name
            / "MD"
            / self.rst_folder
            / f"{file_name}.npt.pmemd.10.rst"
        )

        os.system(f"cp {init_prmtop} {md_prep_path}/{file_name}.prmtop")
        os.system(f"cp {init_inpcrd} {md_prep_path}/{file_name}.restrt0")

        coordfile = f"../prep_files/{file_name}.restrt0"

        if self.slurm_run:
            cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {coordfile}"
        else:
            cmd = f"{run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {coordfile}"
        self.submit_job(
            md_path,
            out_folder,
            output_dir_name,
            cmd,
        )
        print(
            f"Simulation for {self.wat_model}_{self.wat_num}wat-{self.temp} submitted..."
        )
        print("-" * 160)

    def run_equilibrium_before_prod(self):
        run_script = self.run_script / "equi2.run"
        output_dir_name = "2_Equi"
        folder_name = f"{self.case}-{self.temp}"

        exp_name = f"{self.case}-{self.temp}"
        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        coordfile = f"../1_Equi/{file_name}.npt.pmemd.rst"

        if self.slurm_run:
            cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {coordfile}"
        else:
            cmd = f"{run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {coordfile}"
        self.submit_job(
            md_path,
            out_folder,
            output_dir_name,
            cmd,
        )
        print(
            f"Simulation for {self.wat_model}_{self.wat_num}wat-{self.temp} submitted..."
        )
        print("-" * 160)

    def run_prod(self):
        run_script = self.run_script / "prod.run"
        output_dir_name = "3_Prod"
        folder_name = f"{self.case}-{self.temp}"

        exp_name = f"{self.case}-{self.temp}"
        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        coordfile = f"../2_Equi/{file_name}.npt.pmemd.rst"

        min_iter = "1"
        max_iter = "10"

        if self.slurm_run:
            cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {min_iter} {max_iter} {coordfile}"
        else:
            cmd = f"{run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {min_iter} {max_iter} {coordfile}"
        self.submit_job(
            md_path,
            out_folder,
            output_dir_name,
            cmd,
        )
        print(
            f"Simulation for {self.wat_model}_{self.wat_num}wat-{self.temp} submitted..."
        )
        print("-" * 160)

    def run_dipole(self):
        run_script = self.run_script / "prod.run"
        output_dir_name = "3_Prod"
        folder_name = f"{self.case}-{self.temp}"

        exp_name = f"{self.case}-{self.temp}"
        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        coordfile = f"../2_Equi/{file_name}.npt.pmemd.rst"

        min_iter = "1"
        max_iter = "10"

        cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {min_iter} {max_iter} {coordfile}"
        self.submit_job(
            md_path,
            out_folder,
            output_dir_name,
            cmd,
        )
        print(
            f"Simulation for {self.wat_model}_{self.wat_num}wat-{self.temp} submitted..."
        )
        print("-" * 160)

    def execute_stages(self, stages):
        if stages.get("1_Equi", False):
            self.run_equilibrium()
        if stages.get("2_Equi", False):
            self.run_equilibrium_before_prod()
        if stages.get("3_Prod", False):
            self.run_prod()


if __name__ == "__main__":
    proj_path = Path("/home8/yxwu/pGM_water_model")

    # run_md_data = proj_path / "MD_data-pGM-temp-local"
    # run_script = proj_path / "scripts" / "md" / "benchmark-pgm_mpi"

    # run_md_data = proj_path / "MD_data_Benchamrk-pgm_mpi-test"
    # run_script = proj_path / "scripts" / "md" / "benchmark-pgm_mpi-test"

    # run_md_data = proj_path / "MD_data_Benchamrk-pgm_mpi-3"
    run_md_data = proj_path / "MD_data_Benchamrk-long-pgm_mpi"
    run_script = proj_path / "scripts" / "md" / "benchmark-pgm_mpi"

    # run_md_data = proj_path / "MD_data_Benchamrk-pgm_mpi-gpu"
    # run_script = proj_path / "scripts" / "md" / "benchmark-pgm_mpi-gpu"

    # run_md_data = proj_path / "MD_data_Benchamrk-pgm_mpi-gpu-test"
    # run_script = proj_path / "scripts" / "md" / "benchmark-pgm_mpi-gpu-test"

    # md_data_ref_folder = "MD_data_Benchmark"
    md_data_ref_folder = "MD_data_Benchmark-long"
    rst_folder = "4_Prod"

    wat_models = [
        "TIP3P",
        "TIP4P",
        "TIP5P",
        "OPC",
        "OPC3",
        "SPCE",
        "TIP4PEW",
    ]

    md_data_ref = proj_path / md_data_ref_folder

    temps = ["298"] + [str(t) for t in range(240, 381, 5)]
    # temps = [str(t) for t in range(240, 381, 5)]
    # temps = ["298"]

    wat_num = "512"

    testing = False
    # testing = True

    slurm_run = True

    stages_to_run = {
        "1_Equi": True,
        # "2_Equi": True,
        # "3_Prod": True,
    }

    out_folder_path = (
        Path(f"/home8/yxwu/pGM_water_model/md_test/{run_md_data.name}")
        if testing
        else run_md_data
    )
    out_folder_path.mkdir(parents=True, exist_ok=True)

    for wat_model in wat_models:
        case = wat_model
        for temp in temps:
            runner = pGMRunner(
                wat_model,
                wat_num,
                run_script,
                temp,
                out_folder_path,
                run_md_data,
                md_data_ref,
                case,
                rst_folder,
                slurm_run,
            )
            runner.execute_stages(stages_to_run)

            if testing:
                break
        if testing:
            break
