import os
from pathlib import Path


class BenchmarkRunner:
    def __init__(
        self,
        wat_model,
        wat_num,
        bench_ref,
        bench_prep,
        bench_script,
        temp,
        out_folder_path,
        bench_md_data,
    ):
        self.wat_model = wat_model
        self.wat_num = wat_num
        self.bench_ref = bench_ref
        self.bench_prep = bench_prep
        self.bench_script = bench_script
        self.temp = temp
        self.out_folder_path = out_folder_path
        self.bench_md_data = bench_md_data

    def prepare_directories(self, folder_name):
        file_name = f"{self.wat_model.lower()}_{self.wat_num}wat"
        out_folder = self.out_folder_path / self.wat_model / self.wat_num / folder_name
        md_path = out_folder / "MD"
        md_prep_path = md_path / "prep_files"
        md_prep_path.mkdir(parents=True, exist_ok=True)

        return file_name, out_folder, md_path, md_prep_path

    def submit_job(
        self,
        md_path,
        out_folder,
        output_dir_name,
        cmd,
    ):
        output_dir = md_path / output_dir_name
        output_dir.mkdir(parents=True, exist_ok=True)

        slurm_output = out_folder / "slurm"
        slurm_output.mkdir(parents=True, exist_ok=True)

        os.chdir(slurm_output)
        os.system(cmd)

    def run_minimization(self):
        run_script = self.bench_script / "min.run"
        output_dir_name = "1_Min"
        folder_name = f"{self.wat_model}_{self.wat_num}wat"

        exp_name = f"{self.wat_model}_{self.wat_num}wat"
        self.out_folder_path = bench_prep
        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        init_prmtop = (
            self.bench_ref / self.wat_model / self.wat_num / f"{file_name}.prmtop"
        )
        init_inpcrd = (
            self.bench_ref / self.wat_model / self.wat_num / f"{file_name}.inpcrd"
        )
        os.system(f"cp {init_prmtop} {md_prep_path}/{file_name}.prmtop")
        os.system(f"cp {init_inpcrd} {md_prep_path}/{file_name}.restrt0")

        cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name}"
        self.submit_job(
            md_path,
            out_folder,
            output_dir_name,
            cmd,
        )
        print(f"Simulation for {self.wat_model}_{self.wat_num}wat submitted...")
        print("-" * 160)

    def run_equilibrium(self):
        run_script = self.bench_script / "equi.run"
        output_dir_name = "2_Equi"
        folder_name = f"{self.wat_model}_{self.wat_num}wat-{self.temp}"

        exp_name = f"{self.wat_model}_{self.wat_num}wat-{self.temp}"
        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        init_prmtop = (
            self.bench_ref / self.wat_model / self.wat_num / f"{file_name}.prmtop"
        )
        init_inpcrd = (
            self.bench_md_data
            / "prep"
            / self.wat_model
            / self.wat_num
            / f"{self.wat_model}_{wat_num}wat"
            / "MD"
            / "1_Min"
            / f"{self.wat_model.lower()}_{self.wat_num}wat.min.pmemd.rst"
        )
        os.system(f"cp {init_prmtop} {md_prep_path}/{file_name}.prmtop")
        os.system(f"cp {init_inpcrd} {md_prep_path}/{file_name}.restrt0")

        cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp}"
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

    def run_langevin_production(self):
        run_script = self.bench_script / "prod.langevin.run"
        output_dir_name = "3_Langevin_Prod"
        folder_name = f"{self.wat_model}_{self.wat_num}wat-{self.temp}"
        min_iter = "1"
        max_iter = "10"

        exp_name = f"{self.wat_model}_{self.wat_num}wat-{self.temp}"
        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {min_iter} {max_iter}"
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

    def run_production(self):
        run_script = self.bench_script / "prod.run"
        output_dir_name = "4_Prod"
        folder_name = f"{self.wat_model}_{self.wat_num}wat-{self.temp}"
        min_iter = "1"
        max_iter = "10"

        exp_name = f"{self.wat_model}_{self.wat_num}wat-{self.temp}"
        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name} {self.temp} {min_iter} {max_iter}"
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


# Usage
bench_md_data = Path("/home/yxwu/pGM_water_model/MD_data_Benchmark")
bench_prep = bench_md_data / "prep"
bench_ref = Path("/home/yxwu/pGM_water_model/datasets/benchmark")
bench_script = Path("/home/yxwu/pGM_water_model/scripts/md/benchmark")
wat_models = ["TIP3P", "TIP4P", "TIP5P"]

temperatures = set(range(240, 381, 5))

temps = [
    "298",
] + list(temperatures)

wat_num = "512"


for wat_model in wat_models:
    for temp in temps:
        out_folder_path = bench_md_data
        # out_folder_path = Path("/home/yxwu/pGM_water_model/md_test")
        out_folder_path.mkdir(parents=True, exist_ok=True)
        runner = BenchmarkRunner(
            wat_model,
            wat_num,
            bench_ref,
            bench_prep,
            bench_script,
            temp,
            out_folder_path,
            bench_md_data,
        )
        # runner.run_minimization()
        runner.run_equilibrium()
        # runner.run_langevin_production()
