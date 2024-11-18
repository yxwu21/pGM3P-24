import os
from pathlib import Path
import shutil
import logging


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
        self.out_folder_path = Path(out_folder_path)
        self.bench_md_data = bench_md_data
        logging.basicConfig(level=logging.INFO)

    def prepare_directories(self, folder_name):
        file_name = f"{self.wat_model.lower()}_{self.wat_num}wat"
        out_folder = (
            self.out_folder_path / self.wat_model / str(self.wat_num) / folder_name
        )
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

    def copy_files(self, source, destination):
        try:
            shutil.copy(source, destination)
        except IOError as e:
            logging.error(f"Unable to copy file. {e}")
        except:
            logging.error("Unexpected error:", sys.exc_info())

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

    def run_simulation(
        self, stage, script_name, copy_initial_files=True, min_iter=None, max_iter=None
    ):
        run_script = self.bench_script / script_name
        output_dir_name = stage
        folder_name = f"{self.wat_model}_{self.wat_num}wat-{self.temp if 'Min' not in stage else ''}"
        exp_name = folder_name

        file_name, out_folder, md_path, md_prep_path = self.prepare_directories(
            folder_name
        )

        if copy_initial_files:
            folder_name = f"{self.wat_model}_{self.wat_num}wat"
            init_prmtop = (
                self.bench_ref
                / self.wat_model
                / str(self.wat_num)
                / f"{file_name}.prmtop"
            )
            init_inpcrd = (
                self.bench_ref
                / self.wat_model
                / str(self.wat_num)
                / f"{file_name}.inpcrd"
                if stage == "1_Min"
                else self.bench_md_data
                / "prep"
                / self.wat_model
                / str(self.wat_num)
                / folder_name
                / "MD"
                / "1_Min"
                / f"{file_name}.min.pmemd.rst"
            )
            self.copy_files(init_prmtop, md_prep_path / f"{file_name}.prmtop")
            self.copy_files(init_inpcrd, md_prep_path / f"{file_name}.restrt0")

        cmd = f"sbatch -J {exp_name} {run_script} {exp_name} {file_name} {md_path} {output_dir_name}"
        if "Min" not in stage:
            cmd += f" {temp}"
            if min_iter is not None and max_iter is not None:
                cmd += f" {min_iter} {max_iter}"

        self.submit_job(md_path, out_folder, output_dir_name, cmd)
        logging.info(f"Simulation for {exp_name} submitted...")

    def execute_stages(self, stages):
        if stages.get("minimization", False):
            self.run_minimization()
        if stages.get("equilibrium", False):
            self.run_simulation("2_Equi", "equi.run", copy_initial_files=True)
        if stages.get("langevin_production", False):
            self.run_simulation(
                "3_Langevin_Prod",
                "prod.langevin.run",
                copy_initial_files=False,
                min_iter="1",
                max_iter="10",
            )
        if stages.get("production", False):
            self.run_simulation(
                "4_Prod",
                "prod.run",
                copy_initial_files=False,
                min_iter="1",
                max_iter="10",
            )


if __name__ == "__main__":

    # bench_md_data = Path("/home8/yxwu/pGM_water_model/MD_data_Benchmark")
    # bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark")

    # bench_md_data = Path("/home8/yxwu/pGM_water_model/MD_data_Benchmark-MC")
    # bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark-MC")

    # bench_md_data = Path("/home8/yxwu/pGM_water_model/MD_data_Benchmark-MC-Langevin")
    # bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark-MC-Langevin")

    # bench_md_data = Path("/home8/yxwu/pGM_water_model/MD_data_Benchmark-gpu")
    # bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark-gpu")

    # bench_md_data = Path(
    #     "/home8/yxwu/pGM_water_model/MD_data_Benchmark-MC-Langevin-gamma1"
    # )
    # bench_script = Path(
    #     "/home8/yxwu/pGM_water_model/scripts/md/benchmark-MC-Langevin-gamma1"
    # )

    # bench_md_data = Path(
    #     "/home8/yxwu/pGM_water_model/MD_data_Benchmark-MC-Langevin-gamma0.1"
    # )
    # bench_script = Path(
    #     "/home8/yxwu/pGM_water_model/scripts/md/benchmark-MC-Langevin-gamma0.1"
    # )

    # bench_md_data = Path("/home8/yxwu/pGM_water_model/MD_data_Benchmark-cut6")
    # bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark-cut6")

    # bench_md_data = Path("/home8/yxwu/pGM_water_model/MD_data_benchmark-Langevin-MC")
    # bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark-Langevin-MC")

    # bench_md_data = Path("/home8/yxwu/pGM_water_model/MD_data_Benchmark-gpu")
    # bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark-gpu")

    bench_md_data = Path("/home8/yxwu/pGM_water_model/MD_data_Benchmark-long")
    bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark-long")

    bench_prep = bench_md_data / "prep"
    bench_ref = Path("/home8/yxwu/pGM_water_model/datasets/benchmark")
    # bench_script = Path("/home8/yxwu/pGM_water_model/scripts/md/benchmark")

    wat_models = [
        "TIP3P",
        # "TIP4P",
        # "TIP5P",
        # "OPC",
        # "OPC3",
        # "SPCE",
        # "TIP4PEW",
    ]

    # temps = ["298"] + [str(t) for t in range(240, 381, 5)]
    # temps = [str(t) for t in range(240, 381, 5)]
    # temps = ["298"]
    temps = ["298", "245", "265", "330", "335"]

    wat_num = "512"

    # testing = False
    testing = True

    stages_to_run = {
        "minimization": True,
        # "equilibrium": True,
        # "langevin_production": True,
        # "production": True,
    }

    out_folder_path = (
        Path("/home8/yxwu/pGM_water_model/md_test") if testing else bench_md_data
    )
    out_folder_path.mkdir(parents=True, exist_ok=True)

    if stages_to_run.get("minimization", False):
        temps = ["298"]

    for wat_model in wat_models:
        for temp in temps:
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
            runner.execute_stages(stages_to_run)
            # runner.run_simulation("1_Min", "min.run", copy_initial_files=True)
            # runner.run_simulation("2_Equi", "equi.run", copy_initial_files=True)
            # runner.run_simulation(
            #     "3_Langevin_Prod",
            #     "prod.langevin.run",
            #     copy_initial_files=False,
            #     min_iter="1",
            #     max_iter="10",
            # )
            # runner.run_simulation(
            #     "4_Prod",
            #     "prod.run",
            #     copy_initial_files=False,
            #     min_iter="1",
            #     max_iter="10",
            # )

            if testing:
                break
        if testing:
            break
