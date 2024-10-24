import os
import threading

from pathlib import Path
from .constant_dict import time2nstlim


AMBERENVPREFIX = (
    "export AMBERHOME=/home8/yxwu/AMBER/bin/pmemd-install.dipole; source $AMBERHOME/amber.sh;"
    # "export AMBERHOME=/home/rayl/pmemd-install; source $AMBERHOME/amber.sh;"
)

# TESTsander="/home8/larry/amber24"
# DO_PARALLEL="mpirun -np 4"

# AMBERENVPREFIX = (
#     "export PATH=/home8/larry/apps/openmpi/3.1.6/bin/:$PATH;" \
#     "export LD_LIBRARY_PATH=/home8/larry/apps/openmpi/3.1.6/lib/:$LD_LIBRARY_PATH;" \
#     f"export AMBERHOME={TESTsander}; source $AMBERHOME/amber.sh; {DO_PARALLEL}"
# )


def timeout(time_limit):
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Create a thread to run the function
            result = [Exception("Function did not return")]
            thread = threading.Thread(
                target=lambda: result.append(func(*args, **kwargs))
            )
            thread.start()
            thread.join(time_limit)
            if thread.is_alive():
                # Terminate the function if it's still running after the time limit
                thread._stop()
                raise TimeoutError("Function exceeded the time limit")
            else:
                return result[1]

        return wrapper

    return decorator


class SimulationManager:
    def create_mdin_template(
        self,
        ensemble_type,
        nstlim,
        ntt,
        gamma_ln,
        skinnb,
        dipole_print=False,
        additional_parameters=None,
    ):
        # Common parameters for all ensembles
        common_parameters = {
            "imin": 0,
            "nstlim": nstlim,
            "irest": 1,
            "ntx": 5,
            "ntt": ntt,
            "temp0": 298,
            "gamma_ln": gamma_ln,
            "ntc": 2,
            "ntf": 2,
            "cut": 9.0,
            "ntpr": 1,
            "ntwr": 1,
            "ntwx": 1,
            "dt": 0.001,
            "ipgm": 1,
        }

        # Specific parameters for different ensembles
        ensemble_specific_parameters = {
            "npt": {"ntb": 2, "ntp": 1, "tol": f"{0.000001:.6f}"},
            "nvt": {"ntb": 1, "ntp": 0, "ig": 1, "tol": f"{0.000001:.6f}"},
            "nve": {"ntb": 1, "ntp": 0, "ig": 1, "tol": f"{0.0000001:.7f}"}
            # Add more ensembles as needed
        }

        # Other ensemble-dependent parameters
        ensemble_dependent_params = {
            "dipole_scf_init_step": {"npt": "2", "nvt": "4", "nve": "4"},
            "dipole_scf_init_order": {"npt": "3", "nvt": "2", "nve": "2"},
            "scf_solv_opt": {"npt": "3", "nvt": "3", "nve": "2"},
            "dipole_scf_tol": {"npt": "0.01", "nvt": "0.01", "nve": "0.001"},
        }

        if ensemble_type not in ensemble_specific_parameters:
            raise ValueError(f"Invalid ensemble type: {ensemble_type}")

        # Merge common, specific, and additional parameters
        parameters = {
            **common_parameters,
            **ensemble_specific_parameters[ensemble_type],
            **(additional_parameters or {}),
        }

        # Building the configuration string
        cntrl_section = "\n   ".join(
            [f"{key}={value}," for key, value in parameters.items()]
        )

        # Constructing the pol_gauss section based on dipole_print
        if dipole_print:
            pol_gauss_additional = (
                "dipole_print=1,"  # Replace with actual parameter(s) for dipole_print
            )
        else:
            pol_gauss_additional = ""

        # Build the mdin content
        mdin_content = f"""  single point calc.
 &cntrl
   {cntrl_section}
 /
 &ewald
  skinnb={skinnb}.,ew_coeff=0.4,nfft1=50,nfft2=50,nfft3=50,order=8,vdwmeth=0
 /
 &pol_gauss
   pol_gauss_verbose=0,ee_dsum_cut=9.0,
   dipole_scf_tol={ensemble_dependent_params["dipole_scf_tol"][ensemble_type]},dipole_scf_init=3,dipole_scf_init_order={ensemble_dependent_params["dipole_scf_init_order"][ensemble_type]},dipole_scf_init_step={ensemble_dependent_params["dipole_scf_init_step"][ensemble_type]},
   scf_solv_opt={ensemble_dependent_params["scf_solv_opt"][ensemble_type]},scf_sor_niter=100,scf_sor_coefficient=0.65,
   scf_cg_niter=50,scf_local_niter=3,scf_local_cut=4.0,{pol_gauss_additional}
 /
    """

        return mdin_content

    def construct_md_cmd(
        self,
        md_path: Path,
        output_dir_name: str,
        executable: str,
        file_name: str,
        ensemble_type: str,
        nstlim: int,
        ntt: int,
        gamma_ln: float,
        skinnb: float,
        last_iterate: bool = False,
        last_ensemble: str = None,
        last_executable: str = None,
        last_iter_num: int = None,
        restart: bool = False,
        input_dir_name: str = None,
        dipole_print: bool = False,
    ):
        # Validate inputs
        if not isinstance(md_path, Path):
            raise ValueError("md_path must be a Path object")

        # Executable path
        amber_home = Path("$AMBERHOME")  # Adjust as needed
        exe = amber_home / "bin" / executable

        # File paths
        prefix = f"{file_name}.{ensemble_type}.{executable}"
        parmfile = md_path / "prep_files" / f"{file_name}.prmtop"

        if restart:
            if last_iterate:
                coordfile = (
                    md_path
                    / f"{input_dir_name}"
                    / f"{file_name}.{last_ensemble}.{last_executable}.{last_iter_num}.rst"
                )
            else:
                coordfile = (
                    md_path
                    / f"{input_dir_name}"
                    / f"{file_name}.{last_ensemble}.{last_executable}.rst"
                )
        else:
            coordfile = md_path / "prep_files" / f"{file_name}.restrt0"

        mdin = md_path / output_dir_name / f"{ensemble_type}.in"

        # Generate mdin content
        mdin_content = self.create_mdin_template(
            ensemble_type, nstlim, ntt, gamma_ln, skinnb, dipole_print=dipole_print
        )

        # Write mdin content to file
        with mdin.open("w") as file:
            file.write(mdin_content)

        # Construct the command
        cmd = (
            f"{exe} -O -o {prefix}.out "
            f"-r {prefix}.rst "
            f"-c {coordfile} "
            f"-i {mdin} "
            f"-p {parmfile} "
            f"-ref {coordfile} "
            f"-x {prefix}.nc"
        )

        return cmd

    def construct_iter_md_cmd(
        self,
        md_path: Path,
        output_dir_name: str,
        executable: str,
        file_name: str,
        ensemble_type: str,
        nstlim: int,
        ntt: int,
        gamma_ln: float,
        skinnb: float,
        iteration: int,
        last_iterate: bool = False,
        last_ensemble: str = None,
        last_executable: str = None,
        last_iter_num: int = None,
        restart: bool = False,
        input_dir_name: str = None,
        dipole_print: bool = False,
    ):
        # Validate inputs
        if not isinstance(md_path, Path):
            raise ValueError("md_path must be a Path object")

        # Executable path
        amber_home = Path("$AMBERHOME")  # Adjust as needed
        exe = amber_home / "bin" / executable

        # File paths
        prefix = f"{file_name}.{ensemble_type}.{executable}"
        parmfile = md_path / "prep_files" / f"{file_name}.prmtop"

        mdin = md_path / output_dir_name / f"{ensemble_type}.in"

        if iteration == 1:
            # Generate mdin content
            mdin_content = self.create_mdin_template(
                ensemble_type,
                nstlim,
                ntt,
                gamma_ln,
                skinnb,
                dipole_print=dipole_print,
            )

            # Write mdin content to file
            with mdin.open("w") as file:
                file.write(mdin_content)

            if restart:
                if last_iterate:
                    coordfile = (
                        md_path
                        / f"{input_dir_name}"
                        / f"{file_name}.{last_ensemble}.{last_executable}.{last_iter_num}.rst"
                    )

                else:
                    coordfile = (
                        md_path
                        / f"{input_dir_name}"
                        / f"{file_name}.{last_ensemble}.{last_executable}.rst"
                    )
            else:
                coordfile = md_path / "prep_files" / f"{file_name}.restrt0"
        else:
            coordfile = (
                md_path
                / f"{output_dir_name}"
                / f"{file_name}.{ensemble_type}.{executable}.{iteration - 1}.rst"
            )

        # Construct the command
        iter_cmd = (
            f"{exe} -O -o {prefix}.{iteration}.out "
            f"-r {prefix}.{iteration}.rst "
            f"-c {coordfile} "
            f"-i {mdin} "
            f"-p {parmfile} "
            f"-ref {coordfile} "
            f"-x {prefix}.{iteration}.nc"
        )

        return iter_cmd

    def run_or_restart_simulation(
        self,
        md_path: Path,
        comb_filename: str,
        ensemble_type: str,
        executable: str,
        md_length: str,
        id: str,
        ntt,
        gamma_ln,
        skinnb,
        dipole_print: bool = False,
        iterate: bool = False,
        total_iterations: int = 1,
        last_iterate: bool = False,
        last_ensemble: str = None,
        last_executable: str = None,
        last_iter_num: int = None,
        restart: bool = False,
        input_dir_name: str = None,
        timeout_limit: int = None,
    ) -> str:
        """
        Run or restart a molecular dynamics simulation.

        Parameters:
        - md_path (Path): Path to the molecular dynamics directory.
        - comb_filename (str): Name of the comb file.
        - ensemble_type (str): Type of ensemble.
        - executable (str): Name of the executable to use.
        - md_length (str): Duration of the molecular dynamics simulation.
        - id (str): Identifier for the simulation.
        - ntt, gamma_ln, skinnb: Various simulation parameters.
        - dipole_print (bool): Flag for dipole printing, specific to restarts.
        - iterate (bool): Flag to indicate if it's an iterative process.
        - total_iterations (int): Total number of iterations for restarts.
        - last_iterate, last_ensemble, last_executable, last_iter_num: Parameters for restarts.

        Returns:
        - (str): Name of the output directory.
        """
        nstlim = time2nstlim[md_length]

        if dipole_print:
            output_dir_name = f"dipole_{id}-{md_length}_output_{ensemble_type}"
        elif iterate:
            output_dir_name = (
                f"{id}-{md_length}_{total_iterations}iter_output_{ensemble_type}"
            )
        else:
            output_dir_name = f"{id}-{md_length}_output_{ensemble_type}"

        output_dir_path = md_path / output_dir_name
        output_dir_path.mkdir(exist_ok=True)
        cwrk_dir_path = Path.cwd()

        # timeout command
        timeout_cmd = f"timeout {timeout_limit}h" if timeout_limit else ""

        os.chdir(output_dir_path)
        if iterate:
            for iteration in range(1, total_iterations + 1):
                iter_cmd = self.construct_iter_md_cmd(
                    md_path,
                    output_dir_name,
                    executable,
                    comb_filename,
                    ensemble_type,
                    nstlim,
                    ntt,
                    gamma_ln,
                    skinnb,
                    iteration,
                    last_iterate,
                    last_ensemble,
                    last_executable,
                    last_iter_num,
                    restart,
                    input_dir_name,
                    dipole_print,
                )
                os.system(f"{AMBERENVPREFIX} {timeout_cmd} {iter_cmd}")
        else:
            cmd = self.construct_md_cmd(
                md_path,
                output_dir_name,
                executable,
                comb_filename,
                ensemble_type,
                nstlim,
                ntt,
                gamma_ln,
                skinnb,
                last_iterate,
                last_ensemble,
                last_executable,
                last_iter_num,
                restart,
                input_dir_name,
                dipole_print,
            )
            os.system(f"{AMBERENVPREFIX} {timeout_cmd} {cmd}")

        os.chdir(cwrk_dir_path)
        return output_dir_name
