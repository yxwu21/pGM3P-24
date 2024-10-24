import numpy as np

from glob import glob
from typing import Dict
from .md import MDResluts
from .yxwu_lib import read_analysis_data, extract_density, dipole_extract_data


def peak_metrics(peak: np.ndarray, gt: np.ndarray) -> float:
    """A metrics for peak detection.

    :param peak: the peak detected
    :param gt: the ground truth
    :return: a score for the peak detection
    """
    # gt = gt[5::10]
    gt = gt[1::2]
    min_len = min(len(peak), len(gt))
    peak = peak[:min_len]
    gt = gt[:min_len]
    return np.linalg.norm(peak[:, 1] - gt[:, 1])


def eval_md_results(md_results: MDResluts) -> tuple[Dict[str, float], Dict[str, float]]:
    """Evaluate a MD run.

    :param output_dir: the output directory of the MD run
    :return: a dictionary contains all metrics
    """
    # calculate rdf score
    identifier = "#Distance_(Ang)     @O_=>_@O"
    data_file_pattern = f"{md_results.rdf_analysis_folder}/*.out"
    md_rdf, _ = read_analysis_data(
        md_results.rdf_analysis_folder, data_file_pattern, identifier
    )
    exp_rdf, _ = read_analysis_data(
        md_results.experimental_rdf_analysis_folder, data_file_pattern, identifier
    )
    rdf_score = peak_metrics(np.array(md_rdf), np.array(exp_rdf))

    # calculate density score
    density = extract_density(glob(str(md_results.dipole_folder / "*.out"))[0])
    dipole_data = dipole_extract_data(md_results.dipole_folder / "fort.100")
    total_moment = [data[-1] for data in dipole_data]
    avg_total_moment = np.average(total_moment)
    density_score = abs(density - 0.997)
    dipole_score = abs(avg_total_moment - 0.5)

    scores = {
        "rdf_score": rdf_score,
        "density_score": density_score,
        "dipole_score": dipole_score,
    }
    data = {"rdf": np.array(md_rdf), "exp_rdf": np.array(exp_rdf)}
    return scores, data
