<h1 align="center"><strong>pGM-Tri24</strong></h1>

The source code for the paper "Automated Refinement of Property-Specific Polarizable Gaussian Multipole Water Models
Using Bayesian Black Box Optimization".
Our code is implemented using SMAC, and all experiments are executed through the cluster job scheduler slurm.

| Parameters   | pGM3P-24              |
|-------------|--------------------|
| vdw_a       | 622716.375799251  |
| vdw_b       | 600.4123381041812 |
| scale_q     | 1.0683576293005586 |
| scale_p     | 0.7762292363437578 |
| scale_pol   | 0.7713367864674853 |
| scale_rad   | 0.7502257723289764 |
| sigma       | 3.1815599649844395 |
| epsilon     | 0.1447267928698027 |
| r_min       | 3.5711803151155306 |

## File Structures

- `configs`: contains all experiment related configurations
- `scripts`: contains all MD scripts and analysis scripts
- `src`: contains all source code for running automated machine learning searching
- `run.py`: run automated machine learning based searching
- `./*.py`: figure plots

To view all configurable options, run

```bash
python run.py -h
```

## Dependencies

1. Python: 3.11

2. Install other packages

```bash
pip install -r requirements.txt
```

## Acknowledgements

If you have any questions, feel free to reach out to us via email at `yongxian.wu@uci.edu`.
