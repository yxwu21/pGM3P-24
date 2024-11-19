<h1 align="center"><strong>pGM3P-24</strong></h1>

The source code for the paper "Automated Refinement of Property-Specific Polarizable Gaussian Multipole Water Models
Using Bayesian Black Box Optimization".
Our code is implemented using SMAC, and all experiments are executed through the cluster job scheduler slurm.

## pGM3P-24 Water Model Parameters

| Parameter     | pGM3P-24                |
|---------------|----------------------|
| vdW A (kcal/mol * $\AA^{12}$)         | 622716.375799251     |
| vdW B (kcal/mol * $\AA^6$)         | 600.4123381041812    |
| $\sigma_{LJ}$ ($\AA$)         | 3.1815599649844395   |
| $\epsilon_{LJ}$ (kJ/mol)       | 0.1447267928698027   |
| $r_{min}$ ($\AA$)         | 3.5711803151155306   |

| Polarizable Parameter   | O       | H    |
|-------------|--------------------|-------------|
| Charge ($a.u.$)       | -2.040562217277963 | 1.020281087271829 |
| Permenent Dipole moment ($a.u.$)       | -0.3613181603107004 | 0.16228405747594296 |
| Polarizability ($a.u.^3$)     | 7.542285365436364 | 2.2244581584935808 |
| Gaussian radius ($a.u.$)     | 1.1435691447610588 | 1.0133299506847484 |


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
