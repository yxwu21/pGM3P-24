<h1 align="center"><strong>pGM-Tri24</strong></h1>

The source code for the paper "Automated Refinement of Property-Specific Polarizable Gaussian Multipole Water Models
Using Bayesian Black Box Optimization". Our code is implemented using SMAC and all experiments are run through the
cluster job scheduler `slurm`.

## File Structures

- `configs`: contains all experiment related configurations
- `scripts`: contains all MD scripts and analysis scripts
- `src`: contains all source code for running automated machine learning searching
- `run.py`: run automated machine learning based searching
- `./*.py`: figure plots

To see all configurable options, run

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
