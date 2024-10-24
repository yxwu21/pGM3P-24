THERMOSTAT = "langevin"
BAROSTAT = "Berendsen"
# exe = "pmemd"
# BAROSTAT = "constTemp"
# PIPELINE_NAME = "nvt-npt-npt-dipole"
PIPELINE_NAME = "nvt-npt-dipole"

pipeline_dict = {
    "langevin_Berendsen": {
        "nvt-npt-dipole": {
            "MD_LENGTHS": "10ps 200ps 1step",
            "GAMMA_LNS": "100 100 100",
        }
    },
    "langevin_constTemp": {
        "nvt-npt-npt-dipole": {
            "MD_LENGTHS": "10step 10step 10step 1step",
            "GAMMA_LNS": "100 100 0 100",
        }
    },
}
