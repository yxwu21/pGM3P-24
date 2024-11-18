#!/bin/python
# -*- coding:utf-8 -*-

import numpy as np
import pandas as pd
import glob

Epsilon0 = 8.854187817 * 10**-12  # in the unit of F/m
D2Cm = 3.335 * 10**-30  # in the unit of C \cdot m
nm2m = 1 * 10**-9  # in the unit of m
A2nm = 10**-1  # in the unit of nm
kb = 1.380649 * 10**-23  # in the unit of J/K
eA2D = 1 / 0.2081943  # in the unit of Debye


def permittivity(M_squared_mean, M_mean_squared, T, V):
    """
    T: in the unit of K
    V: in the unit of nm^3
    M_squared_mean/ M_mean_squared in the unit of Debye
    """
    epsilon = (
        1
        + 1
        / 3
        / (kb * T * V * nm2m**3)
        / Epsilon0
        * (M_squared_mean - M_mean_squared)
        * D2Cm**2
    )
    return epsilon


def load_fort(filename, numW):
    dipoles = []
    output = {}
    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            else:
                if "@ Step" in line:
                    frame = int(line.split()[-1])
                    output[frame] = {}
                    output[frame]["x"] = []
                    output[frame]["y"] = []
                    output[frame]["z"] = []
                    output[frame]["t"] = []
                    for x in range(numW):
                        line = f.readline()
                        output[frame]["x"].append(float(line.split()[-4]))
                        output[frame]["y"].append(float(line.split()[-3]))
                        output[frame]["z"].append(float(line.split()[-2]))
                        output[frame]["t"].append(float(line.split()[-1]))
                    dipole = pd.DataFrame(output[frame])[["x", "y", "z"]].values
                    dipoles.append(np.sum(dipole, axis=0))
    return np.array(dipoles)


def load_out(filename):
    output = {"T": [], "V": []}
    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            else:
                if "A V E R A G E S   O V E R" in line:
                    break
                elif "VOLUME" in line:
                    output["V"].append(float(line.split()[-1]))
                elif "TEMP(K)" in line:
                    output["T"].append(float(line.split()[-4]))
    return output


if __name__ == "__main__":
    case_path_pattern = "/home8/yxwu/pGM_water_model/MD_data/langevin_Berendsen-295/nvt-npt-dipole/vital-cosmos-120_340/MD/dipole_*-200ps_output_npt"
    case_paths = glob.glob(case_path_pattern)
    for case_path in case_paths[:-1]:
        print(f"{case_path}:")
        filename = f"{case_path}/fort.200"
        dipoles = load_fort(filename, 512)
        print("Dipoles in the unit of e \cdot A")
        print(dipoles)
        print("Mean dipole moment in x, y, z (eA)")
        print(np.mean(dipoles, axis=0))
        print("squared Mean (Debye^2)")
        print(np.mean(np.sum(dipoles**2, axis=1)) * eA2D**2)
        print("Mean squared (Debye^2)")
        print(np.sum(np.mean(dipoles, axis=0) ** 2) * eA2D**2)

        filelog = f"{case_path}/case_1_pgm_512wat.npt.sander.out"
        outT_V = load_out(filelog)
        T = np.array(outT_V["T"])
        V = np.array(outT_V["V"])
        print("Temperature in the unit of K")
        print(T)
        print("Volume in the unit of A^3")
        print(V)

        if dipoles.shape[0] != len(T):
            print("ERROR! The length of dipoles is inconsistent with Temp")
        else:
            for i in range(dipoles.shape[0]):
                accumulated_dxyz = dipoles[: i + 1, :]
                accumulated_T = T[: i + 1]
                accumulated_V = V[: i + 1]

                M_squared_mean = np.mean(np.sum(accumulated_dxyz**2, axis=1))
                M_mean_squared = np.sum(np.mean(accumulated_dxyz, axis=0) ** 2)

                T_mean = np.mean(accumulated_T)
                V_mean = np.mean(accumulated_V)

                epsilon = permittivity(
                    M_squared_mean * eA2D**2,
                    M_mean_squared * eA2D**2,
                    T_mean,
                    V_mean * A2nm**3,
                )
                print(i, M_squared_mean * eA2D**2, M_mean_squared * eA2D**2, epsilon)

        print("-" * 180)
