#!/usr/bin/env python3
import argparse as ap

# import numpy as np
import sys

# import f90nml
# import math
# from scipy.linalg import lu_factor, lu_solve

###### parameters ######
maxq = 8000
natm = 0

###### constants ######
au2A = 0.52917725  # constant converting from a.u. (bohr) to angstrom
au2D = 2.541765  # constant converting from a.u. (e*bohr) to debye

###### files ######
qin, prmtop_in, prmtop_out = [None] * 3
qin_file, prmtop_in_file, prmtop_out_file = [None] * 3

###### runlab ######
title = None

###### orig ######
crd, q, p, pol, rad = [None] * 5  # coordinates, charge, dipole, polarizability, radii


def file_in():
    parser = ap.ArgumentParser(
        usage="amend_pgm_wat.py [-O] -q qin -p prmtop -o prmtop_out"
    )
    parser.add_argument(
        "-O", action="store_true", help="Overwrite output files if they exist"
    )
    parser.add_argument(
        "-q", "--qin", help="type: input, needed; description: pGM qout file"
    )
    parser.add_argument(
        "-p",
        "--prmtop",
        required=True,
        help="type: input, needed; description: input prmtop file",
    )
    parser.add_argument(
        "-o",
        "--prmtop_out",
        required=True,
        help="type: output, always produced; description: output of results",
    )

    args = parser.parse_args()
    qin, prmtop_in, prmtop_out = args.qin, args.prmtop, args.prmtop_out

    return qin, prmtop_in, prmtop_out


def end_sec(readfile):
    line = readfile.readline()
    while line != "\n":
        line = readfile.readline()


def next_sec(readfile):
    line = {}

    line = readfile.readline()
    while line == "\n":
        line = readfile.readline()
    line = line.split()
    while line[0] != "%FLAG":
        line = readfile.readline().split()

    return line


def read_qin():
    global natm, ndip
    global crd, q, p, pol, rad

    crd = {}
    q = {}
    p = {}
    pol = {}
    rad = {}
    line = {}

    qin_file = open(qin, "r")

    read_crd = 0
    read_q = 0
    read_p = 0

    while read_crd == 0 or read_q == 0 or read_p == 0:
        line = next_sec(qin_file)
        if line[1] == "ATOM":
            if line[2] == "CRD":
                line = qin_file.readline()
                line = qin_file.readline().split()
                while line:
                    natm = int(line[0]) - 1
                    crd[natm] = [float(line[1]), float(line[2]), float(line[3])]
                    pol[natm] = float(line[4])
                    rad[natm] = float(line[5])
                    line = qin_file.readline().split()
                read_crd = natm + 1

            elif line[2] == "CHRG:":
                line = qin_file.readline()
                line = qin_file.readline().split()
                while line:
                    ia = int(line[0]) - 1
                    q[ia] = float(line[3])
                    line = qin_file.readline().split()
                read_q = ia + 1

        elif line[1] == "PERM":
            if line[2] == "DIP" and line[3] == "LOCAL:":
                line = qin_file.readline()
                line = qin_file.readline().split()
                while line:
                    ndip = int(line[0]) - 1
                    p[ndip] = [int(line[1]), int(line[2]), float(line[4])]
                    line = qin_file.readline().split()
                read_p = ndip + 1

        else:
            end_sec(qin_file)

    return natm + 1, ndip + 1, crd, q, p, pol, rad


def wrt_int_blk(ofile, nwat, nitm, lst, shift):
    m = 0
    for i in range(nwat):
        for j in range(nitm):
            k = i * shift
            ofile.write("{:8d}".format(lst[j] + k))
            m = m + 1
            if m == 10:
                ofile.write("\n")
                m = 0
    if m != 0:
        ofile.write("\n")
        m = 0


def wrt_flt_blk(ofile, nwat, nitm, lst):
    m = 0
    for i in range(nwat):
        for j in range(nitm):
            ofile.write("{:16.8E}".format(lst[j]))
            # ofile.flush()
            m = m + 1
            if m == 5:
                ofile.write("\n")
                m = 0
    if m != 0:
        ofile.write("\n")
        m = 0


def wrt_out(natm, ndip, prmin, prmout):
    toks = {}
    line = prmin.readline()
    k = 0
    m = 0
    l = 0

    # copy the input prmtop

    while line:
        toks = line.split()
        if line != "\n":
            if toks[0] == "%FLAG":
                if toks[1] == "IPOL":
                    l = -1
                else:
                    l = 0
        if l == 0:
            prmout.write(line)
        if k >= 0:
            toks = line.split()
            if toks[0] == "%FLAG":
                if toks[1] == "POINTERS":
                    k = 1
                else:
                    k = -k
            if k == 1:
                m = m + 1
                if m == 3:
                    natm = int(toks[0])
                elif m == 4:
                    nwat = int(toks[1])
        line = prmin.readline()
    natmpwat = int(natm / nwat)

    if natmpwat != 3:
        print("ERROR: number of atoms per water not equals to 3.\n")
        sys.exit()

    prmout.write("%FLAG IPOL\n")
    prmout.write("%FORMAT(1I8)\n")
    prmout.write("%8d\n" % (0))

    prmout.write("%FLAG POL_GAUSS_FORCEFIELD\n")
    prmout.write(
        "%COMMENT This indicates that this parm file is specific to pGM force field\n"
    )
    prmout.write("%COMMENT This must be present if ipgm (in mdin) is 1\n")
    prmout.write("%COMMENT This must NOT be present if ipgm is 0\n")
    prmout.write("%FORMAT(i5)\n")
    prmout.write("%5d\n" % (1))

    prmout.write("%FLAG POL_GAUSS_COVALENT_POINTERS_LIST\n")
    prmout.write("%COMMENT number of covalent dipoles per atom\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * natmpwat, nwat, natmpwat))
    prmout.write("%FORMAT(10I8)\n")
    nptr = {}
    for i in range(natm):
        nptr[i] = 0
    for i in range(nwat):
        for j in range(ndip):
            if p[j][0] <= 3:
                idx = p[j][0] - 1 + i * 3
                # 				print(i, j, idx)
                nptr[idx] += 1

    wrt_int_blk(prmout, nwat, natmpwat, nptr, 0)

    prmout.write("%FLAG POL_GAUSS_COVALENT_ATOMS_LIST\n")
    prmout.write("%COMMENT   dimension = ")
    m = 0
    for j in range(natmpwat):
        m = m + nptr[j]
    prmout.write("(%d=%d*%d)\n" % (nwat * m, nwat, m))
    prmout.write("%FORMAT(10I8)\n")
    lst = {}
    for j in range(m):
        lst[j] = p[j][1]
    wrt_int_blk(prmout, nwat, m, lst, natmpwat)

    prmout.write("%FLAG POL_GAUSS_COVALENT_DIPOLES_LIST\n")
    prmout.write("%COMMENT   unit: e-Angstrom\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * m, nwat, m))
    prmout.write("%FORMAT(5E16.8)\n")
    for j in range(m):
        lst[j] = p[j][2] * au2A
    wrt_flt_blk(prmout, nwat, m, lst)

    prmout.write("%FLAG POL_GAUSS_MONOPOLES_LIST\n")
    prmout.write("%COMMENT   unit: e\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * natmpwat, nwat, natmpwat))
    prmout.write("%FORMAT(5E16.8)\n")
    wrt_flt_blk(prmout, nwat, natmpwat, q)

    prmout.write("%FLAG POL_GAUSS_RADII_LIST\n")
    prmout.write("%COMMENT   unit: Angstrom\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * natmpwat, nwat, natmpwat))
    prmout.write("%FORMAT(5E16.8)\n")

    for j in range(natmpwat):
        rad[j] = rad[j] * au2A
    wrt_flt_blk(prmout, nwat, natmpwat, rad)

    prmout.write("%FLAG POL_GAUSS_POLARIZABILITY_LIST\n")
    prmout.write("%COMMENT   unit: Angstrom**3\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * natmpwat, nwat, natmpwat))
    prmout.write("%FORMAT(5E16.8)\n")
    for j in range(natmpwat):
        pol[j] = pol[j] * au2A**3
    wrt_flt_blk(prmout, nwat, natmpwat, pol)


def wrt_out_with_replace(natm, ndip, prmin, prmout):
    toks = {}
    line = prmin.readline()
    k = 0
    m = 0
    l = 0

    # copy the input prmtop

    while line:
        toks = line.split()

        # check reach the modify section
        if line != "\n" and l != -1:
            if toks[0] == "%FLAG":
                if toks[1] == "IPOL":
                    l = -1
                else:
                    l = 0
        if l == 0:
            prmout.write(line)

        if k >= 0:
            toks = line.split()
            if toks[0] == "%FLAG":
                if toks[1] == "POINTERS":
                    k = 1
                else:
                    k = -k
            if k == 1:
                m = m + 1
                if m == 3:
                    natm = int(toks[0])
                elif m == 4:
                    nwat = int(toks[1])
        line = prmin.readline()
    natmpwat = int(natm / nwat)

    if natmpwat != 3:
        print("ERROR: number of atoms per water not equals to 3.\n")
        sys.exit()

    prmout.write("%FLAG IPOL\n")
    prmout.write("%FORMAT(1I8)\n")
    prmout.write("%8d\n" % (0))

    prmout.write("%FLAG POL_GAUSS_FORCEFIELD\n")
    prmout.write(
        "%COMMENT This indicates that this parm file is specific to pGM force field\n"
    )
    prmout.write("%COMMENT This must be present if ipgm (in mdin) is 1\n")
    prmout.write("%COMMENT This must NOT be present if ipgm is 0\n")
    prmout.write("%FORMAT(i5)\n")
    prmout.write("%5d\n" % (1))

    prmout.write("%FLAG POL_GAUSS_COVALENT_POINTERS_LIST\n")
    prmout.write("%COMMENT number of covalent dipoles per atom\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * natmpwat, nwat, natmpwat))
    prmout.write("%FORMAT(10I8)\n")
    nptr = {}
    for i in range(natm):
        nptr[i] = 0
    for i in range(nwat):
        for j in range(ndip):
            if p[j][0] <= 3:
                idx = p[j][0] - 1 + i * 3
                # 				print(i, j, idx)
                nptr[idx] += 1

    wrt_int_blk(prmout, nwat, natmpwat, nptr, 0)

    prmout.write("%FLAG POL_GAUSS_COVALENT_ATOMS_LIST\n")
    prmout.write("%COMMENT   dimension = ")
    m = 0
    for j in range(natmpwat):
        m = m + nptr[j]
    prmout.write("(%d=%d*%d)\n" % (nwat * m, nwat, m))
    prmout.write("%FORMAT(10I8)\n")
    lst = {}
    for j in range(m):
        lst[j] = p[j][1]
    wrt_int_blk(prmout, nwat, m, lst, natmpwat)

    prmout.write("%FLAG POL_GAUSS_COVALENT_DIPOLES_LIST\n")
    prmout.write("%COMMENT   unit: e-Angstrom\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * m, nwat, m))
    prmout.write("%FORMAT(5E16.8)\n")
    for j in range(m):
        lst[j] = p[j][2] * au2A
    wrt_flt_blk(prmout, nwat, m, lst)

    prmout.write("%FLAG POL_GAUSS_MONOPOLES_LIST\n")
    prmout.write("%COMMENT   unit: e\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * natmpwat, nwat, natmpwat))
    prmout.write("%FORMAT(5E16.8)\n")
    wrt_flt_blk(prmout, nwat, natmpwat, q)

    prmout.write("%FLAG POL_GAUSS_RADII_LIST\n")
    prmout.write("%COMMENT   unit: Angstrom\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * natmpwat, nwat, natmpwat))
    prmout.write("%FORMAT(5E16.8)\n")

    for j in range(natmpwat):
        rad[j] = rad[j] * au2A
    wrt_flt_blk(prmout, nwat, natmpwat, rad)

    prmout.write("%FLAG POL_GAUSS_POLARIZABILITY_LIST\n")
    prmout.write("%COMMENT   unit: Angstrom**3\n")
    prmout.write("%COMMENT   dimension = ")
    prmout.write("(%d=%d*%d)\n" % (nwat * natmpwat, nwat, natmpwat))
    prmout.write("%FORMAT(5E16.8)\n")
    for j in range(natmpwat):
        pol[j] = pol[j] * au2A**3
    wrt_flt_blk(prmout, nwat, natmpwat, pol)


# -----------------------------------------------------------------------
# The beginning of the main program
# -----------------------------------------------------------------------

###### Get the file names ######
qin, prmtop_in, prmtop_out = file_in()

qin_file = open(qin, "r")
prmtop_in_file = open(prmtop_in, "r")
prmtop_out_file = open(prmtop_out, "w")

###### Read the control parameters ######
natm, ndip, crd, q, p, pol, rad = read_qin()

wrt_out_with_replace(natm, ndip, prmtop_in_file, prmtop_out_file)

# -----------------------------------------------------------------------
# The end of the main program
# -----------------------------------------------------------------------
