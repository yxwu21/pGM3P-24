import os
import re
import math
import numpy as np
import glob

from collections import defaultdict
from timeit import default_timer as timer
from datetime import timedelta


# pdb_prep


def dry_reduce(folder_name, file_name):
    """
    Use the --dry flag to remove crystallographic waters and the --reduce flag to add hydrogen atoms in their optimal locations with pdb4amber.
    """
    print(f"Remove WAT and add H in `{folder_name}`...\n")

    start_time = timer()

    result = os.system(
        f"pdb4amber -i {file_name}.pdb -o {file_name}_dry_H.pdb --dry --reduce"
    )
    if result != 0:
        raise RuntimeError("The simulation is broken.")
    finish_time = timer()
    elapsed_time = timedelta(seconds=finish_time - start_time)

    # report elapsed time
    msg = f"\nProcess finished. Time cost {elapsed_time}."
    print(msg)

    print("-" * 80)

    return {
        "folder": folder_name,
        "file": str(file_name),
        "msg": msg,
        "time": str(elapsed_time),
    }


# lig_del


def read_pdb(file_path):
    print(f"Reading pdb from {file_path}...")
    """
    read pdb file from given path
    """
    file = open(file_path, "r")
    lines = file.readlines()
    return lines


def find_lig(lines, file_name, lig_name):
    print(f"Removing {lig_name} info from {file_name}...")
    """
    remove ligand info from given content and write in new pdb
    """
    with open(f"{file_name}_del_{lig_name}.pdb", "w") as f:
        for line in lines:
            if f"{lig_name}" not in line:
                f.write(line)
    f.close()

    print(f"{lig_name} info has been removed and a new protein pdb has been created")


# cat_protein_lig


def cat_pdbs(folder, protein_pdb, lig_pdb, cat_pdb):
    """
    Combine the protein and ligand PDB files using a cat command
    """
    print(f"Combine {protein_pdb} and {lig_pdb} in {folder}...\n")

    start_time = timer()

    result = os.system(f"cat {protein_pdb}.pdb {lig_pdb}.pdb > {cat_pdb}.pdb")
    if result != 0:
        raise RuntimeError("The simulation is broken.")
    finish_time = timer()
    elapsed_time = timedelta(seconds=finish_time - start_time)

    # report elapsed time
    msg = f"\nProcess finished. Time cost {elapsed_time}."
    print(msg)

    print("-" * 80)

    return {
        "folder": folder,
        "file": str(cat_pdb),
        "msg": msg,
        "time": str(elapsed_time),
    }


# pdb_prep_routine


def prep_pipeline(file_names):
    cwrk = os.path.abspath(os.getcwd())

    """
    Remove WAT and add hydrogen atoms
    """
    for file_name in file_names:
        folder_one = f"/home/yxwu/SHP2/3rd_MD/{file_name}/tleap"
        os.chdir(folder_one)
        dry_reduce(folder_one, file_name)
        os.chdir(cwrk)

        print("*" * 80)

        """
        read pdb file from given path and remove ligand info for a new pdb
        """

        file = f"/home/yxwu/SHP2/3rd_MD/{file_name}/tleap/{file_name}_dry_H.pdb"
        data = read_pdb(file)
        os.chdir(f"/home/yxwu/SHP2/3rd_MD/{file_name}/pdb_prep")
        lig_name = "5OD"
        find_lig(data, file_name, lig_name)
        os.chdir(cwrk)

        print("*" * 80)

        """
        combine the protein and ligand PDB files
        """
        if "unbound" not in file_name:
            folder_two = f"/home/yxwu/SHP2/3rd_MD/{file_name}/pdb_prep"

            os.chdir(folder_two)
            os.system(f"cp /home/yxwu/SHP2/3rd_MD/{file_name}/tleap/Shp099_new.pdb ./")

            protein = f"{file_name}_del_5OD"
            lig = "Shp099_new"
            cat_pdb = f"{file_name}_del_5OD_cat"

            cat_pdbs(folder_two, protein, lig, cat_pdb)
            os.chdir(cwrk)

            print("*" * 80)


# add_ions_calculation


def read_log(file_path):
    """
    read leap.log data from given path
    """
    file = open(file_path, "r")
    lines = file.readlines()
    return lines


def find_volume(lines):
    """
    find volume data from given content
    """
    volumes = []
    for line in lines:
        if "Volume:" in line:
            data = line.split(" ")
            value = float(data[3])
            volumes.append(value)
            print(line)
    return volumes


def find_ions(lines):
    """
    find ions required to neutralize the system
    """
    sodium_ions_required = []
    for line in lines:
        if (
            "Na+ ion required to neutralize." in line
            or "Na+ ions required to neutralize." in line
            or "NA ions required to neutralize." in line
        ):
            data = line.split(" ")
            value = int(data[0])
            sodium_ions_required.append(value)
            print(line)
    return sodium_ions_required


def compute_add_ions(volume, desired_concentration):
    """
    Convert the volume of the systems in  A^3 to liters
    """
    conv_vol = (
        volume
        * math.pow(math.pow(10, 2), 3)
        / (math.pow(math.pow(10, 10), 3) * math.pow(10, 3))
    )
    print("Volume = {} L".format(conv_vol))

    """
    Determine how many chloride ions are present in one liter of solution at a certain concentration
    """
    # desired_conc = float(input("Desired Concentration in mM:"))
    chloride_ions_present = (
        desired_concentration * 6.022 * math.pow(10, 23) / math.pow(10, 3)
    )

    """
    Determine how many chloride ions are needed in the system
    """
    chloride_ions_needed = round(conv_vol * chloride_ions_present)

    """
    Determine the number of sodium ions needed
    """
    sodium_ions_needed = chloride_ions_needed

    print("# Na+ = # Cl- = {}".format(sodium_ions_needed))
    return chloride_ions_needed, sodium_ions_needed


def compute_overall_add_ions(
    sodium_ions_required, chloride_ions_needed, sodium_ions_needed
):
    """
    Compute overall Na+ and Cl- ions needed to add
    """
    sodium_ions_add = sodium_ions_needed + sodium_ions_required
    chloride_ions_add = chloride_ions_needed

    addIons_command = "addIons2 mol Na+ {} Cl- {} ".format(
        sodium_ions_add, chloride_ions_add
    )

    print(addIons_command)
    return addIons_command


def addIons_pipeline(log_file, desired_concentration):
    print("*" * 80)
    print("current computiing: {}".format(log_file))
    data = read_log(log_file)
    volumes = find_volume(data)
    last_volume = volumes[-1]
    sodium_ions = find_ions(data)
    last_sodium_ions = sodium_ions[-1]
    add_ions = compute_add_ions(last_volume, desired_concentration)
    addIons_command = compute_overall_add_ions(
        last_sodium_ions, add_ions[0], add_ions[1]
    )

    print()
    return addIons_command


# balanced_ion_tleap_script


def balanced_script(file_names, desired_concentration):
    print(f"Writing balanced tleap script with add_ions calculation results...")

    for file_name in file_names:
        """
        write balanced tleap input file in given path
        """
        log_path = f"/home/yxwu/SHP2/3rd_MD/{file_name}/pdb_prep/leap.log"

        addIons_command = addIons_pipeline(log_path, desired_concentration)

        print("-" * 80)

        folder_path = f"/home/yxwu/SHP2/3rd_MD/{file_name}/pdb_prep"

        lig_file = "Shp099_gaff2"

        cwrk = os.path.abspath(os.getcwd())

        if "unbound" not in file_name:
            script = [
                "source leaprc.water.tip3p",
                "\n",
                "source leaprc.protein.ff19SB",
                "\n",
                "source leaprc.gaff2",
                "\n",
                f"loadamberparams {lig_file}.frcmod",
                "\n",
                f"5OD=loadmol2 {lig_file}.mol2",
                "\n",
                f"mol =loadpdb {file_name}_del_5OD_cat.pdb",
                "\n",
                f"{addIons_command}",
                "\n",
                "solvateOct mol TIP3PBOX 8.0",
                "\n",
                f"saveamberparm mol {file_name}_del_50D_gaff2.prmtop {file_name}_del_50D_gaff2.inpcrd",
                "\n",
                f"savepdb mol {file_name}_del_50D_gaff2_H.pdb",
                "\n",
                "quit",
            ]

        else:
            script = [
                "source leaprc.water.tip3p",
                "\n",
                "source leaprc.protein.ff19SB",
                "\n",
                "source leaprc.gaff2",
                "\n",
                f"loadamberparams {lig_file}.frcmod",
                "\n",
                f"5OD=loadmol2 {lig_file}.mol2",
                "\n",
                f"mol =loadpdb {file_name}_del_5OD.pdb",
                "\n",
                f"{addIons_command}",
                "\n" "solvateOct mol TIP3PBOX 8.0",
                "\n",
                f"saveamberparm mol {file_name}_del_50D_gaff2.prmtop {file_name}_del_50D_gaff2.inpcrd",
                "\n",
                f"savepdb mol {file_name}_del_50D_gaff2_H.pdb",
                "\n",
                "quit",
            ]

        script_name = f"balanced_{file_name}_tleap"

        os.chdir(folder_path)
        write_script(folder_path, script_name, script)
        os.chdir(cwrk)

        print("*" * 80)


# md_prep


def md_input_prep(file_name, file, location, destination):
    """
    Copy MD input files to destinated folder
    """
    location = f"/home/yxwu/SHP2/1st_MD/{file_name}/md/{file}"
    destination = f"/home/yxwu/SHP2/3rd_MD/{file_name}/MD"

    print(f"Copying {location} to {destination}...")
    os.system(f"cp {location} {destination}")
    print("Copying finished.")


# run_tleap


def start_tleap(folder_name, script):
    print(f"Run tleap in `{folder_name}`...\n")

    start_time = timer()

    result = os.system(f"tleap -s -f {script}")
    if result != 0:
        raise RuntimeError("The simulation is broken.")
    finish_time = timer()
    elapsed_time = timedelta(seconds=finish_time - start_time)

    # report elapsed time
    msg = f"\nRun finished. Time cost {elapsed_time}."
    print(msg)
    return {
        "folder": folder_name,
        "script": str(script),
        "msg": msg,
        "time": str(elapsed_time),
    }


# write_tleap_script


def write_script(folder_path, script_name, script):
    print(f"Writing tleap script in {folder_path}...")
    """
    write tleap input file in given path
    """
    with open(f"{script_name}.in", "w") as f:
        f.writelines(script)
    f.close()

    print(f"{script_name}.in has been created and saved")


# tleap_prep_routine


def tleap_pipeline(file_names):
    cwrk = os.path.abspath(os.getcwd())

    """
    prepare tleap input file
    """

    for file_name in file_names:
        folder_path = f"/home/yxwu/SHP2/3rd_MD/{file_name}/pdb_prep"

        lig_file = "Shp099_gaff2"

        if "unbound" not in file_name:
            script = [
                "source leaprc.water.tip3p",
                "\n",
                "source leaprc.protein.ff19SB",
                "\n",
                "source leaprc.gaff2",
                "\n",
                f"loadamberparams {lig_file}.frcmod",
                "\n",
                f"5OD=loadmol2 {lig_file}.mol2",
                "\n",
                f"mol =loadpdb {file_name}_del_5OD_cat.pdb",
                "\n",
                "addIons2 mol Na+ 0",
                "\n",
                "addIons2 mol Cl- 0",
                "\n",
                "solvateOct mol TIP3PBOX 8.0",
                "\n",
                f"saveamberparm mol {file_name}_del_50D_gaff2.prmtop {file_name}_del_50D_gaff2.inpcrd",
                "\n",
                f"savepdb mol {file_name}_del_50D_gaff2_H.pdb",
                "\n",
                "quit",
            ]

        else:
            script = [
                "source leaprc.water.tip3p",
                "\n",
                "source leaprc.protein.ff19SB",
                "\n",
                "source leaprc.gaff2",
                "\n",
                f"loadamberparams {lig_file}.frcmod",
                "\n",
                f"5OD=loadmol2 {lig_file}.mol2",
                "\n",
                f"mol =loadpdb {file_name}_del_5OD.pdb",
                "\n",
                "addIons2 mol Na+ 0",
                "\n",
                "addIons2 mol Cl- 0",
                "\n",
                "solvateOct mol TIP3PBOX 8.0",
                "\n",
                f"saveamberparm mol {file_name}_del_50D_gaff2.prmtop {file_name}_del_50D_gaff2.inpcrd",
                "\n",
                f"savepdb mol {file_name}_del_50D_gaff2_H.pdb",
                "\n",
                "quit",
            ]

        script_name = f"{file_name}_tleap"

        os.chdir(folder_path)
        os.system(f"cp /home/yxwu/SHP2/3rd_MD/{file_name}/tleap/{lig_file}.frcmod ./")
        os.system(f"cp /home/yxwu/SHP2/3rd_MD/{file_name}/tleap/{lig_file}.mol2 ./")
        write_script(folder_path, script_name, script)
        os.chdir(cwrk)

        print("*" * 80)

        """
        run tleap
        """

        folder_path = f"/home/yxwu/SHP2/3rd_MD/{file_name}/pdb_prep"

        script_run = f"{file_name}_tleap.in"

        os.chdir(folder_path)
        start_tleap(folder_path, script_run)
        os.chdir(cwrk)

        print("*" * 160)


# write_MD_script


def read_file(ref_file_path):
    print(f"Reading file from {ref_file_path}...")
    """
    Read in file from given path
    """

    file = open(ref_file_path, "r")
    lines = file.readlines()

    return lines


def modi_script(lines, folder_path, script_name, content, new_content):
    print(f"Writing new script in {folder_path}...")
    """
    Modify MD script with replacing the target string
    """
    with open(f"{script_name}", "w") as f:
        for line in lines:
            new_script = line.replace(f"{content}", f"{new_content}")
            f.write(new_script)
    f.close()

    print(f"New script in {folder_path} has been created and saved.")


def find_intValue(lines, identifier, index):
    """
    find specific data for identifier from given content
    """
    for line in lines:
        if f"{identifier}" in line:
            line_strip = line.strip()
            datas = line_strip.split(" ")
            # print(datas)

            valid_data = []
            for data in datas:
                if data:
                    valid_data.append(data)
            # print(valid_data)

            index_num = int(f"{index}")
            value = int(valid_data[index_num])
            # print(line)
    return value


def find_floatValue(lines, identifier, index):
    """
    find specific data for identifier from given content
    """
    for line in lines:
        if f"{identifier}" in line:
            line_strip = line.strip()
            datas = line_strip.split(" ")
            # print(datas)

            valid_data = []
            for data in datas:
                if data:
                    valid_data.append(data)
            # print(valid_data)

            index_num = int(f"{index}")
            value = float(valid_data[index_num])
            # print(line)
    return value


def read_dat(file):
    print(f"Reading data from {file}...")
    """
    Read data from given path
    """
    data = []
    with open(f"{file}", "r") as f:
        reader = f.readlines()
        for line in reader:
            row = line.strip()
            get_col = row.split()
            data.append([float(i) for i in get_col])
    return data


def accuracy(predict, target):
    equal_mask = predict == target
    acc = np.mean(equal_mask.astype(np.float32))
    return acc


def get_protein(path, delimiter):
    os.chdir(path)
    dirs = glob.glob(delimiter)

    proteins = []
    for dir in dirs:
        protein = dir.split(".")
        proteins.append(protein[0])

    return proteins


# construct_dssp_file


def write_dssp(num, type, file_name, folder_path):
    os.chdir(folder_path)

    with open(f"{file_name}.dssp", "w") as f:
        f.write(num)
        f.write("\n")
        f.write(f"{type}" * int(num))

    print(f"{file_name}.dssp has been created in {folder_path}")


# construct_inp_file


def write_inp(
    cmpd,
    cmpd_comb,
    cmpd_num,
    atp_num,
    cmpd_inside_box,
    atp_inside_box,
    box_sides,
    file_name,
    folder_path,
):
    print(f"Writing {file_name}_mixture.inp script...")

    os.chdir(folder_path)

    if "ATP" not in file_name:
        script = [
            "#",
            "\n",
            "# A mixture of water and urea",
            "\n",
            "#",
            "\n",
            "tolerance 2.0",
            "\n",
            "filetype pdb",
            "\n",
            f"output {cmpd_comb}_{cmpd_num}_vac.pdb",
            "\n",
            "\n",
            f"structure {cmpd}-CG.pdb",
            "\n",
            f"  number {cmpd_num}",
            "\n",
            f"  inside box {cmpd_inside_box}",
            "\n",
            "end structure",
            "\n",
            "\n",
            f"add_box_sides {box_sides}",
        ]
    else:
        script = [
            "#",
            "\n",
            "# A mixture of water and urea",
            "\n",
            "#",
            "\n",
            "tolerance 2.0",
            "\n",
            "filetype pdb",
            "\n",
            f"output {cmpd_comb}_{cmpd_num}_vac.pdb",
            "\n",
            "\n",
            f"structure {cmpd}-CG.pdb",
            "\n",
            f"  number {cmpd_num}",
            "\n",
            f"  inside box {cmpd_inside_box}",
            "\n",
            "end structure",
            "\n",
            "\n",
            "structure ATP.pdb",
            "\n",
            f"  number {atp_num}",
            "\n",
            f"  inside box {atp_inside_box}",
            "\n",
            "end structure",
            "\n",
            "\n",
            f"add_box_sides {box_sides}",
        ]

    with open(f"{file_name}_mixture.inp", "w") as f:
        f.writelines(script)

    print(f"{file_name}_mixture.inp has been created in {folder_path}")


# water molecule number


def findW(lines):
    """
    find water molecule number
    """
    W_num = []
    for line in lines:
        if "W (   1 atoms):" in line:
            data = line.split(": ")[-1]
            data = data.split(" ")[0]
            value = int(data)
            W_num.append(value)
            print(line)
    return W_num


# Construct top file


def write_top(cmpd, cmpd_num, atp_num, W_num, file_name, folder_path):
    print(f"Writing mixture.inp script...")

    os.chdir(folder_path)

    if "ATP" not in file_name:
        script = [
            '#include "martini_v2.2.itp"',
            "\n",
            "\n",
            "\n",
            '#include "Protein_X.itp""',
            "\n",
            "\n",
            "[ system ]",
            "\n",
            "; name",
            "\n",
            f"Martini system from {cmpd}_relax.pdb",
            "\n",
            "\n",
            "[ molecules ]",
            "\n",
            "; name        number",
            "\n",
            f"Protein_X 	 {cmpd_num}",
            "\n",
            f"W			{W_num} ",
        ]
    else:
        script = [
            '#include "martini_v2.1-dna.itp"',
            "\n",
            "\n",
            "\n",
            '#include "Protein_X.itp"',
            "\n",
            '#include "ATP.itp"',
            "\n",
            "\n",
            "[ system ]",
            "\n",
            "; name",
            "\n",
            f"Martini system from {cmpd}_ATP_relax.pdb",
            "\n",
            "\n",
            "[ molecules ]",
            "\n",
            "; name        number",
            "\n",
            f"Protein_X 	 {cmpd_num}",
            "\n",
            f"ATP          {atp_num}",
            "\n",
            f"W			{W_num} ",
        ]

    with open(f"{file_name}_sol.top", "w") as f:
        f.writelines(script)

    print(f"{file_name}_sol.top has been created in {folder_path}")


# Construct neu_top file


def write_neu_top(file_name, folder_path):
    print(f"Reading {file_name}_sol.top file...")
    os.chdir(folder_path)
    file = open(f"{folder_path}/{file_name}_sol.top", "r")
    lines = file.readlines()

    print("write a new top file for neutralization...")

    with open(f"{file_name}_neu.top", "w") as f:
        for line in lines:
            if "Protein_X.itp" in line:
                f.write(line)
                f.write('#include "martini_v2.0_ions.itp"')
                f.write("\n")
                continue
            f.write(line)

    print(f"{file_name}_neu.top has been created")


# Construct single_top file


def write_sin_top(cmpd, pro_num, folder_path):
    print(f"Writing mixture.inp script...")
    os.chdir(folder_path)

    script = [
        '#include "martini.itp"',
        "\n",
        "\n",
        "\n",
        '#include "Protein_X.itp"',
        "\n",
        "\n",
        "[ system ]",
        "\n",
        "; name",
        "\n",
        f"Martini system from {cmpd}_relax.pdb",
        "\n",
        "\n",
        "[ molecules ]",
        "\n",
        "; name        number",
        "\n",
        f"Protein_X 	 {pro_num}",
        "\n",
    ]

    with open(f"single_{cmpd}.top", "w") as f:
        f.writelines(script)

    print(f"single_{cmpd}.top has been created in {folder_path}")


# 01_modelConstruct pipeline


def modelConstruct01(
    conc,
    num_lys,
    proj_path,
    type,
    pro_num,
    cmpd_num,
    atp_num,
    cmpd_inside_box,
    atp_inside_box,
    box_sides,
):
    """
    01_modelConstruct pipeline
    """
    cmpd = f"K{num_lys}"
    cmpd_comb = f"{cmpd}_ATP"
    cmpd_comb_conc = f"{cmpd_comb}_{conc}M"
    # scripts = f"{proj_path}/scripts"
    ref = f"{proj_path}/references"
    main = f"{proj_path}/main"
    folder = f"{proj_path}/{cmpd_comb_conc}"

    file_name = f"{cmpd_comb}_{cmpd_num}"
    mixture_inp = f"{cmpd_comb}_{cmpd_num}_mixture"

    modelPath = f"{folder}/model_construct"
    # mdPath = f"{folder}/MD"
    minPath = f"{folder}/MD/1_minimization"
    equiPath = f"{folder}/MD/2_equilibration"
    prodPath = f"{folder}/MD/3_product"
    prepPath = f"{folder}/MD/prep_files"

    # cwrk = os.path.abspath(os.getcwd())

    os.makedirs(f"{modelPath}", exist_ok=True)
    # os.system(f"cp {scripts}/{cmpd_comb}_relax.pdb {modelPath}")

    os.makedirs(f"{minPath}", exist_ok=True)
    os.system(f"cp {ref}/minimization.mdp {minPath}")

    os.makedirs(f"{equiPath}", exist_ok=True)
    os.system(f"cp {ref}/equilibration.mdp {equiPath}")

    os.makedirs(f"{prodPath}", exist_ok=True)
    os.system(f"cp {ref}/product.mdp {prodPath}")

    os.makedirs(f"{prepPath}", exist_ok=True)

    os.chdir(modelPath)

    # Model Construction 01

    """
    starting files
    """
    os.system(f"cp {main}/{cmpd}_relax.pdb ./")
    os.system(f"cp {ref}/ATP.pdb ./")
    write_sin_top(cmpd, pro_num, modelPath)

    """
    Martini_prep : prep dssp file
    """
    write_dssp(num_lys, type, cmpd_comb, modelPath)

    """
    Perform Martini
    """
    os.system(
        f"martinize.py -f {cmpd}_relax.pdb -o single_{cmpd}.top -x {cmpd}-CG.pdb -ss {cmpd_comb}.dssp -p backbone -ff martini22"
    )

    """
    packmol
    """
    write_inp(
        cmpd,
        cmpd_comb,
        cmpd_num,
        atp_num,
        cmpd_inside_box,
        atp_inside_box,
        box_sides,
        file_name,
        modelPath,
    )
    os.system(f"packmol < {mixture_inp}.inp")

    """
    Add solvent
    """
    os.chdir(ref)
    os.system(f"cp vac.mdp ions.mdp water.gro {modelPath}")
    os.chdir(modelPath)
    os.system(
        f"gmx solvate -cp {cmpd_comb}_{cmpd_num}_vac.pdb -cs water.gro -radius 0.21 -o {cmpd_comb}_{cmpd_num}_sol.gro &> {cmpd_comb}_{cmpd_num}_sol.gro.log"
    )

    """
    Read gro.log and get water molecule number
    """
    gro_log_file = f"{modelPath}/{cmpd_comb}_{cmpd_num}_sol.gro.log"
    data = read_log(gro_log_file)
    W_num = findW(data)[0]
    write_top(cmpd, cmpd_num, atp_num, W_num, file_name, modelPath)

    """
    prep add ions
    """
    os.chdir(ref)
    os.system(f"cp *.itp {modelPath}")
    os.chdir(modelPath)
    os.system(
        f"gmx grompp -f ions.mdp -c {cmpd_comb}_{cmpd_num}_sol.gro -p {cmpd_comb}_{cmpd_num}_sol.top -o {cmpd_comb}_{cmpd_num}_sol.tpr &> {cmpd_comb}_{cmpd_num}_sol.tpr.log"
    )

    """
    add ions: neutralization
    """
    write_neu_top(file_name, modelPath)
    os.system(
        f"gmx genion -s {cmpd_comb}_{cmpd_num}_sol.tpr -o {cmpd_comb}_{cmpd_num}_neu.gro -p {cmpd_comb}_{cmpd_num}_neu.top -pname NA -nname CL -neutral"
    )

    print("01_modelConstruction completed")


# 02_modelConstruct pipeline


def modelConstruct02(conc, num_lys, proj_path, cmpd_num):
    """
    02_modelConstruct pipeline
    """
    cmpd = f"K{num_lys}"
    cmpd_comb = f"{cmpd}_ATP"
    cmpd_comb_conc = f"{cmpd_comb}_{conc}M"
    # scripts = f"{proj_path}/scripts"
    ref = f"{proj_path}/references"
    # main = f"{proj_path}/main"
    folder = f"{proj_path}/{cmpd_comb_conc}"

    # file_name = f"{cmpd_comb}_{cmpd_num}"
    # mixture_inp = f"{cmpd_comb}_{cmpd_num}_mixture"

    modelPath = f"{folder}/model_construct"
    # mdPath = f"{folder}/MD"
    minPath = f"{folder}/MD/1_minimization"
    equiPath = f"{folder}/MD/2_equilibration"
    prodPath = f"{folder}/MD/3_product"
    prepPath = f"{folder}/MD/prep_files"

    # cwrk = os.path.abspath(os.getcwd())

    os.makedirs(f"{modelPath}", exist_ok=True)
    # os.system(f"cp {scripts}/{cmpd_comb}_relax.pdb {modelPath}")

    os.makedirs(f"{minPath}", exist_ok=True)
    os.system(f"cp {ref}/minimization.mdp {minPath}")

    os.makedirs(f"{equiPath}", exist_ok=True)
    os.system(f"cp {ref}/equilibration.mdp {equiPath}")

    os.makedirs(f"{prodPath}", exist_ok=True)
    os.system(f"cp {ref}/product.mdp {prodPath}")

    os.makedirs(f"{prepPath}", exist_ok=True)

    os.chdir(modelPath)

    # Model Construction 02

    """
    add ions: neutralization - generate tpr file
    """
    os.system(
        f"gmx grompp -f ions.mdp -c {cmpd_comb}_{cmpd_num}_neu.gro -p {cmpd_comb}_{cmpd_num}_neu.top -o {cmpd_comb}_{cmpd_num}_neu.tpr"
    )

    """
    add ions: specified concentration
    """
    os.system(
        f"cp {cmpd_comb}_{cmpd_num}_neu.top ./{cmpd_comb_conc}_{cmpd_num}_conc.top"
    )
    os.system(
        f"gmx genion -s {cmpd_comb}_{cmpd_num}_neu.tpr -o {cmpd_comb_conc}_{cmpd_num}_conc.gro -p {cmpd_comb_conc}_{cmpd_num}_conc.top -pname NA -nname CL -conc {conc}"
    )

    """
    Move files to MD/prepfile
    """
    os.system(f"cp {ref}/*.itp {prepPath}")
    os.system(f"cp {modelPath}/{cmpd_comb_conc}_{cmpd_num}_conc.* {prepPath}")

    print("02_modelConstruction completed")
    print(f"Model construction for {cmpd_comb_conc} completed")


# ALDH


def getNA_num(log_file, cl_num):
    print("*" * 80)
    print("current computiing: {}".format(log_file))
    data = read_log(log_file)
    sodium_ions = find_ions(data)
    last_sodium_ions = sodium_ions[-1]
    na_num = int(cl_num) + last_sodium_ions

    addions_na_command = "addions complex NA {} ".format(na_num)
    addions_cl_command = "addions complex CL {}".format(cl_num)

    print(addions_na_command)
    print(addions_cl_command)

    return addions_na_command, addions_cl_command


def balanced_tleap(
    folder_path, tleap_ref, mutation, addions_na_command, addions_cl_command
):
    print(f"Reading tleap.in template...")
    os.chdir(folder_path)
    file = open(f"{tleap_ref}", "r")
    lines = file.readlines()

    print("write a balanced tleap.in file for mutant...")

    with open(f"tleap_{mutation}.in", "w") as f:
        for line in lines:
            if "addions complex NA" in line:
                f.write(addions_na_command)
                f.write("\n")
                continue
            elif "addions complex CL" in line:
                f.write(addions_cl_command)
                f.write("\n")
                continue
            f.write(line)

    print(f"tleap_{mutation}.in has been created")

    balanced_tleap_name = f"tleap_{mutation}.in"

    return balanced_tleap_name


def ALDH_mt_tleapConstruct(state, Type, mutation, cl_num, folder_path):
    ref_path = "/home/yxwu/ALDH_HanLi/references"

    # folder_path = f"/home/yxwu/ALDH_HanLi/MT/{mutation}/{state}/{Type}/model_construct"
    os.makedirs(f"{folder_path}", exist_ok=True)
    os.chdir(f"{folder_path}")

    tleap_ref = f"{ref_path}/tleap_{Type}.in"
    tleap_lines = read_file(tleap_ref)
    tleap_script = f"tleap_{mutation}_v1.in"
    content = "{Mutation}"
    new_content = mutation
    modified_tleap = modi_script(
        tleap_lines, folder_path, tleap_script, content, new_content
    )

    tleap_ref_2 = f"{folder_path}/{tleap_script}"
    tleap_lines_2 = read_file(tleap_ref_2)
    tleap_script_2 = f"tleap_{mutation}_v2.in"
    content_2 = "{state}"
    new_content_2 = state
    modified_tleap_2 = modi_script(
        tleap_lines_2, folder_path, tleap_script_2, content_2, new_content_2
    )

    os.system(f"tleap -s -f {tleap_script_2}")

    log_file = f"{folder_path}/leap.log"
    addions_na_command = getNA_num(log_file, cl_num)[0]
    addions_cl_command = getNA_num(log_file, cl_num)[1]

    modi_tleap_path = f"{folder_path}/{tleap_script_2}"
    tleap_v2 = balanced_tleap(
        folder_path, modi_tleap_path, mutation, addions_na_command, addions_cl_command
    )

    os.system(f"rm {tleap_script}")
    os.system(f"rm {tleap_script_2}")
    os.system("rm extend*")

    return tleap_v2


def sglMt_pdbConstruct(
    resIndex, ref_pdb_path, original_res, mutate_res, res_length, pdb_name, presAtomType
):
    lines = read_file(ref_pdb_path)

    pat = re.compile(f"ATOM\s+[0-9]+\s+([A-Z0-9]+)\s+([A-Z]+)\sX\s+{resIndex}")

    new_pdb_lines = []

    with open(f"{pdb_name}.pdb", "w") as f:
        for line in lines:
            line: str = line.strip()
            if pat.match(line):
                type = pat.findall(line)[0][0]
                if type in presAtomType:
                    list_of_line = list(line)
                    original_res_left = line.index(original_res)
                    list_of_line[
                        original_res_left : original_res_left + res_length
                    ] = mutate_res
                    line = "".join(list_of_line)
                else:
                    continue
            new_pdb_lines.append(line)
            new_pdb = "\n".join(new_pdb_lines)
        f.write(new_pdb)

    new_pdb = f"{pdb_name}.pdb"
    return new_pdb, original_res_left


def check_modelConstruction(resIndex, new_pdb_path, mutate_res, mutation):
    print("Checking if the model is constructed correctly...")

    resid = str(int(resIndex) + 1)
    # new_pdb_path = f"{folder_path}/{new_pdb}"
    new_pdb_lines = read_file(new_pdb_path)

    new_pat = re.compile(f"ATOM\s+[0-9]+\s+[A-Z0-9]+\s+([A-Z]+)\s+{resid}\s+")

    modelCheckResult = False
    for line in new_pdb_lines:
        line: str = line.strip()
        if new_pat.match(line):
            # list_of_line = list(line)
            res = new_pat.findall(line)
            if len(res) != 1 or res[0] != mutate_res:
                print("Error in validating model construction.")
                modelCheckResult = True
            else:
                print(f"Model for {mutation} system is correctly constructed.")
                modelCheckResult = False

    return modelCheckResult


def find_catWat(sol_pdb_path, sol_pdb_name):
    print("Finding the resid of catalytic water...")

    sol_pdb_lines = read_file(sol_pdb_path)
    sol_pat = re.compile("ATOM\s+[0-9]+\s+[A-Z0-9]+\s+CL+\s+([0-9]+)\s+")

    cl_resid_tot = []
    for line in sol_pdb_lines:
        line: str = line.strip()
        if sol_pat.match(line):
            cl_resid = sol_pat.findall(line)
            cl_resid_tot.append(cl_resid)
    lastCl_resid = cl_resid[-1]
    catWat_resid = int(lastCl_resid) + 1

    print(f"Catalytic water resid in {sol_pdb_name}:", catWat_resid)
    return catWat_resid


def convert_dot2xyz(dot_file, xyz_file):
    lines = None
    with open(dot_file) as f:
        lines = f.readlines()

    atom_num = None
    for line in reversed(lines):
        if norm_line := line.strip():
            terms = norm_line.split()

            if terms[0].startswith("H"):
                atom_num = terms[0][1:]
                break

    lines[0] = f"{atom_num}\n"
    with open(xyz_file, "w") as f:
        f.write("".join(lines))


def geo_sglMt_pdbConstruct(
    geo_pat, ref_pdb_path, original_res, mutate_res, res_length, pdb_name, presAtomType
):
    lines = read_file(ref_pdb_path)

    pat = re.compile(geo_pat)

    new_pdb_lines = []

    with open(f"{pdb_name}.pdb", "w") as f:
        for line in lines:
            line: str = line.strip()
            if pat.match(line):
                type = pat.findall(line)[0][0]
                if type in presAtomType:
                    list_of_line = list(line)
                    original_res_left = line.index(original_res)
                    list_of_line[
                        original_res_left : original_res_left + res_length
                    ] = list(mutate_res)
                    line = "".join(list_of_line)
                else:
                    continue
            new_pdb_lines.append(line)
            new_pdb = "\n".join(new_pdb_lines)
        f.write(new_pdb)

    new_pdb = f"{pdb_name}.pdb"
    return new_pdb, original_res_left


def geo_ALDH_mt_tleapConstruct(
    ref_path, geo_type, state, Type, mutation, cl_num, folder_path
):
    os.makedirs(f"{folder_path}", exist_ok=True)
    os.chdir(f"{folder_path}")

    tleap_ref = f"{ref_path}/tleap_{Type}_geo.in"
    tleap_lines = read_file(tleap_ref)
    tleap_script = f"tleap_{mutation}_v1.in"
    content = "{Mutation}"
    new_content = mutation
    modified_tleap = modi_script(
        tleap_lines, folder_path, tleap_script, content, new_content
    )

    tleap_ref_2 = f"{folder_path}/{tleap_script}"
    tleap_lines_2 = read_file(tleap_ref_2)
    tleap_script_2 = f"tleap_{mutation}_v2.in"
    content_2 = "{state}"
    new_content_2 = state
    modified_tleap_2 = modi_script(
        tleap_lines_2, folder_path, tleap_script_2, content_2, new_content_2
    )

    tleap_ref_3 = f"{folder_path}/{tleap_script_2}"
    tleap_lines_3 = read_file(tleap_ref_3)
    tleap_script_3 = f"tleap_{mutation}_v3.in"
    content_3 = "{geo_type}"
    new_content_3 = geo_type
    modified_tleap_3 = modi_script(
        tleap_lines_3, folder_path, tleap_script_3, content_3, new_content_3
    )

    os.system(f"tleap -s -f {tleap_script_3}")

    log_file = f"{folder_path}/leap.log"
    addions_na_command = getNA_num(log_file, cl_num)[0]
    addions_cl_command = getNA_num(log_file, cl_num)[1]

    modi_tleap_path = f"{folder_path}/{tleap_script_3}"
    tleap_v3 = balanced_tleap(
        folder_path, modi_tleap_path, mutation, addions_na_command, addions_cl_command
    )

    os.system(f"rm {tleap_script}")
    os.system(f"rm {tleap_script_2}")
    os.system(f"rm {tleap_script_3}")

    return tleap_v3


# Water Model


def edit_bond_length(tar_prmtop_path: str, ref_prmtop_path: str, args: list):
    lines = read_file(ref_prmtop_path)

    new_lines = []
    hit = False
    edit_log = {}
    template = "  {:.8E}  {:.8E}\n"
    for line in lines:
        if not hit:
            new_lines.append(line)

            if "%FLAG BOND_EQUIL_VALUE" in line:
                hit = True
        else:
            if "%FORMAT(5E16.8)" in line:
                new_lines.append(line)
                continue

            # edit the number value
            edit_log["old"] = line
            edit_line = template.format(*args)
            new_lines.append(edit_line)
            edit_log["new"] = edit_line
            hit = False

    with open(tar_prmtop_path, "w") as f:
        f.writelines(new_lines)
    return edit_log


def edit_vdw_parms(tar_prmtop_path: str, ref_prmtop_path: str, args: list, factor_type):
    lines = read_file(ref_prmtop_path)

    if factor_type == "a":
        pattern = "%FLAG LENNARD_JONES_ACOEF"
    else:
        pattern = "%FLAG LENNARD_JONES_BCOEF"

    new_lines = []
    hit = False
    edit_log = {}
    template = "  {:.8E}  {:.8E}  {:.8E}\n"
    for line in lines:
        if not hit:
            new_lines.append(line)

            if f"{pattern}" in line:
                hit = True
        else:
            if "%FORMAT(5E16.8)" in line:
                new_lines.append(line)
                continue

            # edit the number value
            edit_log["old"] = line
            edit_line = template.format(*args)
            new_lines.append(edit_line)
            edit_log["new"] = edit_line
            hit = False

    with open(tar_prmtop_path, "w") as f:
        f.writelines(new_lines)
    return edit_log


def get_resIndex(line_pat, ref_pdb_path, cap_type):
    lines = read_file(ref_pdb_path)

    pat = re.compile(line_pat)

    capIndex_list = []

    for line in lines:
        line = line.strip()
        if pat.match(line):
            cap = pat.findall(line)[0][0]
            capIndex = pat.findall(line)[0][1]
            if cap == f"{cap_type}":
                capIndex_list.append(capIndex)
    print(f"{cap_type} Index: {capIndex_list}")
    return capIndex_list


def read_md_data(analysis_path, gtStats, data):
    for gtStat in gtStats:
        data_path = f"{analysis_path}/{gtStat}"
        systems = os.listdir(f"{data_path}")
        for system in tqdm(systems):
            system_path = f"{data_path}/{system}"
            states = os.listdir(f"{system_path}")
            for state in tqdm(states):
                state_path = f"{system_path}/{state}"
                Types = os.listdir(f"{state_path}")
                for Type in tqdm(Types):
                    type_path = f"{state_path}/{Type}"
                    analy_folder = f"{type_path}/analysis"
                    metrics = os.listdir(analy_folder)
                    for metric in metrics:
                        metricAnalysisDir = f"{analy_folder}/{metric}"
                        if "rmsf" in metric:
                            rmsf_folders = os.listdir(metricAnalysisDir)
                            for i in rmsf_folders:
                                if ".dat" in i:
                                    rmsf_data = f"{metricAnalysisDir}/{i}"
                                    with open(rmsf_data) as f:
                                        start = False
                                        file_name = os.path.split(rmsf_data)[-1]
                                        for line in f:
                                            if start:
                                                x, y = line.split()
                                                data[file_name].append(
                                                    [float(x), float(y)]
                                                )
                                                continue
                                            if "@type xy" in line:
                                                start = True
                        elif "distance" in metric:
                            dist_folders = os.listdir(metricAnalysisDir)
                            for i in dist_folders:
                                if ".dat" in i:
                                    dist_data = f"{metricAnalysisDir}/{i}"
                                    with open(dist_data) as f:
                                        start = False
                                        file_name = os.path.split(dist_data)[-1]
                                        for line in f:
                                            if start:
                                                x, y = line.split()
                                                data[file_name].append(
                                                    [float(x), float(y)]
                                                )
                                                continue
                                            if "#Frame" in line:
                                                start = True
    return data


def get_next_alphabet(char):
    """
    Given an alphabet, return the next alphabet. If given 'z' or 'Z', wrap around to 'a' or 'A' respectively.
    """
    if char == "z":
        return "a"
    elif char == "Z":
        return "A"
    return chr(ord(char) + 1)


def list_folders(path, exclude_folders):
    try:
        # List all entries in the directory
        all_entries = os.listdir(path)

        # Filter entries, keeping only the folders (excluding the ones in 'exclude_folders')
        folders = [
            entry
            for entry in all_entries
            if os.path.isdir(os.path.join(path, entry)) and entry not in exclude_folders
        ]

        return folders

    except FileNotFoundError:
        print(f"The system cannot find the specified path: {path}")
        return []
    except PermissionError:
        print(f"Permission denied for the specified path: {path}")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []


def get_traj_info(traj_folders, md_folder_path):
    trajectory = []

    parm_pattern = f"{md_folder_path}/prep_files/*.prmtop"
    parm_path = glob.glob(parm_pattern)[0]
    parmtop = parm_path.split("/")[-1]
    parmtop = f"../../MD/prep_files/{parmtop}"

    for traj_folder in traj_folders:
        traj_path_pattern = f"{md_folder_path}/{traj_folder}/*.nc"

        # Get all matching folders
        traj_paths = sorted(glob.glob(traj_path_pattern))
        if traj_paths:
            for traj_path in traj_paths:
                trajectory.append(traj_path)

    def safe_extract_number(filename):
        parts = filename.split(".")
        # Check if the second last part is a number
        if parts[-2].isdigit():
            return int(parts[-2])
        return 0

    trajin_lists = []
    for i in sorted(
        trajectory, key=lambda x: (x.split("/")[-2], safe_extract_number(x))
    ):
        traj_case = i.split("/")[-2]
        traj_nc = i.split("/")[-1]
        relative_traj_nc = f"trajin ../../MD/{traj_case}/{traj_nc}"
        trajin_lists.append(relative_traj_nc)

    """
    back hack to get the last frame of the last 100ps trajectory
    """
    trajin_content = "\n".join(trajin_lists[-11:-1])
    return trajectory, trajin_content, parmtop


def get_output_info(traj_folders, md_folder_path):
    traj_out = []

    parm_pattern = f"{md_folder_path}/prep_files/*.prmtop"
    parm_path = glob.glob(parm_pattern)[0]
    parmtop = parm_path.split("/")[-1]
    parmtop = f"../../MD/prep_files/{parmtop}"

    for traj_folder in traj_folders:
        traj_path_pattern = f"{md_folder_path}/{traj_folder}/*.out"
        # Get all matching folders
        traj_paths = sorted(glob.glob(traj_path_pattern))
        if traj_paths:
            for traj_path in traj_paths:
                traj_out.append(traj_path)

    return traj_out


def read_analysis_data(analysis_folder_path, data_file_pattern, identifier):
    data_path = glob.glob(data_file_pattern)

    data = []
    datafile_name = None
    if data_path:
        data_file = data_path[0]
        datafile_name = data_file.split("/")[-1]
        # Get all matching folders
        with open(f"{analysis_folder_path}/{datafile_name}") as f:
            start = False
            for line in f:
                if start:
                    x, y = line.split()
                    data.append([float(x), float(y)])
                    continue
                if identifier in line:
                    start = True

    return data, datafile_name


# def modify_and_write_qin(input_qin_path, output_qin_path, q_scale, p_scale, pol_scale, rad_scale):
#     with open(input_qin_path, 'r') as file:
#         lines = file.readlines()

#     new_lines = []
#     hit = False
#     edit_log = {}

#     for line in lines:
#         if not hit:
#             new_lines.append(line)

#             if '%FLAG ATOM CRD POL RAD:' in line:
#                 hit = True
#         else:
#             if 'atm.no' in line:
#                 new_lines.append(line)
#                 continue

#             parts = line.split()
#             if len(parts) == 0:
#                 hit = False
#                 new_lines.append(line)
#                 continue
#             else:
#                 # edit the number value
#                 edit_log['old POL'] = parts[4]
#                 edit_log['old RAD'] = parts[5]

#                 parts[4] = f"{float(parts[4]) * pol_scale:.7E}"
#                 parts[5] = f"{float(parts[5]) * rad_scale:.7E}"

#                 edit_log['new POL'] = parts[4]
#                 edit_log['new RAD'] = parts[5]

#                 edit_line = "{:>4} {:>15} {:>15} {:>15} {:>15} {:>15} {:>5} {:>5}\n".format(*parts)
#                 new_lines.append(edit_line)

#     sec_new_lines = []
#     for line in new_lines:
#         if not hit:
#             sec_new_lines.append(line)

#             if '%FLAG ATOM CHRG:' in line:
#                 hit = True
#         else:
#             if 'atm.no' in line:
#                 sec_new_lines.append(line)
#                 continue

#             parts = line.split()
#             if len(parts) == 0:
#                 hit = False
#                 sec_new_lines.append(line)
#                 continue
#             else:
#                 # edit the number value
#                 edit_log['old q(opt)'] = parts[-1]

#                 parts[-1] = f"{float(parts[-1]) * q_scale:.7E}"

#                 edit_log['new q(opt)'] = parts[-1]

#                 edit_line = "{:>4} {:>10} {:>10} {:>17}\n".format(*parts)
#                 sec_new_lines.append(edit_line)

#     thrd_new_lines = []
#     for line in sec_new_lines:
#         if not hit:
#             thrd_new_lines.append(line)

#             if '%FLAG PERM DIP LOCAL:' in line:
#                 hit = True
#         else:
#             if 'dip.no' in line:
#                 thrd_new_lines.append(line)
#                 continue

#             parts = line.split()
#             if len(parts) == 0:
#                 hit = False
#                 thrd_new_lines.append(line)
#                 continue
#             else:
#                 # edit the number value
#                 edit_log['old p(opt)'] = parts[-1]

#                 parts[-1] = f"{float(parts[-1]) * q_scale:.7E}"

#                 edit_log['new p(opt)'] = parts[-1]

#                 edit_line = "{:>4} {:>8} {:>8} {:>8} {:>17}\n".format(*parts)
#                 thrd_new_lines.append(edit_line)

#     with open(output_qin_path, 'w') as file:
#         file.writelines(thrd_new_lines)


# def modify_and_write_qin(input_qin_path, output_qin_path, q_scale, p_scale, pol_scale, rad_scale):
#     with open(input_qin_path, 'r') as file:
#         lines = file.readlines()

#     new_lines = []
#     section = None

#     for line in lines:
#         if '%FLAG ATOM CRD POL RAD:' in line:
#             section = 'POL_RAD'
#         elif '%FLAG ATOM CHRG:' in line:
#             section = 'CHRG'
#         elif '%FLAG PERM DIP LOCAL:' in line:
#             section = 'PERM_DIP'
#         elif 'atm.no' in line or 'dip.no' in line:
#             section = None

#         if section == 'POL_RAD' and len(line.split()) >= 6:
#             parts = line.split()
#             parts[4] = f"{float(parts[4]) * pol_scale:.7E}"
#             parts[5] = f"{float(parts[5]) * rad_scale:.7E}"

#             line = "{:>4} {:>15} {:>15} {:>15} {:>15} {:>15} {:>5} {:>5}\n".format(*parts)
#         elif section == 'CHRG' and len(line.split()) >= 4:
#             parts = line.split()
#             parts[-1] = f"{float(parts[-1]) * q_scale:.7E}"

#             line = "{:>4} {:>10} {:>10} {:>17}\n".format(*parts)
#         elif section == 'PERM_DIP' and len(line.split()) >= 5:
#             parts = line.split()
#             parts[-1] = f"{float(parts[-1]) * p_scale:.7E}"

#             line = "{:>4} {:>8} {:>8} {:>8} {:>17}\n".format(*parts)

#         new_lines.append(line)

#     with open(output_qin_path, 'w') as file:
#         file.writelines(new_lines)


def modify_and_write_qin(
    input_qin_path, output_qin_path, q_scale, p_scale, pol_scale, rad_scale
):
    with open(input_qin_path, "r") as file:
        lines = file.readlines()

    new_lines = []
    section = None

    for line in lines:
        # Check and set the section based on the header
        if "%FLAG ATOM CRD POL RAD:" in line:
            section = "POL_RAD"
        elif "%FLAG ATOM CHRG:" in line:
            section = "CHRG"
        elif "%FLAG PERM DIP LOCAL:" in line:
            section = "PERM_DIP"

        parts = line.split()

        # Process lines based on the current section
        if section == "POL_RAD" and len(parts) > 5 and parts[0].isdigit():
            parts[4] = f"{float(parts[4]) * pol_scale:.7E}"
            parts[5] = f"{float(parts[5]) * rad_scale:.7E}"
            line = "{:>4} {:>15} {:>15} {:>15} {:>15} {:>15} {:>5} {:>5}\n".format(
                *parts
            )

        elif section == "CHRG" and len(parts) > 3 and parts[0].isdigit():
            parts[-1] = f"{float(parts[-1]) * q_scale:.7E}"
            line = "{:>4} {:>10} {:>10} {:>17}\n".format(*parts)

        elif section == "PERM_DIP" and len(parts) > 4 and parts[0].isdigit():
            parts[-1] = f"{float(parts[-1]) * p_scale:.7E}"
            line = "{:>4} {:>8} {:>8} {:>8} {:>17}\n".format(*parts)

        new_lines.append(line)

    with open(output_qin_path, "w") as file:
        file.writelines(new_lines)


def dipole_extract_data(file_path):
    # Define the pattern to match the relevant lines
    pattern = re.compile(
        r"Water #, induced, permanent, and total moments\s+(\d+)\s+([\d.E-]+)\s+([\d.E-]+)\s+([\d.E-]+)"
    )

    # Initialize an empty list to store the extracted data
    extracted_data = []

    # Open the file for reading
    with open(file_path, "r") as file:
        # Read the file line by line
        for line in file:
            # Use regular expression to find matches
            match = pattern.match(line.strip())
            if match:
                # Extract the required data
                water_number = int(match.group(1))
                induced_moment = float(match.group(2))
                permanent_moment = float(match.group(3))
                total_moment = float(match.group(4))

                # Append the extracted data as a tuple to the list
                extracted_data.append(
                    (water_number, induced_moment, permanent_moment, total_moment)
                )

    return extracted_data


def extract_density(file_path):
    with open(file_path, "r") as file:
        # Flag to indicate if we are in the right section
        in_section = False

        for line in file:
            # Check if the section starts
            if "A V E R A G E S   O V E R" in line:
                in_section = True
            # Check if the section ends
            elif "R M S  F L U C T U A T I O N S" in line:
                in_section = False

            # Extract density if in the right section
            if in_section and "Density    =" in line:
                # Extract and return the density value
                return float(line.split("=")[1].strip())
