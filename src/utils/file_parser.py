import re

from collections import defaultdict


def read_output_file(file_path):
    lines = []

    with open(file_path) as f:
        for line in f:
            line = line.strip()
            if line:
                lines.append(line)

    file_string = " ".join(lines)
    return file_string


def preprocess_step_string(file_string):
    # replace special keys
    # file_string = file_string.replace("1-4 NB", "OneToFour NB")
    # file_string = file_string.replace("1-4 EEL", "OneToFour EEL")

    # truncate string
    try:
        avg_ind = file_string.index("A V E R A G E S   O V E R")
        file_string = file_string[:avg_ind]
    except ValueError:
        pass
    return file_string


def parse_step_section(section_string: str):
    delimeters = ["=", ":"]

    # make sure there have spaces around delimeters in the string
    formated_string = ""
    for c in section_string:
        if c in delimeters:
            formated_string += f" {c} "
        else:
            formated_string += c

    results = {}
    terms = formated_string.split()
    last_end = -1
    for i, term in enumerate(terms):
        if term in delimeters:
            key = " ".join(terms[last_end + 1 : i])
            val = terms[i + 1] if i < len(terms) - 1 else "Nan"
            last_end = i + 1

            try:
                val = float(val)
            except ValueError:
                pass
            results[key] = val
    return results


def parse_step_info(file_path, keys=None):
    file_string = read_output_file(file_path)
    file_string = preprocess_step_string(file_string)

    nstep_inds = [m.start() for m in re.finditer("NSTEP", file_string)]
    parsed_dict = defaultdict(dict)
    for i, ind in enumerate(nstep_inds):
        end_ind = None if i == len(nstep_inds) - 1 else nstep_inds[i + 1]
        sub_file_string = file_string[ind:end_ind]
        deliminator = "------------------------------------------------------------------------------"
        del_ind = sub_file_string.index(deliminator)
        sub_file_string = sub_file_string[:del_ind]

        info_results = parse_step_section(sub_file_string)
        step = str(i + 1)
        for k, v in info_results.items():
            parsed_dict[step][k] = v

    return dict(parsed_dict)


if __name__ == "__main__":
    # f = "/home/yxwu/water_model/MD_data/case_1_2023-v08.wat2-gas_npt_pmemd/MD/1ns_output/case_1_pgm_512wat.npt.pmemd.out"
    f = "/home/yxwu/water_model/MD_data/case_1_2023-v08.wat10-wat_nve_pmemd/MD/1ns_output/case_1_pgm_512wat.nve.pmemd.out"
    res = parse_step_info(f, [])
    print(res["117431"])
