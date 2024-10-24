import math


def hh_bond_distance(angle_val, oh_bond):
    angle_rad = (angle_val / 2) / 180 * math.pi
    hh_bond = oh_bond * math.sin(angle_rad) * 2
    return hh_bond


angle_val_dict = {
    "POL3": 109.47,
    "DA/NDA/DC": 104.52,
}

oh_bond_dict = {
    "POL3": 1,
    "DA/NDA/DC": 0.9572,
}

hh_bond_dict = {
    "POL3": hh_bond_distance(angle_val_dict["POL3"], oh_bond_dict["POL3"]),
    "DA/NDA/DC": hh_bond_distance(
        angle_val_dict["DA/NDA/DC"], oh_bond_dict["DA/NDA/DC"]
    ),
}

bond_val_dict = {}
bond_val_dict["case_1"] = [0.974523, 1.53209]


time2nstlim = {
    "1step": "1",
    "10step": "10",
    "1fs": "100",
    "10ps": "10000",
    "20ps": "20000",
    "50ps": "50000",
    "100ps": "100000",
    "200ps": "200000",
    "1ns": "1000000",
}


def calculate_angle_val(oh_bond, hh_bond):
    # Calculate the angle in radians for half the angle
    half_angle_rad = math.asin(hh_bond / (2 * oh_bond))

    # Convert from radians to degrees and multiply by 2 to get the full angle
    angle_val = 2 * (half_angle_rad * (180 / math.pi))

    return angle_val


if __name__ == "__main__":
    oh_bond, hh_bond = bond_val_dict["case_1"]
    angle_val = calculate_angle_val(oh_bond, hh_bond)
    print(angle_val)
