import numpy as np
from scipy import stats

# Temperature data
T = np.array(
    [
        235.5,
        248.0,
        260.5,
        273.0,
        285.5,
        298.0,
        310.5,
        323.0,
        335.5,
        348.0,
        360.5,
        373.0,
        400.0,
    ]
)

# Cvib data
Cvib = np.array(
    [
        -0.2247,
        -0.1894,
        -0.1559,
        -0.1241,
        -0.0938,
        -0.0651,
        -0.0378,
        -0.0118,
        0.0130,
        0.0365,
        0.0590,
        0.0804,
        0.1235,
    ]
)

# Cni data
Cni = np.array(
    [
        -0.0001,
        -0.0003,
        -0.0007,
        -0.0014,
        -0.0027,
        -0.0048,
        -0.0079,
        -0.0125,
        -0.0190,
        -0.0276,
        -0.0387,
        -0.0527,
        -0.0940,
    ]
)

# Perform linear regression for Cvib
slope_Cvib, intercept_Cvib, r_value_Cvib, p_value_Cvib, std_err_Cvib = stats.linregress(
    T, Cvib
)

# Perform linear regression for Cni
slope_Cni, intercept_Cni, r_value_Cni, p_value_Cni, std_err_Cni = stats.linregress(
    T, Cni
)

print(f"Cvib(T) = {slope_Cvib:.6f} * T + {intercept_Cvib:.6f}")
print(f"R-squared for Cvib: {r_value_Cvib**2:.4f}")

print(f"\nCni(T) = {slope_Cni:.6f} * T + {intercept_Cni:.6f}")
print(f"R-squared for Cni: {r_value_Cni**2:.4f}")


# Function to calculate Cvib for a given temperature
def calculate_Cvib(T):
    return slope_Cvib * T + intercept_Cvib


# Function to calculate Cni for a given temperature
def calculate_Cni(T):
    return slope_Cni * T + intercept_Cni


# Example usage
example_temp = 300
print(f"\nFor T = {example_temp} K:")
print(f"Cvib = {calculate_Cvib(example_temp):.4f} kcal mol^-1")
print(f"Cni = {calculate_Cni(example_temp):.4f} kcal mol^-1")
