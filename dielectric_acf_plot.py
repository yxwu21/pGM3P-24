import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from tqdm import tqdm


def read_diele_tot(filename):
    """Reads the diele_tot file and extracts the second column (Mt)."""
    data = np.loadtxt(filename)
    Mt = data[:, 1]
    return Mt


def calculate_time_autocorrelation(Mt):
    """Calculates the time-autocorrelation function <Mt(t)Mt(0)>/<Mt(0)^2>."""
    N = len(Mt)
    Mt0_squared = np.mean(Mt[:1] ** 2)
    autocorrelation = []

    for t in tqdm(range(N), desc="Calculating Time Autocorrelation"):
        product_sum = np.mean(Mt[: N - t] * Mt[t:N])
        autocorrelation_value = product_sum / Mt0_squared
        autocorrelation.append(autocorrelation_value)

    return autocorrelation


def exponential_decay(x, a):
    """Exponential decay function y = exp(-a * x)."""
    return np.exp(-a * x)


def main():
    filename = "/home8/yxwu/pGM_water_model/MD_data-pGM/vital-cosmos-120_340/512/vital-cosmos-120_340-255/MD/3_Prod/fort.101"
    Mt = read_diele_tot(filename)

    autocorrelation = calculate_time_autocorrelation(Mt)
    np.savetxt("outputs/autocorrelation.txt", autocorrelation)
    autocorrelation = np.loadtxt("outputs/autocorrelation.txt")
    timesteps = np.arange(1, len(autocorrelation) + 1) * 1e-3
    print(timesteps)

    # Fit the exponential decay function to the data
    popt, pcov = curve_fit(exponential_decay, timesteps, autocorrelation)
    fitted_a = popt[0]

    # Generate fitted curve data
    fitted_curve = exponential_decay(timesteps, fitted_a)

    # Plotting
    plt.figure(figsize=(10, 6))
    plt.plot(timesteps, autocorrelation, "bo", label="Autocorrelation Data")
    plt.plot(
        timesteps, fitted_curve, "r-", label=f"Fitted Curve: y=exp(-{fitted_a:.4f} x)"
    )
    plt.xlabel("Timestep")
    plt.ylabel("<Mt(t)Mt(0)> / <Mt(0)^2>")
    plt.title("Time Autocorrelation of Mt and Exponential Fit")
    plt.legend()
    plt.grid(True)
    plt.savefig("time_autocorrelation_fit.png")
    print(
        f"Figure saved as time_autocorrelation_fit.png with fitted parameter a={fitted_a:.4f}"
    )


if __name__ == "__main__":
    main()
