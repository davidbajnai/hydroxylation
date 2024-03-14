# The following code is used to calculate KIE for CO2
# based on the *revised* Table 7 in Christensen et al. (2021)

# >>>>>>>>>

# Import libraries
import pandas as pd
import numpy as np


# Functions that make life easier

def a18OH(T=273.15+22, eq="Z20-X3LYP"):
    if (eq == "Z20-X3LYP"):
        e18_H2O_OH = (-4.4573 + (10.3255 * 10**3) / (T) + (-0.5976 * 10**6) / (T)**2)

    elif (eq == "Z20-MP2"):
        e18_H2O_OH = (-4.0771 + (9.8350 * 10**3) / (T) + (-0.8729 * 10**6) / (T)**2)

    elif (eq == "BH21_original"):
        e18_H2O_OH = -0.034 * (T-273.15) + 43.4

    elif (eq == "BH21"):
        e18_H2O_OH = -0.035 * (T-273.15) + 40.1

    return e18_H2O_OH / 1000 + 1


def B_from_a(a, A):
    return (A + 1000) / a - 1000


# The *revised* Table 7 in Christensen et al. (2021)
df = pd.DataFrame({
    "T": [22, 22, 5, 4, 21, 28, 28, 28, 27, 17],
    "d18O_CO2": [36.6, 11.1, 35.7, 41.5, 41.5, 41.7, 41.7, 41.7, 41.7, 41.5],
    "H2O": [-6.7, -11.9, -9.6, -7, -7, -0.5, -0.9, -0.5, 0.1, -5.7],
    "d18O_mineral": [6.8, -11.2, 3.4, 8.3, 10, 14.3, 13.5, 14.4, 14, 11.4]})

df["d18O_OH_effective"] = B_from_a(a18OH(273.15+df["T"], "BH21"), df["H2O"])

df["d18O_CO2_effective"] = (3/2 * (df["d18O_mineral"] - 1/3 * df["d18O_OH_effective"]))

df["KIE"] = (np.log((df["d18O_CO2_effective"]+1000) / (df["d18O_CO2"]+1000))*1000)

print(f'KIE_CO2 = {df["KIE"].mean():.0f}(±{df["KIE"].std():.0f})‰')