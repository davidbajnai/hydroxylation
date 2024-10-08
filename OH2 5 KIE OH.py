# This code is used to calculate the 18KIE for OH-

# INPUT: OH2 Table S3.csv, OH2 BH21 york error.csv
# OUTPUT: OH2 Figure 3.png

# >>>>>>>>>

# Import libraries
import os
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from functions import *

# Plot parameters
plt.rcParams.update({'font.size': 7})
plt.rcParams['scatter.edgecolors'] = "k"
plt.rcParams['scatter.marker'] = "o"
plt.rcParams['lines.markersize'] = 5
plt.rcParams["lines.linewidth"] = 0.5
plt.rcParams["patch.linewidth"] = 0.5
plt.rcParams["figure.figsize"] = (4, 4)
plt.rcParams["savefig.dpi"] = 800
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams['savefig.transparent'] = False
plt.rcParams['mathtext.default'] = 'regular'

# Functions that make life easier

def a18cal(T):
    # Hayles et al. (2018) - calcite
    B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    return np.exp(B_calcite) / np.exp(B_water)


def theta_cal(T):
    # Hayles et al. (2018) - calcite
    K_calcite = 1.019124E+09 / T**5 + -2.117501E+07 / T**4 + 1.686453E+05 / T**3 + -5.784679E+02 / T**2 + 1.489666E-01 / T + 0.5304852
    B_calcite = 7.027321E+14 / T**7 + -1.633009E+13 / T**6 + 1.463936E+11 / T**5 + -5.417531E+08 / T**4 + -4.495755E+05 / T**3  + 1.307870E+04 / T**2 + -5.393675E-01 / T + 1.331245E-04
    K_water = 7.625734E+06 / T**5 + 1.216102E+06 / T**4 + -2.135774E+04 / T**3 + 1.323782E+02 / T**2 + -4.931630E-01 / T + 0.5306551
    B_water = -6.705843E+15 / T**7 + 1.333519E+14 / T**6 + -1.114055E+12 / T**5 + 5.090782E+09 / T**4 + -1.353889E+07 / T**3 + 2.143196E+04 / T**2 + 5.689300 / T + -7.839005E-03
    a18 = np.exp(B_calcite) / np.exp(B_water)
    return K_calcite + (K_calcite-K_water) * (B_water / np.log(a18))


def a17cal(T):
    return a18cal(T)**theta_cal(T)


def d18Ocal(T, d18Ow):
    return a18cal(T) * (d18Ow+1000) - 1000


def d17Ocal(T, d18Ow):
    return a17cal(T) * (d18Ow+1000) - 1000


def a18OH(T=273.15+22, eq="Z20-X3LYP"):
    if (eq == "Z20-X3LYP"):
        e18_H2O_OH = (-4.4573 + (10.3255 * 10**3) /
                      (T) + (-0.5976 * 10**6) / (T)**2)
    elif (eq == "Z20-MP2"):
        e18_H2O_OH = (-4.0771 + (9.8350 * 10**3) /
                      (T) + (-0.8729 * 10**6) / (T)**2)
    elif (eq == "BH21_original"):
        e18_H2O_OH = -0.034 * (T-273.15) + 43.4
    elif (eq == "BH21"):
        e18_H2O_OH = -0.035 * (T-273.15) + 40.1
    elif (eq == "GZ19"):
        e18_H2O_OH = ((0.04264*(1000/T)+0.89702)-1)*1000

    return e18_H2O_OH / 1000 + 1


def a17OH(T = 273.15+22, eq = "Z20-X3LYP", theta = 0.530):
    return a18OH(T, eq)**theta


# Read in data
df = pd.read_csv(os.path.join(sys.path[0], "OH2 Table S3.csv"))
df_BH21_york_err = pd.read_csv(os.path.join(sys.path[0], "OH2 BH21 york error.csv"))

# Isotope composition of CO2 gas
d18O_CO2 = df.loc[df['SampleName'] == 'KoelnRefCO2-2', 'd18O_CO2'].iloc[0]
d18O_CO2_err = df.loc[df['SampleName'] == 'KoelnRefCO2-2', 'd18O_error'].iloc[0] #/ np.sqrt(df.loc[df['SampleName'] == 'KoelnRefCO2-2', 'Replicates'].iloc[0])
Dp17O_CO2 = df.loc[df['SampleName'] == 'KoelnRefCO2-2', 'Dp17O_CO2'].iloc[0]
Dp17O_CO2_err = df.loc[df['SampleName'] == 'KoelnRefCO2-2', 'Dp17O_error'].iloc[0]
print(f"\nThe composition of the tank CO2 is: d18O = {d18O_CO2:.2f}(±{d18O_CO2_err:.2f})‰, ∆'17O = {Dp17O_CO2:.0f}(±{Dp17O_CO2_err:.0f}) ppm")
df = df[~df['SampleName'].str.contains('Koeln')]

# Isotope composition of the water
d18O_water = -7.92
d18O_water_err = 0.71
Dp17O_water = 25
Dp17O_water_err = 5
print(f"\nThe composition of the water is: d18O = {d18O_water:.2f}(±{d18O_water_err:.2f})‰, ∆'17O = {Dp17O_water:.0f}(±{Dp17O_water_err:.0f}) ppm")


# Figure 3
fig, ax = plt.subplots()

temps = np.linspace(1, 80, 80)
model_Z20_X3LYP = 1000*np.log(a18OH(temps+273.15, "Z20-X3LYP"))
model_Z20_MP2 = 1000*np.log(a18OH(temps+273.15, "Z20-MP2"))
model_BH21 = 1000*np.log(a18OH(temps+273.15, "BH21"))
model_GZ19 = 1000*np.log(a18OH(temps+273.15, "GZ19"))

ax.plot(temps, model_Z20_X3LYP,
        c="k", ls=":", lw=1,
        label="Z20-X3LYP (theoretical)")
ax.plot(temps, model_Z20_MP2,
        c="k", ls="--", lw=1,
        label="Z20-MP2 (theoretical)")
ax.plot(temps, model_GZ19,
        c="k", ls="-.", lw=1,
        label="GZ19 (theoretical)")

ax.fill_between(df_BH21_york_err["x"], df_BH21_york_err["y_low"], df_BH21_york_err["y_up"],
                color='k', alpha=0.3, zorder=-10, ec='none')
ax.plot(temps, model_BH21,
        c="k", ls="-", lw=1,
        label=r"BH21$^\ddag$ (experimental)")

df["d18O_OH_effective"] = (df["d18O_AC"] - 2/3 * d18O_CO2) * 3
df["1000lna18"] = elena(d18O_water, df["d18O_OH_effective"])
ax.scatter(np.full(len(df), 22), df["1000lna18"],
           marker="o", color="#38342F", ec="k", zorder = 2,
           label="experimental data from\nthis study (" + r"$\it{n}$" + f" = {len(df)})")
print(f"\nThe mean value of 1000lna18OH- is {df['1000lna18'].mean():.1f}‰ (1σ = {df['1000lna18'].std():.1f}‰)")

# Show BH21 equation without the acid fractionation correction
ax.text(0, 14,
        r"$^\ddag $10$^3$ ln $\alpha_{H_2O/OH^-}^{18}$ = -0.034($\pm$0.004) $\times$ $\mathit{T}$ + 39.3($\pm$0.2)",
        ha="left", va="top", color="k")

# Calculate and show KIE for OH-
temp = 22
OH_eq = B_from_a(a18OH(temp+273.15, "Z20-X3LYP"), d18O_water)
OH_eff = df["d18O_OH_effective"].mean()

KIE_OH = (np.log((OH_eff+1000) / (OH_eq+1000))*1000)

print(f"\nKIE_OH = {KIE_OH:.1f}‰ (relative to the theoretical equilibrium value from Z20-X3LYP")
ax.annotate("", xy=(temp, np.log(a18OH(temp+273.15, "Z20-X3LYP"))*1000), xytext=(temp, df["1000lna18"].mean()),
            arrowprops=dict(arrowstyle="<|-", color="#EC0016", lw=1.5))

ax.text(23, 30, r"$^{18}KIE_{OH^{-}}$", 
        ha="left", va="center", color="#EC0016")

ax.set_ylim(12, 52)

ax.legend(loc = "upper right")

ax.set_ylabel(r"10$^3$ ln $\alpha_{H_2O/OH^-}^{18}$")
ax.set_xlabel("Temperature (°C)")

plt.savefig(os.path.join(sys.path[0], "OH2 Figure 3.png"))
plt.close("all")