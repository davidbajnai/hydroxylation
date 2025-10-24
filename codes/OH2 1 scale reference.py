"""
This script was used in the original publication
to scale the reference gases based on the accepted values
for IAEA-603 and NBS-18 from Wostbrock et al. (2020).

INPUT:
- OH2 Table S2.csv

OUTPUT:
- OH2 Figure S1.png
"""

# Import libraries
import os
import pandas as pd
import matplotlib.pyplot as plt
from functions import *

# Retrieve directory paths
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, '../data')
figures_dir = os.path.join(script_dir, '../figures')

# Plot parameters
plt.rcParams.update({'font.size': 7})
plt.rcParams['scatter.edgecolors'] = "k"
plt.rcParams['scatter.marker'] = "o"
plt.rcParams["lines.linewidth"] = 0.5
plt.rcParams["patch.linewidth"] = 0.5
plt.rcParams["figure.figsize"] = (9, 4)
plt.rcParams["savefig.dpi"] = 800
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams['savefig.transparent'] = False
plt.rcParams['mathtext.default'] = 'regular'


# Function to print info for scaled samples
def print_info(df, d18O_col, Dp17O_col, sample_name):
    gas_subset = df[df["SampleName"].str.contains(sample_name)].copy()
    d18O_mean = gas_subset[d18O_col].mean()
    d18O_std = gas_subset[d18O_col].std()
    Dp17O_mean = gas_subset[Dp17O_col].mean()
    Dp17O_std = gas_subset[Dp17O_col].std()
    N_gas = len(gas_subset)
    print(f"{sample_name}, N = {N_gas}, d18O = {d18O_mean:.3f}(±{d18O_std:.3f})‰, ∆'17O = {Dp17O_mean:.0f}(±{Dp17O_std:.0f}) ppm", end="")


# Read data from Table S1
df = pd.read_csv(os.path.join(data_dir, "OH2_Table_S2.csv"))
df = df[df["measurementPeriod"].str.contains("2023-09-20 to 2023-10-08")]
df = df[df["SampleName"].str.contains("light|heavy|NBS18|IAEA")]

print("\nNumber of replicates and standard deviations for the references:")
print_info(df, "d18O", "Dp17O", "light"); print("\t<--- unscaled")
print_info(df, "d18O", "Dp17O", "heavy"); print("\t<--- unscaled")
print_info(df, "d18O", "Dp17O", "IAEA603"); print("\t<--- unscaled")
print_info(df, "d18O", "Dp17O", "NBS18"); print("\t<--- unscaled")

# Scale values

# Measured CO2 values
IAEA603_d18O_measured = df[df["SampleName"].str.contains("IAEA603")]["d18O"].mean()
IAEA603_d17O_measured = df[df["SampleName"].str.contains("IAEA603")]["d17O"].mean()
print("\nMeasured values:")
print(f"IAEA-603: d18O = {IAEA603_d18O_measured:.3f}‰, d17O = {IAEA603_d17O_measured:.3f}‰")

NBS18_d18O_measured = df[df["SampleName"].str.contains("NBS18")]["d18O"].mean()
NBS18_d17O_measured = df[df["SampleName"].str.contains("NBS18")]["d17O"].mean()
print(f"NBS-18: d18O = {NBS18_d18O_measured:.3f}‰, d17O = {NBS18_d17O_measured:.3f}‰")

# Accepted CO2 values
IAEA603_d18O_accepted = (to_VSMOW(-2.37) + 1000) * 1.01025 - 1000
IAEA603_Dp17O_accepted = -147 # from Wostbrock et al. (2020)
IAEA603_d17O_accepted = unprime(IAEA603_Dp17O_accepted/1000 + 0.528 * prime(IAEA603_d18O_accepted))
print("\nAccepted values:")
print(f"IAEA-603: d18O = {IAEA603_d18O_accepted:.3f}‰, ∆'17O = {IAEA603_Dp17O_accepted:.0f} ppm")

NBS18_d18O_accepted = (to_VSMOW(-23.2) + 1000) * 1.01025 - 1000
NBS18_Dp17O_accepted = -100 # from Wostbrock et al. (2020)
NBS18_d17O_accepted = unprime(NBS18_Dp17O_accepted/1000 + 0.528 * prime(NBS18_d18O_accepted))
print(f"NBS-18: d18O = {NBS18_d18O_accepted:.3f}‰, ∆'17O = {NBS18_Dp17O_accepted:.0f} ppm")

# Calculate the scaling factors
slope_d18O = (NBS18_d18O_accepted - IAEA603_d18O_accepted) / (NBS18_d18O_measured - IAEA603_d18O_measured)
intercept_d18O = IAEA603_d18O_accepted - slope_d18O * IAEA603_d18O_measured
slope_d17O = (NBS18_d17O_accepted - IAEA603_d17O_accepted) / (NBS18_d17O_measured - IAEA603_d17O_measured)
intercept_d17O = IAEA603_d17O_accepted - slope_d17O * IAEA603_d17O_measured

# Scale the measured values
df["d18O_scaled"] = slope_d18O*df['d18O']+intercept_d18O
df["d17O_scaled"] = slope_d17O*df['d17O']+intercept_d17O
df["Dp17O_scaled"] = Dp17O(df["d17O_scaled"], df["d18O_scaled"])

# Calculate the mean and SD of the replicate measurements
df_gases = df.loc[:, ["SampleName", "d18O_scaled", "Dp17O_scaled"]]
df_gases = df_gases[df_gases["SampleName"].str.contains("light|heavy")]
df_gases_A = df_gases.groupby('SampleName').mean().reset_index()
df_gases_A = df_gases_A.rename(columns={'d18O_scaled': 'd18O_CO2', 'Dp17O_scaled': 'Dp17O_CO2'})
df_gases_SD = df_gases.groupby('SampleName').std().reset_index()
df_gases_SD = df_gases_SD.rename(columns={'d18O_scaled': 'd18O_SD', 'Dp17O_scaled': 'Dp17O_SD'})

print("\nScaled values:")
print(df_gases_A.merge(df_gases_SD, on='SampleName').round({"Dp17O_CO2": 0, "d18O_CO2": 3, "Dp17O_SD": 0, "d18O_SD": 3}))


# Create Figure S1
fig, (ax1, ax2) = plt.subplots(1, 2)

for ax in fig.get_axes():
    i = fig.get_axes().index(ax)
    ax.text(0.025, 0.975, chr(65 + i),
            size=14, weight="bold", ha="left", va="top",
            transform=ax.transAxes)

# Create new dataframes for IAEA-603 and NBS-18
df_IAEA603 = df[df["SampleName"].str.contains("IAEA603")]
df_NBS18 = df[df["SampleName"].str.contains("NBS18")]

# Plot IAEA-603 data
ax1.scatter(df_IAEA603["d18O"], df_IAEA603["Dp17O"],
            label="IAEA-603", c="#1455C0", marker="s")
ax1.errorbar(df_IAEA603["d18O"], df_IAEA603["Dp17O"],
             yerr=df_IAEA603["Dp17OError"], xerr=df_IAEA603["d18OError"],
             fmt="none", color="#cacaca", zorder=0)

ax1.set_title("IAEA-603")
ax1.set_xlabel("$\delta^{18}$O (‰, unscaled)")
ax1.set_ylabel("$\Delta\prime^{17}$O (ppm, unscaled)")

# Plot NBS-18 data
ax2.scatter(df_NBS18["d18O"], df_NBS18["Dp17O"],
            label="NBS-18", c="#EC0016", marker="o")
ax2.errorbar(df_NBS18["d18O"], df_NBS18["Dp17O"],
             yerr=df_NBS18["Dp17OError"], xerr=df_NBS18["d18OError"],
             fmt="none", color="#cacaca", zorder=0)

ax2.set_title("NBS-18")
ax2.set_xlabel("$\delta^{18}$O (‰, unscaled)")
ax2.set_ylabel("$\Delta\prime^{17}$O (ppm, unscaled)")

plt.savefig(os.path.join(figures_dir, "OH2_Figure_S1.png"))
plt.close("all")