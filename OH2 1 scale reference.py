# This code is used to scale the reference gases
# based on the "accepted" values for IAEA-603 and NBS-18
# from Wostbrock et al. (2020)

# INPUT: OH2 Table S1.csv
# OUTPUT: OH2 Figure S1.png

# >>>>>>>>>

# Import libraries
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Plot parameters
plt.rcParams["legend.loc"] = "best"
plt.rcParams.update({'font.size': 7})
plt.rcParams['scatter.edgecolors'] = "k"
plt.rcParams['scatter.marker'] = "o"
plt.rcParams["lines.linewidth"] = 0.5
plt.rcParams["patch.linewidth"] = 0.5
plt.rcParams["figure.figsize"] = (9, 4)
plt.rcParams["savefig.dpi"] = 600
plt.rcParams["savefig.bbox"] = "tight"

# Functions that make life easier
def prime(delta):
    dprime = 1000 * np.log(delta/1000 + 1)
    return dprime


def unprime(dprime):
    delta = (np.exp(dprime/1000) - 1) * 1000
    return delta


def Dp17O(d17O, d18O):
    # value in ppm
    return (prime(d17O) - 0.528 * prime(d18O)) * 1000


def to_VSMOW(VPDB):
    return (VPDB*1.03092)+30.92


# Read data from Table S1
df = pd.read_csv(sys.path[0] + "/" + "OH2 Table S1.csv")
df = df[df["MeasurementPeriod"].str.contains("2023-09-20 to 2023-10-08")]
df = df[df["SampleName"].str.contains("light|heavy|NBS18|IAEA")]

heavy_d18O_SD = df[df["SampleName"].str.contains("heavy")]["d18O"].std()
heavy_Dp17O_SD = df[df["SampleName"].str.contains("heavy")]["Dp17O"].std()
N_heavy = len(df[df["SampleName"].str.contains("heavy")])
print(f"\nHeavy reference gas, N = {N_heavy}, d18O SD = ±{heavy_d18O_SD:.3f}‰, ∆'17O = ±{heavy_Dp17O_SD:.0f} ppm")

light_d18O_SD = df[df["SampleName"].str.contains("light")]["d18O"].std()
light_Dp17O_SD = df[df["SampleName"].str.contains("light")]["Dp17O"].std()
N_light = len(df[df["SampleName"].str.contains("light")])
print(f"Light reference gas, N = {N_light}, d18O: ±{light_d18O_SD:.3f}‰, ∆'17O: ±{light_Dp17O_SD:.0f} ppm")

IAEA603_d18O_SD = df[df["SampleName"].str.contains("IAEA603")]["d18O"].std()
IAEA603_Dp17O_SD = df[df["SampleName"].str.contains("IAEA603")]["Dp17O"].std()
N_IAEA603 = len(df[df["SampleName"].str.contains("IAEA603")])
print(f"IAEA-603, N = {N_IAEA603}, d18O: ±{IAEA603_d18O_SD:.3f}‰, ∆'17O: ±{IAEA603_Dp17O_SD:.0f} ppm")

NBS18_d18O_SD = df[df["SampleName"].str.contains("NBS18")]["d18O"].std()
NBS18_Dp17O_SD = df[df["SampleName"].str.contains("NBS18")]["Dp17O"].std()
N_NBS18 = len(df[df["SampleName"].str.contains("NBS18")])
print(f"NBS-18, N = {N_NBS18}, d18O: ±{NBS18_d18O_SD:.3f}‰, ∆'17O: ±{NBS18_Dp17O_SD:.0f} ppm")

# Do the scaling here
IAEA603_d18O_measured = df[df["SampleName"].str.contains("IAEA603")]["d18O"].mean()
IAEA603_d17O_measured = df[df["SampleName"].str.contains("IAEA603")]["d17O"].mean()

NBS18_d18O_measured = df[df["SampleName"].str.contains("NBS18")]["d18O"].mean()
NBS18_d17O_measured = df[df["SampleName"].str.contains("NBS18")]["d17O"].mean()

# Accepted CO2 values
IAEA603_d18O_accepted = (to_VSMOW(-2.37) + 1000) * 1.01025 - 1000
IAEA603_Dp17O_accepted = -147 # from Wostbrock et al. (2020)
IAEA603_d17O_accepted = unprime(IAEA603_Dp17O_accepted/1000 + 0.528 * prime(IAEA603_d18O_accepted))
print(f"\nAccepted for IAEA603 d18O = {IAEA603_d18O_accepted:.3f}‰, ∆'17O = {IAEA603_Dp17O_accepted:.0f} ppm")

NBS18_d18O_accepted = (to_VSMOW(-23.2) + 1000) * 1.01025 - 1000
NBS18_Dp17O_accepted = -100 # from Wostbrock et al. (2020)
NBS18_d17O_accepted = unprime(NBS18_Dp17O_accepted/1000 + 0.528 * prime(NBS18_d18O_accepted))
print(f"Accepted for NBS18 d18O = {NBS18_d18O_accepted:.3f}‰, ∆'17O = {NBS18_Dp17O_accepted:.0f} ppm")

slope_d18O = (NBS18_d18O_accepted - IAEA603_d18O_accepted) / (NBS18_d18O_measured - IAEA603_d18O_measured)
intercept_d18O = IAEA603_d18O_accepted - slope_d18O * IAEA603_d18O_measured

slope_d17O = (NBS18_d17O_accepted - IAEA603_d17O_accepted) / (NBS18_d17O_measured - IAEA603_d17O_measured)
intercept_d17O = IAEA603_d17O_accepted - slope_d17O * IAEA603_d17O_measured

df["d18O_scaled"] = slope_d18O*df['d18O']+intercept_d18O
df["d17O_scaled"] = slope_d17O*df['d17O']+intercept_d17O
df["Dp17O_scaled"] = (prime(df["d17O_scaled"]) - 0.528 * prime(df["d18O_scaled"]))*1000

# Calculate the mean values from the replicate measurements
df_scaled = df.loc[:, ["SampleName", "d18O_scaled", "Dp17O_scaled"]]
df_mean = df_scaled.groupby('SampleName').mean().reset_index()
df_mean = df_mean.rename(columns={'d18O_scaled': 'd18O_CO2', 'Dp17O_scaled': 'Dp17O_CO2'})

# Calculate the standard deviation from the replicate measurements
df_std = df_scaled.groupby('SampleName').std().reset_index()
df_std = df_std.rename(columns={'d18O_scaled': 'd18O_error',  'Dp17O_scaled': 'Dp17O_error'})

dfMerged = df_mean.merge(df_std, on='SampleName')
dfMerged['Replicates'] = df_scaled.groupby('SampleName').size().reset_index(name='counts')['counts']

print("\nScaled values:\n", dfMerged.round({"Dp17O_CO2": 0, "Dp17O_error": 0}).round(3))

# Create new dataframes for IAEA-603 and NBS-18
df_IAEA603 = df[df["SampleName"].str.contains("IAEA603")]
df_NBS18 = df[df["SampleName"].str.contains("NBS18")]


# Create Figure S1
fig, (ax1, ax2) = plt.subplots(1, 2)

# Plot IAEA-603 data
ax1.scatter(df_IAEA603["d18O"], df_IAEA603["Dp17O"],
            label="IAEA-603", c="#1455C0", marker="s")
ax1.errorbar(df_IAEA603["d18O"], df_IAEA603["Dp17O"],
             yerr=df_IAEA603["Dp17OError"], xerr=df_IAEA603["d18OError"],
             fmt="none", color="#cacaca", zorder=0)

ax1.text(0.98, 0.98, "A", size=14, ha="right", va="top",
         transform=ax1.transAxes, fontweight="bold")
ax1.set_title("IAEA-603")
ax1.set_xlabel("$\delta^{18}$O (‰, VSMOW)")
ax1.set_ylabel("$\Delta^{\prime 17}$O (ppm, unscaled)")

# Plot NBS-18 data
ax2.scatter(df_NBS18["d18O"], df_NBS18["Dp17O"],
            label="NBS-18", c="#EC0016", marker="o")
ax2.errorbar(df_NBS18["d18O"], df_NBS18["Dp17O"],
             yerr=df_NBS18["Dp17OError"], xerr=df_NBS18["d18OError"],
             fmt="none", color="#cacaca", zorder=0)

ax2.text(0.98, 0.98, "B", size=14, ha="right", va="top",
         transform=ax2.transAxes, fontweight="bold")
ax2.set_title("NBS-18")
ax2.set_xlabel("$\delta^{18}$O (‰, VSMOW)")
ax2.set_ylabel("$\Delta^{\prime 17}$O (ppm, unscaled)")

plt.savefig(sys.path[0] + "/" + "OH2 Figure S1.png")
plt.close("all")