# This code is used to plot carbon isotope data

# INPUT: OH2 Table S3.csv
# OUTPUT: OH2 Figure 2.png

# >>>>>>>>>

# Import libraries
import sys
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Plot parameters
plt.rcParams["legend.loc"] = "best"
plt.rcParams.update({'font.size': 7})
plt.rcParams['scatter.edgecolors'] = "w"
plt.rcParams['scatter.marker'] = "o"
plt.rcParams["lines.linewidth"] = 0.5
plt.rcParams["patch.linewidth"] = 0.5
plt.rcParams["figure.figsize"] = (9, 4)
plt.rcParams["savefig.dpi"] = 600
plt.rcParams["savefig.bbox"] = "tight"

# Import data from CSV files
df = pd.read_csv(sys.path[0] + "/OH2 Table S3.csv", sep=",")
df["SampleName"] = df["SampleName"].str.replace("Exp", "")

# Calculate the d13C offset
d13C_CO2 = -4.89
df['d13COffset'] = df['d13C'] - d13C_CO2

df_tripleO = pd.read_csv(sys.path[0] + "/OH2 Table S2.csv")
df_tripleO["SampleName"] = df_tripleO["SampleName"].str.replace("Exp", "")

# Merge the TILDAS and KIEL data
dfMerged = pd.merge(df, df_tripleO, on="SampleName")

# Create Figure 2
fig, (ax1, ax2) = plt.subplots(1, 2, sharex=True)

# Subplot A

# Samples
samples_not_measured = "C|D|F|J|K"
df_included = df[~df['SampleName'].str.contains(samples_not_measured)]
df_excluded = df[df['SampleName'].str.contains(samples_not_measured)]

ax1.scatter(df_included['d13C'], df_included['d18O'],
            color="#38342F", marker='o', s=50, ec="k", zorder=0,
            label="measured for $\Delta^{\prime 17}$O")
ax1.scatter(df_excluded['d13C'], df_excluded['d18O'],
            color="#9C9A8E", marker='o', s=50, ec="k", zorder=-1,
            label="not measured for $\Delta^{\prime 17}$O")
ax1.errorbar(df['d13C'], df["d18O"],
                yerr=df["d18O_error"],
                xerr=df["d13C_error"],
                fmt="none", color="#cacaca", zorder=-2)
for i, txt in enumerate(df['SampleName']):
    ax1.annotate(txt, (df['d13C'][i], df['d18O'][i]),
                 ha="center", va="center", size=5, color="w")

# Refgas d13C
ax1.axvspan(d13C_CO2-0.02, d13C_CO2+0.02, color="#00A099", alpha=0.2, lw=0)
ax1.axvline(d13C_CO2, color="#00A099", ls="dashed", lw=0.8)
ax1.annotate("$\delta^{13}$C of the CO$_2$\nused for the experiments", (d13C_CO2+0.01, 7.5), (d13C_CO2+0.2, 7.8),
            arrowprops=dict(arrowstyle="fancy", color="#BCBBB2", lw=0.5,connectionstyle="arc3,rad=-0.3"),
            ha="left", va="center")

# Arrows showing possible fractionation effects
ax1.annotate("preferential\nprecipitation\nof light CO$_2$", (-5, 6.8), ha="right", va="bottom", color="k")
ax1.annotate("", xy=(-5, 6.7), xytext=(-5.2, 6.7),
             arrowprops=dict(arrowstyle="<-,head_length=0.7,head_width=0.5", color="k", lw=1.5))

ax1.annotate("loss\nof light CO$_2$", (-4.8, 6.8), ha="left", va="bottom", color="k")
ax1.annotate("", xy=(-4.6, 6.7), xytext=(-4.8, 6.7),
             arrowprops=dict(arrowstyle="->,head_length=0.7,head_width=0.5", color="k", lw=1.5))

ax1.text(0.02, 0.98, "A", size=14, ha="left", va="top",
         transform=ax1.transAxes, fontweight="bold")

ax1.set_xlabel("$\delta^{13}$C (‰, VPDB)")
ax1.set_ylabel("$\delta^{18}$O (‰, VSMOW)")


# Subplot B

ax2.scatter(dfMerged['d13C'], dfMerged['Dp17O_AC'],
            color="#38342F", marker='o', s=50, ec="k",
            label=dfMerged['SampleName'])
ax2.errorbar(dfMerged['d13C'], dfMerged["Dp17O_AC"],
             yerr=dfMerged["Dp17O_error"],
             xerr=dfMerged["d13C_error"],
             fmt="none", color="#cacaca", zorder=0)
for i, txt in enumerate(dfMerged['SampleName']):
    ax2.annotate(txt, (dfMerged['d13C'][i], dfMerged['Dp17O_AC'][i]),
                 ha="center", va="center", size=5, color="w")

ax2.set_xlabel("$\delta^{13}$C (‰, VPDB)")
ax2.set_ylabel("$\delta^{18}$O (‰, VSMOW)")

ax2.text(0.02, 0.98, "B", size=14, ha="left", va="top",
         transform=ax2.transAxes, fontweight="bold")

plt.savefig(sys.path[0] + "/OH2 Figure 2.png")
plt.close("all")


# Plot d18O from IRMS vs d18O from TILDAS
# fig, ax1 = plt.subplots(figsize=(4.5, 4.5))

# plt.scatter(dfMerged["d18O"], dfMerged["d18O_AC"],
#             marker="o", color="#38342F", ec="k")

# x = np.linspace(5.5, 6.5, 10)
# plt.plot(x, x, c = "k", ls="dashed", zorder = -1)
# ax1.text(6.4, 6.35, "1:1", ha="center", va="center", rotation=45)

# diff = dfMerged["d18O"] - dfMerged["d18O_AC"]
# print(f"The average difference between the IRMS and TILDAS is {np.mean(diff):.2f}‰")

# plt.savefig(sys.path[0] + "/OH2 Figure SX.png", dpi=300, bbox_inches='tight')
# plt.close()