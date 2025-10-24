"""
This code is used to calculate the hydroxylation thetas for two scenarios

OUTPUT:
- OH2 Figure 5.png
"""

# Import libraries
import os
import numpy as np
import matplotlib.pyplot as plt
from functions import *

# Retrieve directory paths
script_dir = os.path.dirname(os.path.abspath(__file__))
data_dir = os.path.join(script_dir, '../data')
figures_dir = os.path.join(script_dir, '../figures')

# Plotting parameters
plt.rcParams.update({'font.size': 7})
plt.rcParams['scatter.edgecolors'] = "k"
plt.rcParams['scatter.marker'] = "o"
plt.rcParams['lines.markersize'] = 5
plt.rcParams["lines.linewidth"] = 0.5
plt.rcParams["patch.linewidth"] = 0.5
plt.rcParams["figure.figsize"] = (9, 4)
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

    # Daeron et al. (2019) - calcite
    # hydroxylation theta for seawater becomes 0.533(+0.001), for lakewater 0.541(+0.002)
    # return np.exp((17.57 * 1000 / T - 29.13) / 1000)


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


def d17Ocal(T, d17Ow):
    return a17cal(T) * (d17Ow+1000) - 1000


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

    return e18_H2O_OH / 1000 + 1


def a17OH(T=273.15+22, eq="Z20-X3LYP", theta=0.5296):
    return a18OH(T, eq)**theta


# Create Figure 5
fig, (ax1, ax2) = plt.subplots(1, 2)

for ax in fig.get_axes():
    i = fig.get_axes().index(ax)
    ax.text(0.975, 0.975, chr(65 + i),
            size=14, weight="bold", ha="right", va="top",
            transform=ax.transAxes)

# Subpolot A: lakewater

# Lakewater
d18O_water = -6
Dp17O_water = 35
d17O_water = d17O(d18O_water, Dp17O_water)

# CO2
d18O_CO2 = 41.5
Dp17O_CO2 = -200
d17O_CO2 = d17O(d18O_CO2, Dp17O_CO2)

# CO2 KIE
CO2_KIE_shift = -3
CO2_KIE_theta = (np.log((12+16+16)/(12+17+16)))/(np.log((12+16+16)/(12+18+16)))
d18O_CO2_KIE = d18O_CO2 + CO2_KIE_shift
Dp17O_CO2_KIE = apply_theta(d18O_CO2, Dp17O_CO2, shift_d18O=CO2_KIE_shift, theta=CO2_KIE_theta)

# Line between CO2 and CO2 KIE
ax1.text((prime(d18O_CO2) + prime(d18O_CO2_KIE))/2 + 5, (Dp17O_CO2 + Dp17O_CO2_KIE)/2,
         r"$\theta_{CO_2}^{KIE}$ = " + f"{CO2_KIE_theta:.3f}",
         ha="left", va="center", color="k",
         bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax1.annotate("",
             xy=(prime(d18O_CO2), Dp17O_CO2),
             xytext=(prime(d18O_CO2_KIE), Dp17O_CO2_KIE),
             arrowprops=dict(arrowstyle="<|-", color="k", lw=1))

# Equilibrium OH-
d18O_OH_eq = B_from_a(a18OH(), d18O_water)
d17O_OH_eq = B_from_a(a17OH(), d17O_water)
Dp17O_OH_eq = Dp17O(d17O_OH_eq, d18O_OH_eq)

# Effective OH-
d18O_OH_eff = B_from_a(a18OH(eq="BH21"), d18O_water)
d17O_OH_eff = B_from_a(a17OH(eq="BH21", theta=0.523), d17O_water)
Dp17O_OH_eff = Dp17O(d17O_OH_eff, d18O_OH_eff)

# Line between effective OH- equilibrium OH-
ax1.text((prime(d18O_OH_eff) + prime(d18O_OH_eq))/2+5, (Dp17O_OH_eff + Dp17O_OH_eq)/2,
         r"$\theta_{OH^-}^{KIE}$ = 0.514",
         ha="left", va="center", color="k",
         bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax1.annotate("", xy=(prime(d18O_OH_eff), Dp17O_OH_eff), xycoords='data',
             xytext=(prime(d18O_OH_eq), Dp17O_OH_eq), textcoords='data',
             arrowprops=dict(arrowstyle="-|>", color="k", lw=1))

# Carbonate precipitating in equilibrium
d18Occ = d18Ocal(22+273.15, d18O_water)
d17Occ = d17Ocal(22+273.15, d17O(d18O_water, Dp17O_water))
Dp17Occ = Dp17O(d17Occ, d18Occ)

# Hydroxylation endmember
mixdfKIE = mix_d17O(d18O_OH_eff, D17O_A=Dp17O_OH_eff,
                    d18O_B=d18O_CO2_KIE, D17O_B=Dp17O_CO2_KIE)
ax1.plot(prime(mixdfKIE["mix_d18O"]), mixdfKIE["mix_Dp17O"],
         color="k", lw=.5, ls=":", zorder=-10)
d18Occ_OHeff_KIE = mixdfKIE["mix_d18O"].iloc[67]
Dp17Occ_OHeff_KIE = mixdfKIE["mix_Dp17O"].iloc[67]

# calculate theta between equilibrium carbonate and KIE carbonates
theta_OHeff = calculate_theta(d18Occ, Dp17Occ, d18Occ_OHeff_KIE, Dp17Occ_OHeff_KIE).round(3)
print(f"\nHydroxylation theta with effective OH-: {theta_OHeff:.3f} (LAKEWATER)")

# Plot the points
ax1.scatter(prime(d18O_water), Dp17O_water, marker="*", fc="#1455C0", ec="k", s = 50, zorder=3, label=f"H$_2$O ({round(d18O_water, 1)}‰, {round(Dp17O_water)} ppm)")
ax1.scatter(prime(d18O_CO2), Dp17O_CO2, marker="D", fc="w", ec="k", zorder=3, label=f"CO$_2$ ({round(d18O_CO2, 1)}‰, {round(Dp17O_CO2)} ppm)")
ax1.scatter(prime(d18O_CO2_KIE), Dp17O_CO2_KIE, marker="D", fc="k", ec="k", zorder=3, label=f"CO$_2$ ({round(d18O_CO2, 1)}‰, {round(Dp17O_CO2)} ppm)")
ax1.scatter(prime(d18O_OH_eq), Dp17O_OH_eq, marker="s", fc="w", ec="k", zorder=3, label="OH$^{-}_{eq}$")
ax1.scatter(prime(d18O_OH_eff), Dp17O_OH_eff, marker="s", fc="k", ec="k", zorder=3, label="OH$^{-}_{effective}$")
ax1.scatter(prime(d18Occ), Dp17Occ, marker="o", fc="w", ec="k", zorder=3, label="calcite (eq)")
ax1.scatter(prime(d18Occ_OHeff_KIE), Dp17Occ_OHeff_KIE, marker="o", fc="k", ec="k", zorder=3, label="hydroxylation endmember \w OH$^{-}_{effective}$")

# Add labels
ax1.text(prime(d18O_CO2)+1, Dp17O_CO2-5,
         f"dissolved CO$_2$\n({round(d18O_CO2, 1)}‰, {round(Dp17O_CO2)} ppm)",
         ha="left", va="top", color="k")
ax1.text(prime(d18O_CO2_KIE)+1, Dp17O_CO2_KIE+5,
         "CO$_2$ + KIE",
         ha="left", va="bottom", color="k")
ax1.text(prime(d18O_water), Dp17O_water+10,
         f"lakewater\n({round(d18O_water, 1)}‰, {round(Dp17O_water)} ppm)",
         ha="center", va="bottom", color="#1455C0")
ax1.text(prime(d18O_OH_eq)+3, Dp17O_OH_eq,
        r"OH$_{eq}^{-}$",
        ha="left", va="center", color="k",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax1.text(prime(d18O_OH_eff)+3, Dp17O_OH_eff,
        r"OH$^{-}$ $\plus$ KIE",
        ha="left", va="center", color="k")
ax1.text(prime(d18Occ)-2, Dp17Occ+5,
         "equilibrium\ncarbonate",
         ha="right", va="bottom", color="k")
ax1.text(prime(d18Occ_OHeff_KIE), Dp17Occ_OHeff_KIE-10,
         "hydroxylation\nendmember",
         ha="center", va="top", color="k")

# Plot the line between the equilibirum carbonate and the hydroxylation enmember
ax1.text(prime(d18Occ_OHeff_KIE)+5, (Dp17Occ + Dp17Occ_OHeff_KIE)/2,
        r"$\theta^{}_{hydrox.}$ = " + f"{theta_OHeff:.3f}",
        ha="right", va="center", color="#EC0016",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax1.annotate("", xy=(prime(d18Occ), Dp17Occ), xycoords='data',
            xytext=(prime(d18Occ_OHeff_KIE), Dp17Occ_OHeff_KIE), textcoords='data',
            arrowprops=dict(arrowstyle="<|-", color="#EC0016", lw=2), zorder=-1)

# Axis properties
ax1.set_xlabel("$\delta\prime^{18}$O (‰, VSMOW)")
ax1.set_ylabel("$\Delta\prime^{17}$O (ppm)")
ax1.set_xlim(-55, 85)
ax1.set_ylim(-260, 260)


# Subplot B: seawater

# Seawater
d18O_water = 0
Dp17O_water = -11
d17O_water = d17O(d18O_water, Dp17O_water)

# CO2
d18O_CO2 = A_from_a(1.042077731,d18O_water)
d17O_CO2 = A_from_a(1.042077731**0.524573353, d17O_water)
Dp17O_CO2 = Dp17O(d17O_CO2, d18O_CO2)

# CO2 KIE
d18O_CO2_KIE = d18O_CO2 + CO2_KIE_shift
Dp17O_CO2_KIE = apply_theta(d18O_CO2, Dp17O_CO2, shift_d18O = CO2_KIE_shift, theta = CO2_KIE_theta)

# Line between CO2 and CO2 KIE
ax2.text((prime(d18O_CO2) + prime(d18O_CO2_KIE))/2 + 5, (Dp17O_CO2 + Dp17O_CO2_KIE)/2,
        r"$\theta_{CO_2}^{KIE}$ = " + f"{CO2_KIE_theta:.3f}",
        ha="left", va="center", color="k",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax2.annotate("", xy=(prime(d18O_CO2), Dp17O_CO2), xycoords='data',
            xytext=(prime(d18O_CO2_KIE), Dp17O_CO2_KIE), textcoords='data',
            arrowprops=dict(arrowstyle="<|-", color="k", lw=1))

# Equilibrium OH-
d18O_OH_eq = B_from_a(a18OH(), d18O_water)
d17O_OH_eq = B_from_a(a17OH(), d17O_water)
Dp17O_OH_eq = Dp17O(d17O_OH_eq, d18O_OH_eq)

# Effective OH-
d18O_OH_eff = B_from_a(a18OH(eq = "BH21"), d18O_water)
d17O_OH_eff = B_from_a(a17OH(eq = "BH21", theta = 0.523), d17O_water)
Dp17O_OH_eff = Dp17O(d17O_OH_eff, d18O_OH_eff)

# Line between effective OH- equilibrium OH-
ax2.text((prime(d18O_OH_eff) + prime(d18O_OH_eq))/2+5, (Dp17O_OH_eff + Dp17O_OH_eq)/2,
        r"$\theta_{OH^-}^{KIE}$ = 0.514",
        ha="left", va="center", color="k",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax2.annotate("", xy=(prime(d18O_OH_eff), Dp17O_OH_eff), xycoords='data',
            xytext=(prime(d18O_OH_eq), Dp17O_OH_eq), textcoords='data',
            arrowprops=dict(arrowstyle="-|>", color="k", lw=1))

# Carbonate precipitating in equilibrium
d18Occ = d18Ocal(22+273.15, d18O_water)
d17Occ = d17Ocal(22+273.15, d17O(d18O_water, Dp17O_water))
Dp17Occ = Dp17O(d17Occ, d18Occ)

# Hydroxylation endmember carbonate
mixdfKIE = mix_d17O(d18O_OH_eff, D17O_A=Dp17O_OH_eff, d18O_B=d18O_CO2_KIE, D17O_B=Dp17O_CO2_KIE)
ax2.plot(prime(mixdfKIE["mix_d18O"]), mixdfKIE["mix_Dp17O"],
        color="k", lw=.5, ls=":", zorder=-10)

d18Occ_OHeff_KIE = mixdfKIE["mix_d18O"].iloc[67]
Dp17Occ_OHeff_KIE = mixdfKIE["mix_Dp17O"].iloc[67]

# calculate theta between equilibrium carbonate and KIE carbonates
theta_OHeff = calculate_theta(d18Occ, Dp17Occ, d18Occ_OHeff_KIE, Dp17Occ_OHeff_KIE)
print(f"\nHydroxylation theta with effective OH-: {theta_OHeff:.3f} (SEAWATER)")

# Plot the points
ax2.scatter(prime(d18O_water), Dp17O_water, marker="*", fc="#1455C0", ec="k", s = 50, label=f"H$_2$O ({round(d18O_water, 1)}‰, {round(Dp17O_water)} ppm)")
ax2.scatter(prime(d18O_CO2), Dp17O_CO2, marker="D", fc="w", ec="k", zorder=3, label=f"CO$_2$ ({round(d18O_CO2, 1)}‰, {round(Dp17O_CO2)} ppm)")
ax2.scatter(prime(d18O_CO2_KIE), Dp17O_CO2_KIE, marker="D", fc="k", ec="k", zorder=3, label=f"CO$_2$ ({round(d18O_CO2, 1)}‰, {round(Dp17O_CO2)} ppm)")
ax2.scatter(prime(d18O_OH_eq), Dp17O_OH_eq, marker="s", fc="w", ec="k", zorder=3, label="OH$^{-}_{eq}$")
ax2.scatter(prime(d18O_OH_eff), Dp17O_OH_eff, marker="s", fc="k", ec="k", zorder=3, label="OH$^{-}_{effective}$")
ax2.scatter(prime(d18Occ), Dp17Occ, marker="o", fc="w", ec="k", zorder=3, label="calcite (eq)")
ax2.scatter(prime(d18Occ_OHeff_KIE), Dp17Occ_OHeff_KIE, marker="o", fc="k", ec="k", zorder=3, label="hydroxylation endmember \w OH$^{-}_{effective}$")

# Add labels
ax2.text(prime(d18O_CO2)+1, Dp17O_CO2-5,
         f"dissolved CO$_2$\n({round(d18O_CO2, 1)}‰, {round(Dp17O_CO2)} ppm)",
         ha="left", va="top", color="k")
ax2.text(prime(d18O_CO2_KIE)+1, Dp17O_CO2_KIE+5,
         "CO$_2$ + KIE",
         ha="left", va="bottom", color="k")
ax2.text(prime(d18O_water), Dp17O_water+10,
         f"seawater\n({round(d18O_water, 1)}‰, {round(Dp17O_water)} ppm)",
         ha="center", va="bottom", color="#1455C0")
ax2.text(prime(d18O_OH_eq)+3, Dp17O_OH_eq,
        r"OH$_{eq}^{-}$",
        ha="left", va="center", color="k",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax2.text(prime(d18O_OH_eff)+3, Dp17O_OH_eff,
        r"OH$^{-}$ $\plus$ KIE",
        ha="left", va="center", color="k")
ax2.text(prime(d18Occ)-2, Dp17Occ+5,
         "equilibrium\ncarbonate",
         ha="right", va="bottom", color="k")
ax2.text(prime(d18Occ_OHeff_KIE), Dp17Occ_OHeff_KIE-10,
         "hydroxylation\nendmember",
         ha="center", va="top", color="k")

# Plot the lines between the equilibirum carbonate and the hydroxylation endmember
ax2.text(prime(d18Occ_OHeff_KIE)+5, (Dp17Occ + Dp17Occ_OHeff_KIE)/2,
        r"$\theta^{}_{hydrox.}$ = " + f"{theta_OHeff:.3f}",
        ha="right", va="center", color="#EC0016",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax2.annotate("", xy=(prime(d18Occ), Dp17Occ), xycoords='data',
            xytext=(prime(d18Occ_OHeff_KIE), Dp17Occ_OHeff_KIE), textcoords='data',
            arrowprops=dict(arrowstyle="<|-", color="#EC0016", lw=2), zorder=-1)

# Axis properties
ax2.set_xlabel("$\delta\prime^{18}$O (‰, VSMOW)")
ax2.set_ylabel("$\Delta\prime^{17}$O (ppm)")
ax2.set_xlim(-55, 85)
ax2.set_ylim(-260, 260)

# save the figure
plt.savefig(os.path.join(figures_dir, "OH2_Figure_5.png"))
plt.close("all")