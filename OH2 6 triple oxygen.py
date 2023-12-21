# This code creates Figures 3 and 4 of the manuscript

# INPUT: OH2 Table S1.csv
# OUTPUT: OH2 Figure 3.png, OH2 Figure 4.png

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
plt.rcParams['lines.markersize'] = 5
plt.rcParams["lines.linewidth"] = 0.5
plt.rcParams["patch.linewidth"] = 0.5
plt.rcParams["figure.figsize"] = (4, 4)
plt.rcParams["savefig.dpi"] = 600
plt.rcParams["savefig.bbox"] = "tight"

# Functions that make life easier


def prime(x):
    return 1000 * np.log(x / 1000 + 1)


def apply_prime_to_list(lst):
    new_lst = []
    for item in lst:
        result = prime(item)
        new_lst.append(result)
    return new_lst

def unprime(x):
    return (np.exp(x / 1000) - 1) * 1000


def Dp17O(d17O, d18O):
    return ((1000 * np.log(d17O / 1000 + 1)) - 0.528 * (1000 * np.log(d18O / 1000 + 1))) * 1000


def d17O(d18O, Dp17O):
    return unprime(Dp17O/1000 + 0.528 * prime(d18O))


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


def plot_calcite_equilibrium(Dp17Ow, d18Ow, Tmin, Tmax, ax, fluid_name="precipitating fluid", color="k"):

    d17Ow = unprime(0.528 * prime(d18Ow) + Dp17Ow/1000)

    ax.scatter(prime(d18Ow), Dp17O(d17Ow, d18Ow),
               marker="D", fc=color, ec="k", zorder=3, label=fluid_name)

    # equilibrium, entire T range
    toInf = np.arange(0, 330, 1) + 273.15
    d18O_mineral = d18Ocal(toInf, d18Ow)
    d17O_mineral = d17Ocal(toInf, d17Ow)
    mineral_equilibrium = np.array(
        [d18O_mineral, Dp17O(d17O_mineral, d18O_mineral), toInf]).T
    ax.plot(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
            ":", c=color, zorder=3)

    # equilibrium, highlight range
    equilibrium_temperatures = np.arange(Tmin, Tmax, 0.5) + 273.15
    colors = np.linspace(0, 1, len(equilibrium_temperatures))
    d18O_mineral = d18Ocal(equilibrium_temperatures, d18Ow)
    d17O_mineral = d17Ocal(equilibrium_temperatures, d17Ow)
    mineral_equilibrium = np.array([d18O_mineral, Dp17O(
        d17O_mineral, d18O_mineral), equilibrium_temperatures]).T
    ax.scatter(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
               marker=".", c=colors, cmap='coolwarm', zorder=3)

    # equilibrium, highlight range, marker every 10 °C
    equilibrium_temperatures = np.arange(Tmin, Tmax+1, 10) + 273.15
    d18O_mineral = d18Ocal(equilibrium_temperatures, d18Ow)
    d17O_mineral = d17Ocal(equilibrium_temperatures, d17Ow)
    mineral_equilibrium = np.array([d18O_mineral, Dp17O(
        d17O_mineral, d18O_mineral), equilibrium_temperatures]).T
    ax.scatter(prime(mineral_equilibrium[:, 0]), mineral_equilibrium[:, 1],
               s=15, marker="o", fc="white", ec=color, zorder=3, label="calcite equilibrium (" + str(Tmin) + "–" + str(Tmax) + " °C)")

    # Return equilibrium data as a dataframe
    equilibrium_df = pd.DataFrame(mineral_equilibrium)
    equilibrium_df[2] = equilibrium_df[2]-273.15
    equilibrium_df = equilibrium_df.rename(
        columns={0: 'd18O', 1: 'Dp17O', 2: 'temperature'})
    return equilibrium_df


def mix_d17O(d18O_A, d17O_A=None, D17O_A=None, d18O_B=None, d17O_B=None, D17O_B=None, step=100):
    ratio_B = np.arange(0, 1+1/step, 1/step)

    if d17O_A is None:
        d17O_A = unprime(D17O_A/1000 + 0.528 * prime(d18O_A))

    if d17O_B is None:
        d17O_B = unprime(D17O_B/1000 + 0.528 * prime(d18O_B))

    mix_d18O = ratio_B * float(d18O_B) + (1 - ratio_B) * float(d18O_A)
    mix_d17O = ratio_B * float(d17O_B) + (1 - ratio_B) * float(d17O_A)
    mix_D17O = Dp17O(mix_d17O, mix_d18O)
    xB = ratio_B * 100

    df = pd.DataFrame(
        {'mix_d17O': mix_d17O, 'mix_d18O': mix_d18O, 'mix_Dp17O': mix_D17O, 'xB': xB})
    return df


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


def a17OH(T = 273.15+22, eq = "Z20-X3LYP", theta = 0.530):
    return a18OH(T, eq)**theta


def B_from_a(a, A):
    return (A + 1000) / a - 1000


def A_from_a(a, B):
    return (B + 1000) * a - 1000


def epsilon(d18O_A, d18O_B):
    epsilon = ((d18O_A + 1000) / (d18O_B + 1000) - 1) * 1000
    return epsilon

def elena(d18O_A, d18O_B):
    elena = 1000*np.log((d18O_A + 1000) / (d18O_B + 1000))
    return elena

def calculate_theta(d18O_A, Dp17O_A, d18O_B, Dp17O_B):

    a18 = (d18O_B + 1000) / (d18O_A + 1000)
    a17 = (d17O(d18O_B, Dp17O_B) + 1000) / (d17O(d18O_A, Dp17O_A) + 1000)

    theta = round(np.log(a17) / np.log(a18), 4)

    return theta

def calculate_OH(d18O_CO2, Dp17O_CO2, d18O_precipitate, Dp17O_precipitate):
    d18O_OH = (d18O_precipitate - 2/3 * d18O_CO2) * 3
    d17O_OH = (d17O(d18O_precipitate, Dp17O_precipitate) - 2/3 * d17O(d18O_CO2, Dp17O_CO2)) * 3

    return d18O_OH, Dp17O(d17O_OH, d18O_OH)

monte_carlo_iterations = 10**3

# Read in TILDAS data
df = pd.read_csv(sys.path[0] + "/OH2 Table S2.csv", sep=",")

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

# Average of the precipitates
d18O_precipitate = df['d18O_AC'].mean()
d18O_precipitate_err = df['d18O_AC'].std()
Dp17O_precipitate = df['Dp17O_AC'].mean()
Dp17O_precipitate_err = df['Dp17O_AC'].std()
print(f"\nThe mean composition of the precipitates is: d18O = {d18O_precipitate:.2f}(±{d18O_precipitate_err:.2f})‰, ∆'17O = {Dp17O_precipitate:.0f}(±{Dp17O_precipitate_err:.0f}) ppm")

# Calculate the composition of the KIE OH- using Monte Carlo simulations
d18O_OH_lst = []
Dp17O_OH_lst = []
d18O_CO2_lst = []
d18O_precipitate_lst = []
Dp17O_CO2_lst = []
Dp17O_precipitate_lst = []
for _ in range(monte_carlo_iterations):

    sample_d18O_CO2 = np.random.normal(d18O_CO2, d18O_CO2_err)
    sample_Dp17O_CO2 = np.random.normal(Dp17O_CO2, Dp17O_CO2_err)
    sample_d18O_precipitate = np.random.normal(d18O_precipitate, d18O_precipitate_err)
    sample_Dp17O_precipitate = np.random.normal(Dp17O_precipitate, Dp17O_precipitate_err)
    
    result_d18O, result_Dp17O = calculate_OH(sample_d18O_CO2, sample_Dp17O_CO2, sample_d18O_precipitate, sample_Dp17O_precipitate)

    d18O_OH_lst.append(result_d18O)
    Dp17O_OH_lst.append(result_Dp17O)
    d18O_CO2_lst.append(sample_d18O_CO2)
    d18O_precipitate_lst.append(sample_d18O_precipitate)
    Dp17O_CO2_lst.append(sample_Dp17O_CO2)
    Dp17O_precipitate_lst.append(sample_Dp17O_precipitate)

d18O_OH = np.mean(d18O_OH_lst)
d18O_OH_err = np.std(d18O_OH_lst)
Dp17O_OH = np.mean(Dp17O_OH_lst)
Dp17O_OH_err = np.std(Dp17O_OH_lst)
print(f"\nThe calcualted composition of the OH- is: d18O = {d18O_OH:.2f}(±{d18O_OH_err:.2f})‰, Dp17O = {Dp17O_OH:.0f}(±{Dp17O_OH_err:.0f}) ppm")

# Calculate the effective H2O/OH- theta using Monte Carlo simulations
theta_effective_lst = []
for _ in range(monte_carlo_iterations):
    sample_d18O_OH = np.random.normal(d18O_OH, d18O_OH_err)
    sample_Dp17O_OH = np.random.normal(Dp17O_OH, Dp17O_OH_err)
    sample_d18O_water = np.random.normal(d18O_water, d18O_water_err)
    sample_Dp17O_water = np.random.normal(Dp17O_water, Dp17O_water_err)

    result = calculate_theta(sample_d18O_OH, sample_Dp17O_OH, sample_d18O_water, sample_Dp17O_water)
    
    theta_effective_lst.append(result)

theta_effective = np.round(np.mean(theta_effective_lst), 3)
theta_effective_err = np.round(np.std(theta_effective_lst), 3)
print(f"\nThe effective H2O/OH- theta is: {theta_effective}(±{theta_effective_err})")

# Calculate the OH-(KIE)/OH-(equilibrium) theta using Monte Carlo simulations
theta_eq_H2O_OH = 0.5296
e18_theoretical_model = "Z20-X3LYP"  # Z20-X3LYP or Z20-MP2
d18O_OH_eq = B_from_a(a18OH(T=273.15+22, eq=e18_theoretical_model), d18O_water)
Dp17O_OH_eq = Dp17O(B_from_a(a17OH(T=273.15+22, eq=e18_theoretical_model, theta=theta_eq_H2O_OH), d17O(d18O_water, Dp17O_water)), d18O_OH_eq)

theta_kinetic_lst = []
for _ in range(monte_carlo_iterations):
    sample_d18O_OH = np.random.normal(d18O_OH, d18O_OH_err)
    sample_Dp17O_OH = np.random.normal(Dp17O_OH, Dp17O_OH_err)
    sample_d18O_OH_eq = d18O_OH_eq
    sample_Dp17O_OH_eq = Dp17O(B_from_a(a17OH(T = 273.15+22, eq = e18_theoretical_model, theta = theta_eq_H2O_OH), d17O(d18O_water, Dp17O_water)), sample_d18O_OH_eq)

    result = calculate_theta(sample_d18O_OH, sample_Dp17O_OH, sample_d18O_OH_eq, sample_Dp17O_OH_eq)    
    theta_kinetic_lst.append(result)

theta_kinetic = np.round(np.mean(theta_kinetic_lst),3)
theta_kinetic_err = np.round(np.std(theta_kinetic_lst),3)
print(f"The kinetic OH-/OH- theta is: {theta_kinetic}(±{theta_kinetic_err})")


fig, ax = plt.subplots()

# Monte Carlo results
# ax.scatter(apply_prime_to_list(d18O_OH_lst), Dp17O_OH_lst,
#               marker=".", fc="#cacaca", alpha=0.1, zorder=-2)
# ax.scatter(apply_prime_to_list(d18O_CO2_lst), Dp17O_CO2_lst,
#               marker=".", fc="#cacaca", alpha=0.1, zorder=-2)
# ax.scatter(apply_prime_to_list(d18O_precipitate_lst), Dp17O_precipitate_lst,
#               marker=".", fc="#cacaca", alpha=0.1, zorder=-2)

# water
ax.scatter(prime(d18O_water), Dp17O_water,
           marker="*", fc="#1455C0", ec="k", s = 50, label="H$_2$O")
ax.errorbar(prime(d18O_water), Dp17O_water, xerr=d18O_water_err, yerr=Dp17O_water_err,
            fmt="None", ecolor="#1455C0", zorder=-1)
ax.text(prime(d18O_water)+2, Dp17O_water, "H$_2$O",
        ha="left", va="center", color="#1455C0")

# KIE OH-
ax.scatter(prime(d18O_OH), Dp17O_OH,
           marker="s", fc="#EC0016", ec="k", label="OH$^-$")
ax.errorbar(prime(d18O_OH), Dp17O_OH, xerr=d18O_OH_err, yerr=Dp17O_OH_err,
            fmt="None", ecolor="#EC0016", zorder=-1)
ax.text(prime(d18O_OH)+2, Dp17O_OH,
        r"OH$^{-}$ $\plus$ KIE",
        ha="left", va="center", color="#EC0016")

# Line between effective OH- and water
ax.text((prime(d18O_OH) + prime(d18O_water))/2+25, (Dp17O_OH + Dp17O_water)/2,
        r"$\theta_{H_2O/OH^-}^{effective}$ = " + f"{theta_effective}  \n(±{theta_effective_err})",
        ha="right", va="center", color="#814997")
ax.annotate("", xy=(prime(d18O_OH), Dp17O_OH), xycoords='data',
            xytext=(prime(d18O_water), Dp17O_water), textcoords='data',
            arrowprops=dict(arrowstyle="<|-|>", color="#814997", lw=1.5), zorder = -1)

# OH- equilibrium
ax.scatter(prime(d18O_OH_eq), Dp17O_OH_eq,
           marker="s", fc="#FFBB00", ec="k", label=r"OH$_{eq}^{-}$")
ax.text(prime(d18O_OH_eq)-2, Dp17O_OH_eq,
        r"OH$_{eq}^{-}$",
        ha="right", va="center", color="#FFBB00")
print(f"Equilibrium OH- d18O = {d18O_OH_eq:.2f}")
print(f"The 18KIE_OH- is: {(d18O_OH - d18O_OH_eq):.2f}")

# Line between equilibrium OH- and water
ax.text(0, -10,
        r"$\theta_{H_2O/OH^-}^{eq}$ = " + f"{theta_eq_H2O_OH}",
        ha="right", va="center", color="#63A615")
ax.annotate("", xy=(prime(d18O_OH_eq), Dp17O_OH_eq), xycoords='data',
            xytext=(prime(d18O_water), Dp17O_water), textcoords='data',
            arrowprops=dict(arrowstyle="<|-|>", color="#63A615", lw=1.5), zorder = -1)

# Line between effective OH- equilibrium OH-
grahams_law = np.round((np.log((16+1)/(17+1)))/(np.log((16+1)/(18+1))),3)
print(f"The KIE theta is: {theta_kinetic}. The expected value based on Graham's law is: {grahams_law}")
ax.text((prime(d18O_OH) + prime(d18O_OH_eq))/2+10, (Dp17O_OH + Dp17O_OH_eq)/2,
        r"$\theta_{OH^-}^{KIE}$ = " + f"{theta_kinetic}  \n(±{theta_kinetic_err})",
        ha="right", va="top", color="#FF7A00",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))
ax.annotate("", xy=(prime(d18O_OH), Dp17O_OH), xycoords='data',
            xytext=(prime(d18O_OH_eq), Dp17O_OH_eq), textcoords='data',
            arrowprops=dict(arrowstyle="-|>", color="#FF7A00", lw=1.5), zorder = -1)

# CO2
ax.scatter(prime(d18O_CO2), Dp17O_CO2,
           marker="D", fc="#00A099", ec="k", label="CO$_2$")
ax.errorbar(prime(d18O_CO2), Dp17O_CO2, xerr=d18O_CO2_err, yerr=Dp17O_CO2_err,
            fmt="None", ecolor="#00A099", elinewidth=0.5, zorder=-1)
ax.text(prime(d18O_CO2)-2, Dp17O_CO2,
        "CO$_2$",
        ha="right", va="center", color="#00A099")

# Precipitate
ax.scatter(prime(d18O_precipitate), Dp17O_precipitate,
           marker="o", c="#38342F", ec="k", label="precipitates")
ax.errorbar(prime(d18O_precipitate), Dp17O_precipitate, xerr=d18O_precipitate_err, yerr=Dp17O_precipitate_err,
            fmt="None", ecolor="#38342F", elinewidth=0.5, zorder=-1)
ax.text(prime(d18O_precipitate), (Dp17O_precipitate+20), "witherite\nprecipitates",
        ha="center", va="bottom", color="#38342F",
        bbox=dict(fc='white', ec="None", alpha=0.8, pad=0.1))

# BaCO3 in equilibrium
d18OBaCO3 = d18Ocal(22+273.15, d18O_water)
d17OBaCO3 = d17Ocal(22+273.15, d17O(d18O_water, Dp17O_water))
print(f"The difference between the precipitates and carbonate equilibrium is: {round(Dp17O_precipitate - Dp17O(d17OBaCO3, d18OBaCO3))}")
ax.scatter(prime(d18OBaCO3), Dp17O(d17OBaCO3, d18OBaCO3),
           marker="o", c="w", ec = "k", label="equilibrium carbonate")
ax.text(prime(d18OBaCO3), Dp17O(d17OBaCO3, d18OBaCO3)+10, "equilibrium\ncarbonate",
        ha="center", va="bottom", color="#878C96")

# Mixing curve between KIE OH- and CO2
mixdf = mix_d17O(d18O_OH, D17O_A=Dp17O_OH, d18O_B=d18O_CO2, D17O_B=Dp17O_CO2)
ax.plot(prime(mixdf["mix_d18O"]), mixdf["mix_Dp17O"],
        color="#282D37", lw=.5, ls=":", zorder=-2)

ax.set_ylabel("$\Delta^{\prime 17}$O (ppm)")
ax.set_xlabel("$\delta^{\prime 18}$O (‰, VSMOW)")

# Save figure
plt.savefig(sys.path[0] + "/OH2 Figure 4.png")
plt.close("all")


# Create abstract graphics
plt.rcParams["figure.figsize"] = (13/2, 5/2)
plt.rcParams.update({'font.size': 14})
plt.rcParams["lines.linewidth"] = 1  # error bar width
plt.rcParams["patch.linewidth"] = 1  # marker edge width

fig, ax = plt.subplots()

fig.patch.set_facecolor('#1455C0')
ax.set_facecolor('#1455C0')

# change axis colors to white
ax.spines['bottom'].set_color('w')
ax.spines['top'].set_color('w')
ax.spines['left'].set_color('w')
ax.spines['right'].set_color('w')
ax.xaxis.label.set_color('w')
ax.yaxis.label.set_color('w')
ax.tick_params(axis='x', colors='w')
ax.tick_params(axis='y', colors='w')

# water
ax.scatter(prime(d18O_water), Dp17O_water,
           fc="w", ec="#1455C0", marker="o", s=50, zorder=3)
ax.text(prime(d18O_water), Dp17O_water-15,
        "H$_2$O",
        ha="center", va="top", color="w")

# KIE OH-
ax.scatter(prime(d18O_OH), Dp17O_OH,
           fc="w", ec="#1455C0", marker="o", s=50, zorder=3)
ax.text(prime(d18O_OH), Dp17O_OH+15,
        r"OH$_{KIE}^{-}$",
        ha="center", va="bottom", color="w")

# Line between effective OH- and water
ax.text((prime(d18O_OH) + prime(d18O_water))/2+13, (Dp17O_OH + Dp17O_water)/2+20,
        r"$\theta_{H_2O/OH^-}^{effective}$ = " + f"{theta_effective}",
        ha="right", va="center", color="w",
        fontsize=12)
ax.annotate("", xy=(prime(d18O_OH), Dp17O_OH), xycoords='data',
            xytext=(prime(d18O_water), Dp17O_water), textcoords='data',
            arrowprops=dict(arrowstyle="<|-|>", color="w", lw=1.5), zorder = -1)

# OH- equilibrium
ax.scatter(prime(d18O_OH_eq), Dp17O_OH_eq,
           fc="w", ec="#1455C0", marker="o", s=50, zorder=3)
ax.text(prime(d18O_OH_eq), Dp17O_OH_eq-15,
        r"OH$_{eq.}^{-}$",
        ha="center", va="top", color="w")

# Line between equilibrium OH- and water
ax.text((prime(d18O_OH_eq) + prime(d18O_water))/2, (Dp17O_OH_eq + Dp17O_water)/2-20,
        r"$\theta_{H_2O/OH^-}^{equilibrium}$ = " + f"{theta_eq_H2O_OH}",
        ha="center", va="top", color="w",
        fontsize=12)
ax.annotate("", xy=(prime(d18O_OH_eq), Dp17O_OH_eq), xycoords='data',
            xytext=(prime(d18O_water), Dp17O_water), textcoords='data',
            arrowprops=dict(arrowstyle="<|-|>", color="w", lw=1.5), zorder = -1)

# Line between effective OH- equilibrium OH-
grahams_law = np.round((np.log((16+1)/(17+1)))/(np.log((16+1)/(18+1))),3)
ax.text((prime(d18O_OH) + prime(d18O_OH_eq))/2, (Dp17O_OH + Dp17O_OH_eq)/2-30,
        r"$\theta_{OH^-}^{KIE}$ = " + f"{theta_kinetic}",
        ha="right", va="center", color="w",
        fontsize=12)
ax.annotate("", xy=(prime(d18O_OH), Dp17O_OH), xycoords='data',
            xytext=(prime(d18O_OH_eq), Dp17O_OH_eq), textcoords='data',
            arrowprops=dict(arrowstyle="-|>", color="w", lw=1.5), zorder = -1)

# Axis parameters
ax.set_ylabel("$\Delta^{\prime 17}$O")
ax.set_xlabel("$\delta^{\prime 18}$O")
ax.set_ylim(-90, 290)
ax.set_xlim(-51, -6)
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.xaxis.set_ticks_position('none')
ax.yaxis.set_ticks_position('none')

# Save figure
plt.savefig(sys.path[0] + "/OH2 Graphical Abstract.png")
plt.close("all")