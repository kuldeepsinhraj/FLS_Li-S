import pandas as pd
import numpy as np

from ebatt.visualization.structural import (
    edit_plot,
    save_fig,
    scatter_plot,
    line_plot,
    create_fig,
)

dol_dme_xls = pd.ExcelFile("../../files/results/elec_dol_dme.xls")
dol_tte_xls = pd.ExcelFile("../../files/results/elec_dol_tte.xls")
dol_dme_fls0_xls = pd.ExcelFile("../../files/results/elec_dol_dme_fls0.xls")
dol_dme_fls18_xls = pd.ExcelFile("../../files/results/elec_dol_dme_fls18.xls")

dol_dme_df = pd.read_excel(dol_dme_xls, "Cycle", skiprows=1)
dol_tte_df = pd.read_excel(dol_tte_xls, "Cycle", skiprows=1)
dol_dme_fls0_df = pd.read_excel(dol_dme_fls0_xls, "Cycle", skiprows=1)
dol_dme_fls18_df = pd.read_excel(dol_dme_fls18_xls, "Cycle", skiprows=1)

dol_dme_vol = pd.read_excel(dol_dme_xls, "Record__1", skiprows=1)
dol_tte_vol = pd.read_excel(dol_tte_xls, "Record__1", skiprows=1)
dol_dme_fls0_vol = pd.read_excel(dol_dme_fls0_xls, "Record__1", skiprows=1)
dol_dme_fls18_vol = pd.read_excel(dol_dme_fls18_xls, "Record__1", skiprows=1)

dol_dme_vol = dol_dme_vol[dol_dme_vol["Cycle ID"] == 10]
dol_tte_vol = dol_tte_vol[dol_tte_vol["Cycle ID"] == 10]
dol_dme_fls0_vol = dol_dme_fls0_vol[dol_dme_fls0_vol["Cycle ID"] == 10]
dol_dme_fls18_vol = dol_dme_fls18_vol[dol_dme_fls18_vol["Cycle ID"] == 10]


# Discharge Capacity
fig, ax = create_fig(
    figsize=(12, 6.5),
    set_axis=True,
    axis="both",
    xlabel="Cycle Number",
    ylabel="Discharge\nCapacity (mAh/g)",
    fontsize=28,
    labelsize=28,
)

scatter_plot(
    ax,
    x=dol_dme_df["Cycle ID"],
    y=dol_dme_df["Cap_DChg"] * 1000,
    color="#777DF2",
    marker="o",
    size=50,
    label="DOL/DME",
)
scatter_plot(
    ax,
    x=dol_tte_df["Cycle ID"],
    y=dol_tte_df["Cap_DChg"] * 1000,
    color="#F9554A",
    marker="o",
    size=50,
    label="DOL/TTE",
)
scatter_plot(
    ax,
    x=dol_dme_fls0_df["Cycle ID"],
    y=dol_dme_fls0_df["Cap_DChg"] * 1000,
    color="#F2C57D",
    marker="o",
    size=50,
    label="DOL/DME/FLS0",
)
scatter_plot(
    ax,
    x=dol_dme_fls18_df["Cycle ID"],
    y=dol_dme_fls18_df["Cap_DChg"] * 1000,
    color="#7DF2D2",
    marker="o",
    size=50,
    label="DOL/DME/FLS18",
)

ax.set_xlim(-5, 205)
ax.set_ylim(10, 1100)
edit_plot(
    ax,
    loc="best",
    linewidth=2,
    legend=["DOL/DME", "DOL/TTE", "DOL/DME/FLS0", "DOL/DME/FLS18"],
)
save_fig(fig, working_dir="../../files/results", name="discharge_capacity.pdf")

# Coulombic Efficiency
fig, ax = create_fig(
    figsize=(12, 6.5),
    set_axis=True,
    axis="both",
    xlabel="Cycle Number",
    ylabel="Coulombic\nEfficiency (%)",
    fontsize=28,
    labelsize=28,
)

scatter_plot(
    ax,
    x=dol_dme_df["Cycle ID"],
    y=dol_dme_df["Unnamed: 6"],
    color="#777DF2",
    marker="o",
    size=50,
    label="DOL/DME",
)
scatter_plot(
    ax,
    x=dol_tte_df["Cycle ID"],
    y=dol_tte_df["Unnamed: 5"],
    color="#F9554A",
    marker="o",
    size=50,
    label="DOL/TTE",
)
scatter_plot(
    ax,
    x=dol_dme_fls0_df["Cycle ID"],
    y=dol_dme_fls0_df["Unnamed: 6"],
    color="#F2C57D",
    marker="o",
    size=50,
    label="DOL/DME/FLS0",
)
scatter_plot(
    ax,
    x=dol_dme_fls18_df["Cycle ID"],
    y=dol_dme_fls18_df["Unnamed: 5"],
    color="#7DF2D2",
    marker="o",
    size=50,
    label="DOL/DME/FLS18",
)

ax.set_xlim(-5, 205)
ax.set_ylim(85, 115)
edit_plot(
    ax,
    loc="best",
    linewidth=2,
    legend=["DOL/DME", "DOL/TTE", "DOL/DME/FLS0", "DOL/DME/FLS18"],
)
save_fig(fig, working_dir="../../files/results", name="coulombic_efficiency.pdf")

# Voltage Profiles
fig, ax = create_fig(
    figsize=(12, 6.5),
    set_axis=True,
    axis="both",
    xlabel="Capacity (mAh/g)",
    ylabel="Voltage\n($\mathbf{V}$ - vs. $\mathbf{Li}/\mathbf{Li}^+$)",
    fontsize=28,
    labelsize=28,
)

line_plot(
    ax,
    x=dol_dme_vol[dol_dme_vol["Step ID"] == 19]["CmpEng"],
    y=dol_dme_vol[dol_dme_vol["Step ID"] == 19]["Vol"],
    fmt="-",
    markerfacecolor=None,
    markeredgecolor=None,
    color="#777DF2",
    markersize=0,
    label="DOL/DME",
    linewidth=3.5,
)
line_plot(
    ax,
    x=dol_dme_vol[dol_dme_vol["Step ID"] == 20]["CmpEng"],
    y=dol_dme_vol[dol_dme_vol["Step ID"] == 20]["Vol"],
    fmt="-",
    markerfacecolor=None,
    markeredgecolor=None,
    color="#777DF2",
    markersize=0,
    label=None,
    linewidth=3.5,
)
line_plot(
    ax,
    x=dol_tte_vol[dol_tte_vol["Step ID"] == 19]["CmpEng"],
    y=dol_tte_vol[dol_tte_vol["Step ID"] == 19]["Vol"],
    fmt="-",
    markerfacecolor=None,
    markeredgecolor=None,
    color="#F9554A",
    markersize=0,
    label="DOL/TTE",
    linewidth=3.5,
)
line_plot(
    ax,
    x=dol_tte_vol[dol_tte_vol["Step ID"] == 20]["CmpEng"],
    y=dol_tte_vol[dol_tte_vol["Step ID"] == 20]["Vol"],
    fmt="-",
    markerfacecolor=None,
    markeredgecolor=None,
    color="#F9554A",
    markersize=0,
    label=None,
    linewidth=3.5,
)
line_plot(
    ax,
    x=dol_dme_fls0_vol[dol_dme_fls0_vol["Step ID"] == 19]["CmpEng"],
    y=dol_dme_fls0_vol[dol_dme_fls0_vol["Step ID"] == 19]["Vol"],
    fmt="-",
    markerfacecolor=None,
    markeredgecolor=None,
    color="#F2C57D",
    markersize=0,
    label="DOL/DME/FLS0",
    linewidth=3.5,
)
line_plot(
    ax,
    x=dol_dme_fls0_vol[dol_dme_fls0_vol["Step ID"] == 20]["CmpEng"],
    y=dol_dme_fls0_vol[dol_dme_fls0_vol["Step ID"] == 20]["Vol"],
    fmt="-",
    markerfacecolor=None,
    markeredgecolor=None,
    color="#F2C57D",
    markersize=0,
    label=None,
    linewidth=3.5,
)
line_plot(
    ax,
    x=dol_dme_fls18_vol[dol_dme_fls18_vol["Step ID"] == 19]["CmpEng"],
    y=dol_dme_fls18_vol[dol_dme_fls18_vol["Step ID"] == 19]["Vol"],
    fmt="-",
    markerfacecolor=None,
    markeredgecolor=None,
    color="#7DF2D2",
    markersize=0,
    label="DOL/DME/FLS18",
    linewidth=3.5,
)
line_plot(
    ax,
    x=dol_dme_fls18_vol[dol_dme_fls18_vol["Step ID"] == 20]["CmpEng"],
    y=dol_dme_fls18_vol[dol_dme_fls18_vol["Step ID"] == 20]["Vol"],
    fmt="-",
    markerfacecolor=None,
    markeredgecolor=None,
    color="#7DF2D2",
    markersize=0,
    label=None,
    linewidth=3.5,
)
ax.axvline(x=15, color="black", lw=2, linestyle="--")
ax.axvline(x=900, color="black", lw=2, linestyle="--")

ax.axvspan(-10, 40, color="gray", alpha=0.3)
ax.axvspan(870, 930, color="lightblue", alpha=0.3)

edit_plot(
    ax,
    loc="best",
    linewidth=2,
    legend=["DOL/DME", None, "DOL/TTE", None, "DOL/DME/FLS0", None, "DOL/DME/FLS18"],
)
legend = ax.legend(fontsize=18, bbox_to_anchor=(0.46, 0.57), frameon=False)
ax.set_xlim(-100, 2100)
ax.set_ylim(1.7, 2.9)
ax.set_yticks([i for i in [i for i in np.arange(1.8, 3.0, 0.2)]])
ax.yaxis.set_major_formatter("{:.1f}".format)
save_fig(fig, working_dir="../../files/results", name="voltage_profiles.pdf")

# Typical Voltages
solvents = [
    "DOL/DME/FLS0",
    "DOL/DME",
    "DOL/DME/FLS18",
    "DOL/TTE",
]
properties = [[2.6488, 2.112], [2.7, 2.086], [2.6801, np.nan], [2.6187, 2.054]]
properties_transposed = list(map(list, zip(*properties)))
fig, ax = create_fig(
    figsize=(12, 6.5),
    set_axis=True,
    axis="both",
    xlabel="",
    ylabel="Voltage\n($\mathbf{V}$ - vs. $\mathbf{Li}/\mathbf{Li}^+$)",
    fontsize=28,
    labelsize=28,
)
index = np.arange(len(solvents))
bar_width = 0.2
opacity = 1

bar1 = ax.bar(
    index,
    properties_transposed[0],
    bar_width,
    alpha=opacity,
    label="15 mAh/g",
    color="#777DF2",
    edgecolor="black",
)
bar2 = ax.bar(
    index + bar_width,
    properties_transposed[1],
    bar_width,
    alpha=opacity,
    label="900 mAh/g",
    color="#F9554A",
    edgecolor="black",
)

ax.set_xticks(index + bar_width / 2)
ax.set_xticklabels(solvents, rotation=0, fontsize=22)
ax.set_ylim(1.96, 2.84)
ax.grid(
    axis="y",
    linestyle="--",
    linewidth=1,
    zorder=10,
    color="gray",
    which="minor",
    alpha=0.2,
)
ax.grid(
    axis="y",
    linestyle="--",
    linewidth=1,
    zorder=10,
    color="gray",
    which="major",
    alpha=0.5,
)
ax.yaxis.set_major_formatter("{:.1f}".format)
edit_plot(ax, loc="best", linewidth=2, legend=["15 mAh/g", "900 mAh/g"])
save_fig(fig, working_dir="../../files/results", name="typical_voltages.pdf")

# EIS
