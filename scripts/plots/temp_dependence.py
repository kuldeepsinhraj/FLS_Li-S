import json
import numpy as np
import pandas as pd

from ebatt.config.solvents_config import radical_cids

from ebatt.visualization.structural import (
    _set_axis,
    edit_plot,
    save_fig,
    create_fig,
    line_plot,
)

df = pd.read_csv("../../files/results/temp_dependence.csv")

# Density
density_error_fls0 = np.round(
    np.mean(
        np.abs(df["Exp_Density"][0:7] - df["Comp_Density"][0:7])
        * 100
        / df["Exp_Density"][0:7]
    ),
    1,
)
density_error_fls18 = np.round(
    np.mean(
        np.abs(df["Exp_Density"][7:] - df["Comp_Density"][7:])
        * 100
        / df["Exp_Density"][7:]
    ),
    1,
)

fig, ax = create_fig(
    figsize=(10, 6.2),
    set_axis=True,
    axis="both",
    xlabel="Temperature (°C)",
    ylabel="Density\n(g/cm$^3$)",
    fontsize=28,
    labelsize=28,
    title=f"MAE = FLS0: {density_error_fls0}%, FLS18: {density_error_fls18}%",
)

line_plot(
    ax,
    df["Temperature"][0:7],
    df["Exp_Density"][0:7],
    fmt="-o",
    markerfacecolor="white",
    markeredgecolor="#F9554A",  # 777DF2
    color="black",
    markersize=10,
    # label="Exp",
)
line_plot(
    ax,
    df["Temperature"][0:7],
    df["Comp_Density"][0:7],
    fmt="-o",
    markerfacecolor="#F9554A",
    markeredgecolor="black",
    color="black",
    markersize=10,
    label="FLS0",
    show_errorbars=False,
    yerr=df["Density_STD"][0:7],
)


line_plot(
    ax,
    df["Temperature"][7:],
    df["Exp_Density"][7:],
    fmt="-o",
    markerfacecolor="white",
    markeredgecolor="#777DF2",  # 777DF2
    color="black",
    markersize=10,
    # label="Exp",
)
line_plot(
    ax,
    df["Temperature"][7:],
    df["Comp_Density"][7:],
    fmt="-o",
    markerfacecolor="#777DF2",
    markeredgecolor="black",
    color="black",
    markersize=10,
    label="FLS18",
    show_errorbars=False,
    yerr=df["Comp_Viscosity_STD"][0:7] * 1000,
)
ax.set_xticks([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_ylim(0.95, 1.45)
# ax.set_yticks([1, 2, 3])
ax.set_xlim(3, 42)
edit_plot(ax, legend=None, loc="lower left", linewidth=2.5)
save_fig(fig, working_dir="../../files/results", name="density.pdf")

# Viscosity
viscosity_error_fls0 = np.round(
    np.mean(
        np.abs(df["Exp_Viscosity"][0:7] - df["Comp_Viscosity"][0:7] * 1000)
        * 100
        / df["Exp_Viscosity"][0:7]
    ),
    1,
)
viscosity_error_fls18 = np.round(
    np.mean(
        np.abs(df["Exp_Viscosity"][7:] - df["Comp_Viscosity"][7:] * 1000)
        * 100
        / df["Exp_Viscosity"][7:]
    ),
    1,
)
fig, ax = create_fig(
    figsize=(10, 6.2),
    set_axis=True,
    axis="both",
    xlabel="Temperature (°C)",
    ylabel="Viscosity\n(mPa.s)",
    fontsize=28,
    labelsize=28,
    title=f"MAE = FLS0: {viscosity_error_fls0}%, FLS18: {viscosity_error_fls18}%",
)

line_plot(
    ax,
    df["Temperature"][0:7],
    df["Exp_Viscosity"][0:7],
    fmt="-o",
    markerfacecolor="white",
    markeredgecolor="#F9554A",  # 777DF2
    color="black",
    markersize=10,
    # label="Exp",
)
line_plot(
    ax,
    df["Temperature"][0:7],
    df["Comp_Viscosity"][0:7] * 1000,
    fmt="-o",
    markerfacecolor="#F9554A",
    markeredgecolor="black",
    color="black",
    markersize=10,
    label="FLS0",
    show_errorbars=False,
    yerr=df["Comp_Viscosity_STD"][0:7] * 1000,
)


line_plot(
    ax,
    df["Temperature"][7:],
    df["Exp_Viscosity"][7:],
    fmt="-o",
    markerfacecolor="white",
    markeredgecolor="#777DF2",  # 777DF2
    color="black",
    markersize=10,
    # label="Exp",
)
line_plot(
    ax,
    df["Temperature"][7:],
    df["Comp_Viscosity"][7:] * 1000,
    fmt="-o",
    markerfacecolor="#777DF2",
    markeredgecolor="black",
    color="black",
    markersize=10,
    label="FLS18",
    show_errorbars=False,
    yerr=df["Comp_Viscosity_STD"][0:7] * 1000,
)
ax.set_xticks([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_ylim(0.8, 3.5)
ax.set_yticks([1, 2, 3])
ax.set_xlim(3, 42)
edit_plot(ax, legend=[None, "FLS0", None, "FLS18"], loc="best", linewidth=2.5)
save_fig(fig, working_dir="../../files/results", name="viscosity.pdf")


# Conductivity
conductivity_error_fls0 = np.round(
    np.mean(
        np.abs(df["Exp_Conductivity"][0:7] - df["Comp_Conductivity"][0:7] * 10)
        * 100
        / df["Exp_Conductivity"][0:7]
    ),
    1,
)
conductivity_error_fls18 = np.round(
    np.mean(
        np.abs(df["Exp_Conductivity"][7:] - df["Comp_Conductivity"][7:] * 10)
        * 100
        / df["Exp_Conductivity"][7:]
    ),
    1,
)
fig, ax = create_fig(
    figsize=(10, 6.2),
    set_axis=True,
    axis="both",
    xlabel="Temperature (°C)",
    ylabel="Ionic Conductivity\n(mS/cm)",
    fontsize=28,
    labelsize=28,
    title=f"MAE = FLS0: {conductivity_error_fls0}%, FLS18: {conductivity_error_fls18}%",
)

line_plot(
    ax,
    df["Temperature"][0:7],
    df["Exp_Conductivity"][0:7],
    fmt="-o",
    markerfacecolor="white",
    markeredgecolor="#F9554A",  # 777DF2
    color="black",
    markersize=10,
    # label="Exp",
)
line_plot(
    ax,
    df["Temperature"][0:7],
    df["Comp_Conductivity"][0:7] * 10,
    fmt="-o",
    markerfacecolor="#F9554A",
    markeredgecolor="black",
    color="black",
    markersize=10,
    label="FLS0",
    show_errorbars=False,
    yerr=df["Comp_Viscosity_STD"][0:7] * 1000,
)


line_plot(
    ax,
    df["Temperature"][7:],
    df["Exp_Conductivity"][7:],
    fmt="-o",
    markerfacecolor="white",
    markeredgecolor="#777DF2",  # 777DF2
    color="black",
    markersize=10,
    # label="Exp",
)
line_plot(
    ax,
    df["Temperature"][7:],
    df["Comp_Conductivity"][7:] * 10,
    fmt="-o",
    markerfacecolor="#777DF2",
    markeredgecolor="black",
    color="black",
    markersize=10,
    label="FLS18",
    show_errorbars=False,
    yerr=df["Comp_Viscosity_STD"][0:7] * 1000,
)
ax.set_xticks([5, 10, 15, 20, 25, 30, 35, 40])
ax.set_ylim(5, 13)
# ax.set_yticks([1, 2, 3])
ax.set_xlim(3, 42)
edit_plot(ax, legend=None, loc="best", linewidth=2.5)
save_fig(fig, working_dir="../../files/results", name="conductivity.pdf")
