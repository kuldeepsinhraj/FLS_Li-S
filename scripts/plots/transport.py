import pandas as pd

from ebatt.visualization.structural import (
    edit_plot,
    save_fig,
    create_fig,
)

df = pd.read_csv("../../files/results/transport.csv")

fig, ax = create_fig(
    figsize=(12, 6.5),
    set_axis=True,
    axis="both",
    xlabel="Dynamical Property",
    ylabel="Electrolyte",
    fontsize=28,
    labelsize=28,
)
viscosity_bars = ax.barh(
    df["System"],
    df["Viscosity"] * 1000,
    color="#777DF2",
    edgecolor="black",
    linewidth=2,
    height=0.6,
    label="Viscosity (cP)",
    xerr=df["Viscosity_STD"] * 1000,
    capsize=5,
)
conductivity_bars = ax.barh(
    df["System"],
    -df["Conductivity"] * 10,
    color="#F9554A",
    edgecolor="black",
    linewidth=2,
    height=0.6,
    label="Conductivity (mS/cm)",
    xerr=df["Conductivity_STD"] * 10,
    capsize=5,
)
ax.scatter(
    df["Viscosity_Exp"],
    df["System"],
    color="black",
    marker="*",
    edgecolors="black",
    s=300,
    label="Experimental",
)
ax.scatter(
    -df["Conductivity_Exp"],
    df["System"],
    color="black",
    marker="*",
    s=300,
    label=None,
)

ax.grid(axis="x", linestyle="--", linewidth=1, zorder=10, color="gray")
ax.axvline(x=0, color="black", lw=2)
edit_plot(ax, loc="best", linewidth=2)
ax.set_xlim([-13, 3.05])
ax.set_xticks([-12, -10, -8, -6, -4, -2, 0, 2])
ax.set_yticklabels(df["Label"], fontsize=22)
legend = ax.legend(fontsize=18, loc="center left", bbox_to_anchor=(0.05, 0.36))
legend.get_frame().set_facecolor("white")
legend.get_frame().set_linewidth(0.0)
save_fig(fig, working_dir="../../files/results", name="transport.pdf")
