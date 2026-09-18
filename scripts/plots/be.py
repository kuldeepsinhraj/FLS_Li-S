import json
import pandas as pd

from ebatt.config.solvents_config import radical_cids

from ebatt.visualization.structural import (
    _set_axis,
    edit_plot,
    save_fig,
    create_fig,
    scatter_plot,
)

# exclude cids of radicals from screening
df = pd.read_csv("../../files/results/solvents.csv")
df = df[~df["CID"].isin(radical_cids)]
df["BE_Salt_R"] = df["BE_Salt"] - df.loc[df["CID"] == 8071]["BE_Salt"].values[0]
df["BE_PS_R"] = df["BE_PS"] - df.loc[df["CID"] == 8071]["BE_PS"].values[0]

# read cids of literature fluorinated solvents
with open("../../files/lit_cids.txt", "r") as f:
    lit_cids = json.loads(f.read())

# read cids of screened fluorinated solvents
with open("../../files/cids_s2.txt", "r") as f:
    screened_cids = json.loads(f.read())

# plot 2D BE
fig, ax = create_fig(
    figsize=(8, 7),
    set_axis=True,
    axis="both",
    xlabel="Relative BE - LiTFSI (kcal/mol)",
    ylabel="Relative BE - LiPS (kcal/mol)",
    fontsize=28,
    labelsize=28,
)
ax.axvline(0, color="gray", linestyle="--")
ax.axhline(0, color="gray", linestyle="--")
scatter_plot(
    ax, df["BE_Salt_R"], df["BE_PS_R"], color="#777DF2", marker="o", label="All FLS"
)
# scatter_plot(
#     ax,
#     df[df["CID"].isin(screened_cids)]["BE_Salt_R"],
#     df[df["CID"].isin(screened_cids)]["BE_PS_R"],
#     color="#2BBF81",
#     size=200,
#     marker="^",
#     label="Screened FLS",
# )
scatter_plot(
    ax,
    df[df["CID"].isin(lit_cids)]["BE_Salt_R"],
    df[df["CID"].isin(lit_cids)]["BE_PS_R"],
    color="#FFFF99",
    marker="s",
    size=100,
    label="Literature FLS",
)
ax.set_xticks([-14, -7, 0, 7, 14])
ax.set_yticks([-14, -7, 0, 7, 14])
ax.set_xlim(-17, 17)
ax.set_ylim(-17, 17)
ax.annotate(
    "",
    xy=(0.0, 0.0),
    xycoords="data",
)
ax.annotate(
    "Desired \nQuadrant Space",
    xy=(2, 1),
    xytext=(-8.5, 10),
    fontsize=26,
    ha="center",
    color="#F9554A",
)
edit_plot(ax, linewidth=2.5)
ax.legend(
    frameon=False,
    fontsize=28,
    bbox_to_anchor=(0.58, 1.14),
    loc="upper center",
    ncol=2,
    mode="expand",
    borderaxespad=0.0,
)
save_fig(fig, working_dir="../../files/results", name="relative_be.pdf")

# plot FDeg vs BE/IP
fig, ax = create_fig(
    figsize=(9, 7),
    set_axis=True,
    axis="both",
    xlabel="Fluorination Degree (%)",
    ylabel="Ionization Potential (eV)",
    fontsize=28,
    labelsize=28,
    labelpad=20,
)
ax2 = ax.twinx()
_set_axis(ax2, labelsize=28)
scatter_plot(
    ax, df["FDeg"], df["IP"], color="#777DF2", marker="o", size=100, label="IP"
)
scatter_plot(ax2, df["FDeg"], df["BE_Salt"], color="#F2F2F2", marker="o", label="BE")
ax.set_yticks([4, 6, 8])
ax2.set_yticks([-40, -20, 0, 20, 40])
ax.set_xlim(-7, 107)
ax.set_ylim(3, 9)
ax2.set_ylim(-50, 30)
ax.set_xlabel("Fluorination Degree (%)", fontsize=28, weight="bold", labelpad=10)
ax2.set_ylabel("BE - LiTFSI (kcal/mol)", fontsize=28, weight="bold", labelpad=20)
edit_plot(ax, linewidth=2.5)
edit_plot(ax2, linewidth=2.5)
ax.legend(
    frameon=False,
    fontsize=28,
    bbox_to_anchor=(0.1, 1.12),
    loc="upper center",
    ncol=2,
    mode="expand",
    borderaxespad=0.0,
)
ax2.legend(
    frameon=False,
    fontsize=28,
    bbox_to_anchor=(0.5, 1.12),
    loc="upper center",
    ncol=2,
    mode="expand",
    borderaxespad=0.0,
)
save_fig(fig, working_dir="../../files/results", name="be.pdf")
