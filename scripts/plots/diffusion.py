import json

import pandas as pd
import ptitprince as pt
import matplotlib.pyplot as plt
from ebatt.visualization.structural import (
    save_fig,
    edit_plot,
    _set_axis,
)

colors = ["#C7CAFA", "#FFFFFF", "#F9554A"]
diff_dir = "/Users/rashatwi/Desktop/fls/files/analysis/diffusion"
main_dir = "/Users/rashatwi/Desktop/fls/files/md"
diff_folders = ["dme"] + ["tte"] + [f"fls{i}" for i in range(0, 19)]

fig, axs = plt.subplots(1, 3, figsize=(20, 18))
for idx, sys in enumerate(["0.25_li2s2", "0.25_li2s4", "0.25_li2s8"]):
    df = pd.DataFrame()
    for sol in diff_folders:
        f = open(f"{main_dir}/{sys}/{sol}/nvt/FW.json")
        data = json.load(f)
        diff_names = {}
        for i, j in enumerate(data["spec"]["num_mols_list"]):
            for k, v in data["spec"]["nmols"].items():
                if v == j:
                    if k == sol:
                        k = "solvent"
                    diff_names[i + 1] = k
        msd_int = pd.read_csv(f"{diff_dir}/{sys}/{sol}/msd_int.csv").reset_index()
        msd_int = msd_int[["type", "diff"]]
        msd_int["solvent"] = sol
        msd_int["type"] = msd_int["type"].replace(diff_names)
        df = pd.concat([df, msd_int], ignore_index=True)
    mask = df["type"].isin(["solvent", "li"])
    df = df.loc[mask]
    df["type"] = df["type"].replace({"solvent": "Solvent", "li": "Li"})
    df["diff"] = df["diff"] * 10**9
    pt.RainCloud(
        x="solvent",
        y="diff",
        data=df,
        hue="type",
        orient="h",
        width_viol=1,
        alpha=0.9,
        dodge=True,
        palette=["#777DF2", "#F9554A"],
        linecolor="black",
        width_box=0.25,
        box_showfliers=False, box_showmeans=False,
        ax=axs[idx],
    )
    _set_axis(axs[idx], labelsize=32)
    axs[idx].set_xlabel(
        "Diffusivity\n ({} {})".format(r"$\mathtt{1O^{-9}}$", r"$\mathtt{m^2/s}$"),
        fontsize=32,
        weight="bold",
        labelpad=10,
    )
    axs[idx].get_legend().remove()
    if idx == 0:
        axs[idx].set_yticklabels(
            [i.upper() for i in diff_folders], fontsize=32, rotation="horizontal"
        )
        axs[idx].set_ylabel("Solvent", fontsize=32, weight="bold", labelpad=20)
    else:
        axs[idx].set_yticklabels([])
        axs[idx].set_ylabel(None)
    axs[idx].set_xlim(-0.1, 2.1)
    axs[idx].tick_params(axis="y", which="minor", length=0)
    edit_plot(axs[idx], linewidth=3.5)
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(
    handles,
    labels,
    loc="upper center",
    bbox_to_anchor=(0.5, 1.05),
    frameon=False,
    prop={"weight": "bold", "size": 36},
    ncol=2,
)
save_fig(fig, working_dir="../../files/results", name="diff_dist.pdf")
