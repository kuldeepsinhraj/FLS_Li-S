import json

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as grid_spec

from ebatt.screening.core import s1_mol
from ebatt.config.solvents_config import PROPS, RULES
from ebatt.visualization.structural import (
    save_fig,
    dist_plot,
    edit_plot,
    create_fig,
    scatter_plot,
)

palette = sns.color_palette("Set1")

# read original solvents library and statistics
df = pd.read_csv("../files/results/solvents.csv")
df_stats = pd.read_csv("../files/results/stats.csv", index_col=0)

# screen original library using the deterministic approach
new_df, new_df_stats = s1_mol(
    df,
    df_stats,
    method="det",
    eps=1.5,
    step_size=0.01,
    desired_ratio=0.75,
    lib_size=430,
)

# read literature cids file
with open("../files/lit_cids.txt", "r") as f:
    lit_cids = json.loads(f.read())

# add missing literature cids to screened library
existing_cids = new_df.query(f"CID in {lit_cids}")["CID"].values
missing_cids = list(set(lit_cids).difference(existing_cids))
new_df = new_df.append([df[df["CID"] == i] for i in missing_cids])

new_df.to_csv("../files/results/s1_solvents.csv", index=False)
# TODO: update stats after adding missing CIDs
new_df_stats.to_csv("../files/results/s1_stats.csv")

# visualize distribution of calculated properties versus original properties
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="x")
n_groups = len(PROPS)
ncols = 4
nrows = int(np.ceil(n_groups / ncols))
gs = grid_spec.GridSpec(nrows, ncols)
fig.gca().axes.yaxis.set_visible(False)
ax_objs = []

for ind, prop in enumerate(df_stats.columns):
    d = fig.add_subplot(gs[int(ind % nrows), int(np.floor(ind / nrows))])
    dist_plot(d, df[prop], bins=60, label=prop, color=palette.as_hex()[1], alpha=0.2)
    dist_plot(
        d, new_df[prop], bins=60, label=prop, color=palette.as_hex()[2], alpha=0.2
    )
    # mark lit_props in a different color
    lit_prop = new_df[new_df["CID"].isin(lit_cids)][prop].values
    # 18413652, 12100158 (255), 19063496
    # selected_prop = [df.iloc[588][prop], df.iloc[756][prop], df.iloc[25][prop]]
    scatter_plot(
        d, lit_prop, [d.get_ylim()[1] / 2] * len(lit_prop), color=palette.as_hex()[3]
    )
    # scatter_plot(
    #     d,
    #     selected_prop,
    #     [d.get_ylim()[1] / 2] * len(selected_prop),
    #     color=palette.as_hex()[4],
    # )
    d.grid(axis="x", linestyle="--", color="gray")
    d.set_yticklabels([])
    d.set_ylabel("")
    d.tick_params(left=False)
    d.patch.set_alpha(0)
ax.axis("off")
edit_plot(fig, legend=["Original", "S1"], loc="lower right")
gs.update(hspace=0.5)
save_fig(fig, working_dir="../files/results", name="s1_dist_prop.png")
