import json
import pandas as pd
import seaborn as sns
from ebatt.config.solvents_config import RULES, radical_cids
from ebatt.visualization.structural import _set_axis

from ebatt.visualization.structural import (
    save_fig,
    create_fig,
    scatter_plot,
)

palette = sns.color_palette("Set1")

# exclude cids of radicals from screening
df = pd.read_csv("../files/results/solvents.csv")
df = df[~df["CID"].isin(radical_cids)]
df["BE_Salt_R"] = df["BE_Salt"] - df.loc[df["CID"] == 8071]["BE_Salt"].values[0]
df["BE_PS_R"] = df["BE_PS"] - df.loc[df["CID"] == 8071]["BE_PS"].values[0]

# plot 2D BE
fig, ax = create_fig(
    figsize=(8, 6),
    set_axis=True,
    axis="both",
    xlabel="Relative BE - Salt (kcal/mol)",
    ylabel="Relative BE - PS (kcal/mol)",
)
scatter_plot(ax, df["BE_Salt_R"], df["BE_PS_R"], color=palette.as_hex()[3], marker="o")
ax.axvline(0, color="k", linestyle="--")
ax.axhline(0, color="k", linestyle="--")
ax.set_ylim(-20)
save_fig(fig, working_dir="../files/results", name="relative_be.png")

# plot 2D ionization potential
fig, ax = create_fig(
    figsize=(8, 6),
    set_axis=True,
    axis="both",
    xlabel="Fluorination Degree (%)",
    ylabel="IP (eV)",
)
ax.set_ylabel("IP (eV)", fontsize=18, color=palette.as_hex()[3])
ax2 = ax.twinx()
_set_axis(ax2)
scatter_plot(ax, df["FDeg"], df["IP"], color=palette.as_hex()[3], marker="o")
scatter_plot(ax2, df["FDeg"], df["BE_Salt"], color=palette.as_hex()[2], marker="o")
ax2.set_ylabel("BE - Salt (kcal/mol)", fontsize=18, color=palette.as_hex()[2])
save_fig(fig, working_dir="../files/results", name="corr.png")

# remove DME and DOL from df for screening
df = df[(df["CID"] != 8071) & (df["CID"] != 12586)]

# find solvents that meet the quantitative limits
qry = " and ".join(["{} >= {}".format(k, v) for k, v in RULES.items()])
df_s1 = df.query(qry).reset_index(drop=True)

# write cids to file
with open("../files/cids_s1.txt", "w") as f:
    f.write(json.dumps(df_s1["CID"].to_list()))

# select solvents with lower affinity towards PS compared to DME but
# higher or equal affinity towards salt
df_s2 = df_s1[(df_s1["BE_PS_R"] >= -0.5) & (df_s1["BE_Salt_R"] <= 0.5)].sort_values(
    "BE_PS_R", ascending=False
)

# write cids to file
with open("../files/cids_s2.txt", "w") as f:
    f.write(json.dumps(df_s2["CID"].to_list()))
