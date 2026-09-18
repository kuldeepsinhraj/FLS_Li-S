import numpy as np
import json
import pandas as pd
import seaborn as sns
import matplotlib.gridspec as grid_spec

from ebatt.config.solvents_config import PROPS, radical_cids
from ebatt.visualization.structural import (
    save_fig,
    _set_axis,
    create_fig,
    scatter_plot,
    edit_plot,
    dist_plot,
)

palette = sns.color_palette("Set1")

# read dataframe containing pubchem, structural, cosmo, and dft properties
df = pd.read_csv("../../files/results/solvents.csv")
df = df[~df["CID"].isin(radical_cids)]

# read cids of literature fluorinated solvents
with open("../../files/lit_cids.txt", "r") as f:
    lit_cids = json.loads(f.read())

# read cids of screened fluorinated solvents
with open("../../files/cids_s2.txt", "r") as f:
    screened_cids = json.loads(f.read())

# name the properties to use for labels
prop_names = {
    "MolecularWeight": ["Molecular Weight (g/mol)", "Structural"],
    "Complexity": ["Structure Complexity", "Structural"],
    "R": ["Radius (Å)", "Structural"],
    "GyrationR": ["Gyration Radius (Å)", "Structural"],
    "vdWV": ["Van der Waals Volume (Å³)", "Structural"],
    "ProximityO": ["Proximity of Fluorine to Oxygen (Å)", "Structural"],
    "FDeg": ["Degree of Fluorination (%)", "Structural"],
    "ShapeInd": ["Shape Index", "Structural"],
    "StericEffect": ["Steric Effect", "Structural"],
    "ElectroCount": ["Electronegativity Count", "Structural"],
    "B1": ["B1 (Å)", "Structural"],
    "B5": ["B5 (Å)", "Structural"],
    "BP": ["Boiling Point (K)", "COSMO-RS"],
    "FP": ["Flash Point (K)", "COSMO-RS"],
    "Solubility_LiTFSI": ["Solubility of LiTFSI (mol/L)", "COSMO-RS"],
    "Solubility_LiFSI": ["Solubility of LiFSI (mol/L)", "COSMO-RS"],
    "BE_Salt": ["Binding Energy - LiTFSI (kcal/mol)", "DFT"],
    "BE_PS": ["Binding Energy - LiPS (kcal/mol)", "DFT"],
    "EW": ["Electochemical Window (eV)", "DFT"],
    "Polarizability": ["Polarizability (Å³)", "DFT"],
    "Dipole": ["Dipole Moment (Debye)", "DFT"],
}

# visualize distribution of calculated properties
fig, ax = create_fig(figsize=(24, 28), set_axis=False, axis="x")
n_groups = len(PROPS)
ncols = 3
nrows = int(np.ceil(n_groups / ncols))
gs = grid_spec.GridSpec(nrows, ncols)
fig.gca().axes.yaxis.set_visible(False)
ax_objs = []

for ind, prop in enumerate(prop_names.keys()):
    d = fig.add_subplot(gs[int(ind % nrows), int(np.floor(ind / nrows))])
    dist_plot(d, df[prop], bins=60, label="Property", color="#777DF2", alpha=0.6)
    d.grid(axis="x", linestyle="--", color="gray")
    d.set_yticklabels([])
    d.set_ylabel("")
    d.set_xlabel(prop_names[prop][0], fontsize=24, weight="bold")
    d.set_title(
        prop_names[prop][1],
        loc="right",
        pad=2,
        x=0.95,
        y=0.8,
        fontsize=24,
        backgroundcolor="#F2F2F2",
    )
    d.tick_params(left=False)
    # d.patch.set_alpha(0)
    _set_axis(d, "x", labelsize=24)
    dme = df[df["CID"] == 8071][prop].values
    dol = df[df["CID"] == 12586][prop].values

    # mark literature fls props in a different color
    lit_prop = df[df["CID"].isin(lit_cids)][prop].values
    scatter_plot(
        d,
        lit_prop,
        [d.get_ylim()[1] / 2] * len(lit_prop),
        color="#FFFF99",
        marker="s",
        size=200,
        label="Literature FLS",
    )

    # mark screened props in a different color
    screened_prop = df[df["CID"].isin(screened_cids)][prop].values
    scatter_plot(
        d,
        screened_prop,
        [d.get_ylim()[1] / 2] * len(screened_prop),
        color="#2BBF81",
        size=200,
        marker="^",
        label="Screened FLS",
    )

    # mark DME in a different color
    scatter_plot(
        d,
        dme,
        [d.get_ylim()[1] / 2] * len(dme),
        color="#F9554A",
        marker="o",
        size=200,
        label="DME",
    )

    # mark DOL in a different color
    scatter_plot(
        d,
        dol,
        [d.get_ylim()[1] / 2] * len(dol),
        color="#F2F2F2",
        marker="o",
        size=200,
        label="DOL",
    )

    edit_plot(d, linewidth=3)
d.legend(
    frameon=False,
    fontsize=28,
    bbox_to_anchor=(-1.9, 10.2, 2.5, -0.1),
    loc="lower left",
    ncol=6,
    mode="expand",
    borderaxespad=0.0,
)
ax.axis("off")
gs.update(
    hspace=0.5,
    wspace=0.2,
)
save_fig(fig, working_dir="../../files/results", name="dist_prop.pdf")
