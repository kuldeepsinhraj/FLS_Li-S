import json
import time
import itertools

from itertools import combinations

import numpy as np
import seaborn as sns
import matplotlib.patches as patches
import matplotlib.gridspec as grid_spec

from matplotlib.ticker import FormatStrFormatter

from pymatgen.core.structure import Molecule

from ebatt.config.solvents_config import ABBREVIATIONS
from ebatt.descriptors.structural import similarity_measure
from ebatt.visualization.structural import (
    _set_axis,
    save_fig,
    dist_plot,
    ecdf_plot,
    create_fig,
    edit_plot,
)

palette = sns.color_palette("Set1")
iter_palette = itertools.cycle(palette)

# read literature cids file
with open("../files/lit_cids.txt", "r") as f:
    lit_cids = json.loads(f.read())

# read cids file
with open("../files/cids.txt", "r") as f:
    cids = json.loads(f.read())

molecules = {}
for cid in cids:
    try:
        molecules[cid] = Molecule.from_file(f"../files/xyz/optimized/{cid}.xyz")
    except Exception as e:
        print(e)

# calculate similarity measure between each literature fluorinated solvent and
# all other solvents in the dataset
lit_sm = {}
for ind, lit_cid in enumerate(lit_cids):
    print(f"Calculating similarity measure with {ind}-{lit_cid}")
    start = time.time()
    sm = []
    pattern_mol = molecules[lit_cid]
    for cid in cids:
        if cid != lit_cids[0]:
            try:
                target_mol = molecules[cid]
                score = similarity_measure(pattern_mol, target_mol, fptype="FP2")
                if score <= 0.1:
                    print(f"low score with {cid}")
                sm.append(score)
            except Exception as e:
                print(e)
    end = time.time()
    print(f"Similarity measure calculation took: {end - start}")
    lit_sm[lit_cid] = sm

# visualize distribution of calculated similarity measures
fig, ax = create_fig(figsize=(12.8, 8), set_axis=False)
n_groups = len(lit_sm.keys())
ncols = 2
nrows = int(np.ceil(n_groups / ncols))
gs = grid_spec.GridSpec(nrows, ncols)
fig.gca().axes.yaxis.set_visible(False)
ax_objs = []
for ind, (key, abb) in enumerate(zip(lit_sm, ABBREVIATIONS)):
    d = fig.add_subplot(gs[int(ind % nrows), int(np.floor(ind / nrows))])
    dist_plot(d, lit_sm[key], bins=60, label=key, color="#777DF2", alpha=0.6)
    d.grid(axis="x", linestyle="--", color="gray")
    d.set_xlim(-0.1, 1.1)
    _set_axis(d, "x")
    d.xaxis.set_major_formatter(FormatStrFormatter("%0.1f"))
    d.set_yticklabels([])
    d.set_ylabel(abb, fontsize=18, weight="bold")
    d.patch.set_alpha(0)
    if ind in [6, 13]:
        d.tick_params(bottom=True, left=False)
        d.set_xlabel("Similarity Score", fontsize=18, weight="bold")
    else:
        d.set(xticklabels=[])
        d.tick_params(bottom=False, left=False)
    edit_plot(d)
ax.axis("off")
gs.update(hspace=0.1)
ax.set_xlabel("Similarity score", fontsize=18)
save_fig(fig, working_dir="../files/results", name="lit_sm.pdf")

# calculate pairwise similarity measure in the solvents dataset, results in
# 483,636 calculated scores
solvent_pairs = list(combinations(cids, 2))
all_sm = []
start = time.time()
for ind, pair in enumerate(solvent_pairs):
    print(f"Calculating similarity measure for pair {ind}")
    try:
        all_sm.append(
            similarity_measure(molecules[pair[0]], molecules[pair[1]], fptype="FP2")
        )
    except Exception as e:
        print(e)
end = time.time()
print(f"Similarity measure calculation took: {end - start}")

# visualize CDF of calculated similarity measures
fig, ax = create_fig(figsize=(6, 3.75), xlabel="Similarity Score", ylabel="CDF")
ecdf_plot(ax, all_sm, color="#777DF2")
ax.vlines(0.85, 0, 1, linestyle="--", color="black", linewidth=2)
y_text = 1.01
for ind, threshold in enumerate([0.25, 0.5, 0.75]):
    y_text -= 0.12
    ax.annotate(
        f"Pr (S $\leq$ {threshold}) = "
        f"{round(len(np.where(np.array(all_sm) <= threshold)[0]) / len(all_sm), 2)}",
        xy=(0.08, y_text),
        fontsize=16,
    )
rect = patches.Rectangle(
    (0.08, 0.65), 0.52, 0.75, facecolor="#C7CAFA", alpha=0.4, transform=fig.transFigure
)
ax.xaxis.set_major_formatter(FormatStrFormatter("%0.1f"))
ax.yaxis.set_major_formatter(FormatStrFormatter("%0.1f"))
ax.add_patch(rect)
edit_plot(ax)
save_fig(fig, working_dir="../files/results", name="all_sm.pdf")
