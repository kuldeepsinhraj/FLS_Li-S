import numpy as np
import pandas as pd
import seaborn as sns

# cluster size
from matplotlib.colors import LinearSegmentedColormap
from scipy.interpolate import splev, splrep

from ebatt.visualization.structural import save_fig, dist_plot, edit_plot, create_fig

# palette = sns.color_palette("Set1")
colors = ["#424953", "#D94452", "#5E9CEA", "#8AC054", "#FCD277", "#967ADA"]

project_dir = "/Users/rashaatwi/Desktop/lis/ebatt/md"

results_dir = f"{project_dir}/Results"
analysis_dir = f"{project_dir}/Analysis"

sys1_dir = f"{project_dir}/1-DOL-DME/Analysis"
sys2_dir = f"{project_dir}/2-DOL-DME-PS/Analysis"
sys3_dir = f"{project_dir}/3-DOL-37027-PS/Analysis"
sys4_dir = f"{project_dir}/4-DOL-18335957-PS/Analysis"
sys5_dir = f"{project_dir}/5-DOL-92551-PS/Analysis"
sys6_dir = f"{project_dir}/6-DOL-2776662-PS/Analysis"
sys7_dir = f"{project_dir}/7-DOL-90838605-PS/Analysis"
sys8_dir = f"{project_dir}/8-DOL-139775149-PS/Analysis"
sys9_dir = f"{project_dir}/9-DOL-19992615-PS/Analysis"
sys10_dir = f"{project_dir}/10-DOL-88850722-PS/Analysis"
sys11_dir = f"{project_dir}/11-DOL-91607507-PS/Analysis"
sys12_dir = f"{project_dir}/12-DOL-149597093-PS/Analysis"
sys13_dir = f"{project_dir}/13-DOL-89114499-PS/Analysis"

rdf_sys1 = pd.read_csv(f"{sys1_dir}/RDF/rdf.csv")
rdf_sys2 = pd.read_csv(f"{sys2_dir}/RDF/rdf.csv")
rdf_sys3 = pd.read_csv(f"{sys3_dir}/RDF/rdf.csv")
rdf_sys4 = pd.read_csv(f"{sys4_dir}/RDF/rdf.csv")
rdf_sys5 = pd.read_csv(f"{sys5_dir}/RDF/rdf.csv")
rdf_sys6 = pd.read_csv(f"{sys6_dir}/RDF/rdf.csv")
rdf_sys7 = pd.read_csv(f"{sys7_dir}/RDF/rdf.csv")
rdf_sys8 = pd.read_csv(f"{sys8_dir}/RDF/rdf.csv")
rdf_sys9 = pd.read_csv(f"{sys9_dir}/RDF/rdf.csv")
rdf_sys10 = pd.read_csv(f"{sys10_dir}/RDF/rdf.csv")
rdf_sys11 = pd.read_csv(f"{sys11_dir}/RDF/rdf.csv")
rdf_sys12 = pd.read_csv(f"{sys12_dir}/RDF/rdf.csv")
rdf_sys13 = pd.read_csv(f"{sys13_dir}/RDF/rdf.csv")

sys1_atoms = {
    "O (DOL)": [10, 11],
    "O (DME)": [12, 16],
    "N (TFSI)": [28],
    "O (TFSI)": [30, 31, 37, 38],
    "Li": [43],
}

sys2_atoms = {
    "O (DOL)": [10, 11],
    "O (DME)": [12, 16],
    "O (NO3)": [28, 29, 30],
    "N (TFSI)": [32],
    "O (TFSI)": [34, 35, 41, 42],
    "St (S4)": [47, 50],
    "Si (S4)": [48, 49],
    "Li": [51],
}
sys3_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15, 16, 17, 18],
    "O (FLS)": [19],
    "O (NO3)": [27, 28, 29],
    "N (TFSI)": [31],
    "O (TFSI)": [33, 34, 40, 41],
    "St (S4)": [46, 49],
    "Si (S4)": [47, 48],
    "Li": [50],
}

sys4_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12],
    "O (FLS)": [13],
    "O (NO3)": [27, 28, 29],
    "N (TFSI)": [31],
    "O (TFSI)": [33, 34, 40, 41],
    "St (S4)": [46, 49],
    "Si (S4)": [47, 48],
    "Li": [50],
}

sys5_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15, 16, 17, 18],
    "O (FLS)": [19],
    "O (NO3)": [27, 28, 29],
    "N (TFSI)": [31],
    "O (TFSI)": [33, 34, 40, 41],
    "St (S4)": [46, 49],
    "Si (S4)": [47, 48],
    "Li": [50],
}

sys6_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15, 16, 17, 18, 19],
    "O (FLS)": [20],
    "O (NO3)": [30, 31, 32],
    "N (TFSI)": [34],
    "O (TFSI)": [36, 37, 44, 45],
    "St (S4)": [49, 52],
    "Si (S4)": [50, 51],
    "Li": [53],
}

sys7_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15],
    "O (FLS)": [16, 17],
    "O (NO3)": [28, 29, 30],
    "N (TFSI)": [32],
    "O (TFSI)": [34, 35, 41, 42],
    "St (S4)": [47, 50],
    "Si (S4)": [48, 49],
    "Li": [51],
}

sys8_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15, 16],
    "O (FLS)": [17],
    "O (NO3)": [27, 28, 29],
    "N (TFSI)": [31],
    "O (TFSI)": [33, 34, 40, 41],
    "St (S4)": [46, 49],
    "Si (S4)": [47, 48],
    "Li": [50],
}

sys9_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15, 16, 17],
    "O (FLS)": [18],
    "O (NO3)": [27, 28, 29],
    "N (TFSI)": [31],
    "O (TFSI)": [33, 34, 40, 41],
    "St (S4)": [46, 49],
    "Si (S4)": [47, 48],
    "Li": [50],
}

sys10_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15, 16],
    "O (FLS)": [17],
    "O (NO3)": [30, 31, 32],
    "N (TFSI)": [34],
    "O (TFSI)": [36, 37, 43, 44],
    "St (S4)": [49, 52],
    "Si (S4)": [50, 51],
    "Li": [53],
}

sys11_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15],
    "O (FLS)": [15],
    "O (NO3)": [30, 31, 32],
    "N (TFSI)": [34],
    "O (TFSI)": [36, 37, 43, 44],
    "St (S4)": [49, 52],
    "Si (S4)": [50, 51],
    "Li": [53],
}

sys12_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15, 16],
    "O (FLS)": [17],
    "O (NO3)": [30, 31, 32],
    "N (TFSI)": [34],
    "O (TFSI)": [36, 37, 43, 44],
    "St (S4)": [49, 52],
    "Si (S4)": [50, 51],
    "Li": [53],
}

sys13_atoms = {
    "O (DOL)": [10, 11],
    "F (FLS)": [12, 13, 14, 15, 16],
    "O (FLS)": [17],
    "O (NO3)": [30, 31, 32],
    "N (TFSI)": [34],
    "O (TFSI)": [36, 37, 43, 44],
    "St (S4)": [49, 52],
    "Si (S4)": [50, 51],
    "Li": [53],
}

sys1_li = sys1_atoms["Li"]
sys2_li = sys2_atoms["Li"]
sys3_li = sys3_atoms["Li"]
sys4_li = sys4_atoms["Li"]
sys5_li = sys5_atoms["Li"]
sys6_li = sys6_atoms["Li"]
sys7_li = sys7_atoms["Li"]
sys8_li = sys8_atoms["Li"]
sys9_li = sys9_atoms["Li"]
sys10_li = sys10_atoms["Li"]
sys11_li = sys11_atoms["Li"]
sys12_li = sys12_atoms["Li"]
sys13_li = sys13_atoms["Li"]

# Li-Solvent RDF
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.plot(
    rdf_sys1["r ($\AA$)"].values,
    sum([rdf_sys1[f"g_{sys1_li[0]}-{i}"].values for i in sys1_atoms["O (DME)"]])
    / len(sys1_atoms["O (DME)"]),
    label="DOL/DME",
    color=colors[0],
    linewidth=2,
)
ax.plot(
    rdf_sys2["r ($\AA$)"].values,
    sum([rdf_sys2[f"g_{sys2_li[0]}-{i}"].values for i in sys2_atoms["O (DME)"]])
    / len(sys2_atoms["O (DME)"]),
    label="DOL/DME (PS)",
    color=colors[1],
    linewidth=2,
)
ax.plot(
    rdf_sys3["r ($\AA$)"].values,
    sum([rdf_sys3[f"g_{sys3_li[0]}-{i}"].values for i in sys3_atoms["O (FLS)"]])
    / len(sys3_atoms["O (FLS)"]),
    label="DOL/FLS",
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys4["r ($\AA$)"].values,
    sum([rdf_sys4[f"g_{sys4_li[0]}-{i}"].values for i in sys4_atoms["O (FLS)"]])
    / len(sys4_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys5["r ($\AA$)"].values,
    sum([rdf_sys5[f"g_{sys5_li[0]}-{i}"].values for i in sys5_atoms["O (FLS)"]])
    / len(sys5_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys6["r ($\AA$)"].values,
    sum([rdf_sys6[f"g_{sys6_li[0]}-{i}"].values for i in sys6_atoms["O (FLS)"]])
    / len(sys6_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys7["r ($\AA$)"].values,
    sum([rdf_sys7[f"g_{sys7_li[0]}-{i}"].values for i in sys7_atoms["O (FLS)"]])
    / len(sys7_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys8["r ($\AA$)"].values,
    sum([rdf_sys8[f"g_{sys8_li[0]}-{i}"].values for i in sys8_atoms["O (FLS)"]])
    / len(sys8_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys8["r ($\AA$)"].values,
    sum([rdf_sys8[f"g_{sys8_li[0]}-{i}"].values for i in sys8_atoms["O (FLS)"]])
    / len(sys8_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys9["r ($\AA$)"].values,
    sum([rdf_sys9[f"g_{sys9_li[0]}-{i}"].values for i in sys9_atoms["O (FLS)"]])
    / len(sys9_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys10["r ($\AA$)"].values,
    sum([rdf_sys10[f"g_{sys10_li[0]}-{i}"].values for i in sys10_atoms["O (FLS)"]])
    / len(sys10_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys11["r ($\AA$)"].values,
    sum([rdf_sys11[f"g_{sys11_li[0]}-{i}"].values for i in sys11_atoms["O (FLS)"]])
    / len(sys11_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys12["r ($\AA$)"].values,
    sum([rdf_sys12[f"g_{sys12_li[0]}-{i}"].values for i in sys12_atoms["O (FLS)"]])
    / len(sys12_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys13["r ($\AA$)"].values,
    sum([rdf_sys13[f"g_{sys13_li[0]}-{i}"].values for i in sys13_atoms["O (FLS)"]])
    / len(sys13_atoms["O (FLS)"]),
    color=colors[-1],
    linewidth=2,
)
ax.set_xlim(0, 5)
ax.set_xlabel("r ($\AA$)", fontsize=18)
ax.set_ylabel("g (r)", fontsize=18)
ax.legend(
    loc="upper right", frameon=False, fontsize=18,
)
save_fig(fig, working_dir=results_dir, name="rdf_solvents.png")

# Li-Anion RDF
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.plot(
    rdf_sys1["r ($\AA$)"].values,
    sum([rdf_sys1[f"g_{sys1_li[0]}-{i}"].values for i in sys1_atoms["O (TFSI)"]])
    / len(sys1_atoms["O (TFSI)"]),
    label="DOL/DME",
    color=colors[0],
    linewidth=2,
)
ax.plot(
    rdf_sys2["r ($\AA$)"].values,
    sum([rdf_sys2[f"g_{sys2_li[0]}-{i}"].values for i in sys2_atoms["O (TFSI)"]])
    / len(sys2_atoms["O (TFSI)"]),
    label="DOL/DME (PS)",
    color=colors[1],
    linewidth=2,
)
ax.plot(
    rdf_sys3["r ($\AA$)"].values,
    sum([rdf_sys3[f"g_{sys3_li[0]}-{i}"].values for i in sys3_atoms["O (TFSI)"]])
    / len(sys2_atoms["O (TFSI)"]),
    label="DOL/FLS",
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys4["r ($\AA$)"].values,
    sum([rdf_sys4[f"g_{sys4_li[0]}-{i}"].values for i in sys4_atoms["O (TFSI)"]])
    / len(sys4_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys5["r ($\AA$)"].values,
    sum([rdf_sys5[f"g_{sys5_li[0]}-{i}"].values for i in sys5_atoms["O (TFSI)"]])
    / len(sys5_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys6["r ($\AA$)"].values,
    sum([rdf_sys6[f"g_{sys6_li[0]}-{i}"].values for i in sys6_atoms["O (TFSI)"]])
    / len(sys6_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys7["r ($\AA$)"].values,
    sum([rdf_sys7[f"g_{sys7_li[0]}-{i}"].values for i in sys7_atoms["O (TFSI)"]])
    / len(sys7_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys8["r ($\AA$)"].values,
    sum([rdf_sys8[f"g_{sys8_li[0]}-{i}"].values for i in sys8_atoms["O (TFSI)"]])
    / len(sys8_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys9["r ($\AA$)"].values,
    sum([rdf_sys9[f"g_{sys9_li[0]}-{i}"].values for i in sys9_atoms["O (TFSI)"]])
    / len(sys9_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys10["r ($\AA$)"].values,
    sum([rdf_sys10[f"g_{sys10_li[0]}-{i}"].values for i in sys10_atoms["O (TFSI)"]])
    / len(sys10_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys11["r ($\AA$)"].values,
    sum([rdf_sys11[f"g_{sys11_li[0]}-{i}"].values for i in sys11_atoms["O (TFSI)"]])
    / len(sys11_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys12["r ($\AA$)"].values,
    sum([rdf_sys12[f"g_{sys12_li[0]}-{i}"].values for i in sys12_atoms["O (TFSI)"]])
    / len(sys12_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys13["r ($\AA$)"].values,
    sum([rdf_sys13[f"g_{sys13_li[0]}-{i}"].values for i in sys13_atoms["O (TFSI)"]])
    / len(sys13_atoms["O (TFSI)"]),
    color=colors[-1],
    linewidth=2,
)
ax.set_xlim(0, 8)
ax.set_yticks(np.arange(0, 40, 10.0))
ax.set_xlabel("r ($\AA$)", fontsize=18)
ax.set_ylabel("g (r)", fontsize=18)
ax.legend(loc="upper right", frameon=False, fontsize=18)
save_fig(fig, working_dir=results_dir, name="rdf_anion.png")

# Li-Additive RDF
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.plot(
    rdf_sys2["r ($\AA$)"].values,
    sum([rdf_sys2[f"g_{sys2_li[0]}-{i}"].values for i in sys2_atoms["O (NO3)"]])
    / len(sys2_atoms["O (NO3)"]),
    label="DOL/DME (PS)",
    color=colors[0],
    linewidth=2,
)
ax.plot(
    rdf_sys3["r ($\AA$)"].values,
    sum([rdf_sys3[f"g_{sys3_li[0]}-{i}"].values for i in sys3_atoms["O (NO3)"]])
    / len(sys3_atoms["O (NO3)"]),
    label="DOL/FLS",
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys4["r ($\AA$)"].values,
    sum([rdf_sys4[f"g_{sys4_li[0]}-{i}"].values for i in sys4_atoms["O (NO3)"]])
    / len(sys4_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys5["r ($\AA$)"].values,
    sum([rdf_sys5[f"g_{sys5_li[0]}-{i}"].values for i in sys5_atoms["O (NO3)"]])
    / len(sys5_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys6["r ($\AA$)"].values,
    sum([rdf_sys6[f"g_{sys6_li[0]}-{i}"].values for i in sys6_atoms["O (NO3)"]])
    / len(sys6_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys7["r ($\AA$)"].values,
    sum([rdf_sys7[f"g_{sys7_li[0]}-{i}"].values for i in sys7_atoms["O (NO3)"]])
    / len(sys7_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys8["r ($\AA$)"].values,
    sum([rdf_sys8[f"g_{sys8_li[0]}-{i}"].values for i in sys8_atoms["O (NO3)"]])
    / len(sys8_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys9["r ($\AA$)"].values,
    sum([rdf_sys9[f"g_{sys9_li[0]}-{i}"].values for i in sys9_atoms["O (NO3)"]])
    / len(sys9_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys10["r ($\AA$)"].values,
    sum([rdf_sys10[f"g_{sys10_li[0]}-{i}"].values for i in sys10_atoms["O (NO3)"]])
    / len(sys10_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys11["r ($\AA$)"].values,
    sum([rdf_sys11[f"g_{sys11_li[0]}-{i}"].values for i in sys11_atoms["O (NO3)"]])
    / len(sys11_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys12["r ($\AA$)"].values,
    sum([rdf_sys12[f"g_{sys12_li[0]}-{i}"].values for i in sys12_atoms["O (NO3)"]])
    / len(sys12_atoms["O (NO3)"]),
    color=colors[-1],
    linewidth=2,
)

ax.set_xlim(0, 8)
# ax.set_yticks(np.arange(0, 50, 10.0))
ax.set_xlabel("r ($\AA$)", fontsize=18)
ax.set_ylabel("g (r)", fontsize=18)
ax.legend(loc="upper right", frameon=False, fontsize=18)
save_fig(fig, working_dir=results_dir, name="rdf_additive.png")

# Li-PS RDF
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.plot(
    rdf_sys2["r ($\AA$)"].values,
    sum([rdf_sys2[f"g_{sys2_li[0]}-{i}"].values for i in sys2_atoms["St (S4)"]])
    / len(sys2_atoms["St (S4)"]),
    label="DOL/DME (PS)",
    color=colors[1],
    linewidth=2,
)
ax.plot(
    rdf_sys3["r ($\AA$)"].values,
    sum([rdf_sys3[f"g_{sys3_li[0]}-{i}"].values for i in sys3_atoms["St (S4)"]])
    / len(sys3_atoms["St (S4)"]),
    label="DOL/FLS",
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys4["r ($\AA$)"].values,
    sum([rdf_sys4[f"g_{sys3_li[0]}-{i}"].values for i in sys4_atoms["St (S4)"]])
    / len(sys4_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys5["r ($\AA$)"].values,
    sum([rdf_sys5[f"g_{sys5_li[0]}-{i}"].values for i in sys5_atoms["St (S4)"]])
    / len(sys5_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys6["r ($\AA$)"].values,
    sum([rdf_sys6[f"g_{sys6_li[0]}-{i}"].values for i in sys6_atoms["St (S4)"]])
    / len(sys6_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys7["r ($\AA$)"].values,
    sum([rdf_sys7[f"g_{sys7_li[0]}-{i}"].values for i in sys7_atoms["St (S4)"]])
    / len(sys7_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys8["r ($\AA$)"].values,
    sum([rdf_sys8[f"g_{sys8_li[0]}-{i}"].values for i in sys8_atoms["St (S4)"]])
    / len(sys8_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys9["r ($\AA$)"].values,
    sum([rdf_sys9[f"g_{sys9_li[0]}-{i}"].values for i in sys9_atoms["St (S4)"]])
    / len(sys9_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys10["r ($\AA$)"].values,
    sum([rdf_sys10[f"g_{sys10_li[0]}-{i}"].values for i in sys10_atoms["St (S4)"]])
    / len(sys10_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys11["r ($\AA$)"].values,
    sum([rdf_sys11[f"g_{sys11_li[0]}-{i}"].values for i in sys11_atoms["St (S4)"]])
    / len(sys11_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys12["r ($\AA$)"].values,
    sum([rdf_sys12[f"g_{sys12_li[0]}-{i}"].values for i in sys12_atoms["St (S4)"]])
    / len(sys12_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys13["r ($\AA$)"].values,
    sum([rdf_sys13[f"g_{sys13_li[0]}-{i}"].values for i in sys13_atoms["St (S4)"]])
    / len(sys13_atoms["St (S4)"]),
    color=colors[-1],
    linewidth=2,
)

ax.plot(
    rdf_sys2["r ($\AA$)"].values,
    sum([rdf_sys2[f"g_{sys2_li[0]}-{i}"].values for i in sys2_atoms["Si (S4)"]])
    / len(sys2_atoms["Si (S4)"]),
    color=colors[1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys3["r ($\AA$)"].values,
    sum([rdf_sys3[f"g_{sys3_li[0]}-{i}"].values for i in sys3_atoms["Si (S4)"]])
    / len(sys3_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys4["r ($\AA$)"].values,
    sum([rdf_sys4[f"g_{sys3_li[0]}-{i}"].values for i in sys4_atoms["Si (S4)"]])
    / len(sys4_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys5["r ($\AA$)"].values,
    sum([rdf_sys5[f"g_{sys5_li[0]}-{i}"].values for i in sys5_atoms["Si (S4)"]])
    / len(sys5_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys6["r ($\AA$)"].values,
    sum([rdf_sys6[f"g_{sys6_li[0]}-{i}"].values for i in sys6_atoms["Si (S4)"]])
    / len(sys6_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys7["r ($\AA$)"].values,
    sum([rdf_sys7[f"g_{sys7_li[0]}-{i}"].values for i in sys7_atoms["Si (S4)"]])
    / len(sys7_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys8["r ($\AA$)"].values,
    sum([rdf_sys8[f"g_{sys8_li[0]}-{i}"].values for i in sys8_atoms["Si (S4)"]])
    / len(sys8_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys9["r ($\AA$)"].values,
    sum([rdf_sys9[f"g_{sys9_li[0]}-{i}"].values for i in sys9_atoms["Si (S4)"]])
    / len(sys9_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys10["r ($\AA$)"].values,
    sum([rdf_sys10[f"g_{sys10_li[0]}-{i}"].values for i in sys10_atoms["Si (S4)"]])
    / len(sys10_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys11["r ($\AA$)"].values,
    sum([rdf_sys11[f"g_{sys11_li[0]}-{i}"].values for i in sys11_atoms["Si (S4)"]])
    / len(sys11_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys12["r ($\AA$)"].values,
    sum([rdf_sys12[f"g_{sys12_li[0]}-{i}"].values for i in sys12_atoms["Si (S4)"]])
    / len(sys12_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)
ax.plot(
    rdf_sys13["r ($\AA$)"].values,
    sum([rdf_sys13[f"g_{sys13_li[0]}-{i}"].values for i in sys13_atoms["Si (S4)"]])
    / len(sys13_atoms["Si (S4)"]),
    color=colors[-1],
    linewidth=2,
    linestyle="--",
)

ax.set_xlim(0, 8)
ax.set_yticks(np.arange(0, 140, 30.0))
ax.set_xlabel("r ($\AA$)", fontsize=18)
ax.set_ylabel("g (r)", fontsize=18)
ax.legend(loc="upper right", frameon=False, fontsize=18)
save_fig(fig, working_dir=results_dir, name="rdf_ps.png")

# Li-Cation RDF
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.plot(
    rdf_sys1["r ($\AA$)"].values,
    sum([rdf_sys1[f"g_{sys1_li[0]}-{i}"].values for i in sys1_atoms["Li"]])
    / len(sys1_atoms["Li"]),
    label="DOL/DME",
    color=colors[0],
    linewidth=2,
)
ax.plot(
    rdf_sys2["r ($\AA$)"].values,
    sum([rdf_sys2[f"g_{sys2_li[0]}-{i}"].values for i in sys2_atoms["Li"]])
    / len(sys2_atoms["Li"]),
    label="DOL/DME (PS)",
    color=colors[1],
    linewidth=2,
)
ax.plot(
    rdf_sys3["r ($\AA$)"].values,
    sum([rdf_sys3[f"g_{sys3_li[0]}-{i}"].values for i in sys3_atoms["Li"]])
    / len(sys3_atoms["Li"]),
    label="DOL/FLS",
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys4["r ($\AA$)"].values,
    sum([rdf_sys4[f"g_{sys4_li[0]}-{i}"].values for i in sys4_atoms["Li"]])
    / len(sys4_atoms["Li"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys5["r ($\AA$)"].values,
    sum([rdf_sys5[f"g_{sys5_li[0]}-{i}"].values for i in sys5_atoms["Li"]])
    / len(sys5_atoms["Li"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys6["r ($\AA$)"].values,
    sum([rdf_sys6[f"g_{sys6_li[0]}-{i}"].values for i in sys6_atoms["Li"]])
    / len(sys6_atoms["Li"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys7["r ($\AA$)"].values,
    sum([rdf_sys7[f"g_{sys7_li[0]}-{i}"].values for i in sys7_atoms["Li"]])
    / len(sys7_atoms["Li"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys8["r ($\AA$)"].values,
    sum([rdf_sys8[f"g_{sys8_li[0]}-{i}"].values for i in sys8_atoms["Li"]])
    / len(sys8_atoms["Li"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys9["r ($\AA$)"].values,
    sum([rdf_sys9[f"g_{sys9_li[0]}-{i}"].values for i in sys9_atoms["Li"]])
    / len(sys9_atoms["Li"]),
    linewidth=2,
    color=colors[-1],
)
ax.plot(
    rdf_sys10["r ($\AA$)"].values,
    sum([rdf_sys10[f"g_{sys10_li[0]}-{i}"].values for i in sys10_atoms["Li"]])
    / len(sys10_atoms["Li"]),
    linewidth=2,
    color=colors[-1],
)
ax.plot(
    rdf_sys11["r ($\AA$)"].values,
    sum([rdf_sys11[f"g_{sys11_li[0]}-{i}"].values for i in sys11_atoms["Li"]])
    / len(sys11_atoms["Li"]),
    linewidth=2,
    color=colors[-1],
)
ax.plot(
    rdf_sys12["r ($\AA$)"].values,
    sum([rdf_sys12[f"g_{sys12_li[0]}-{i}"].values for i in sys12_atoms["Li"]])
    / len(sys12_atoms["Li"]),
    linewidth=2,
    color=colors[-1],
)
ax.plot(
    rdf_sys13["r ($\AA$)"].values,
    sum([rdf_sys13[f"g_{sys13_li[0]}-{i}"].values for i in sys13_atoms["Li"]])
    / len(sys13_atoms["Li"]),
    linewidth=2,
    color=colors[-1],
)

# ax.set_xlim(0, 8)
ax.set_xlabel("r ($\AA$)", fontsize=18)
ax.set_ylabel("g (r)", fontsize=18)
ax.legend(loc="upper right", frameon=False, fontsize=18)
save_fig(fig, working_dir=results_dir, name="rdf_cation.png")

# Li-DOL RDF
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.plot(
    rdf_sys1["r ($\AA$)"].values,
    sum([rdf_sys1[f"g_{sys1_li[0]}-{i}"].values for i in sys1_atoms["O (DOL)"]])
    / len(sys1_atoms["O (DOL)"]),
    label="DOL/DME",
    color=colors[0],
    linewidth=2,
)
ax.plot(
    rdf_sys2["r ($\AA$)"].values,
    sum([rdf_sys2[f"g_{sys2_li[0]}-{i}"].values for i in sys2_atoms["O (DOL)"]])
    / len(sys2_atoms["O (DOL)"]),
    label="DOL/DME (PS)",
    color=colors[1],
    linewidth=2,
)
ax.plot(
    rdf_sys3["r ($\AA$)"].values,
    sum([rdf_sys3[f"g_{sys3_li[0]}-{i}"].values for i in sys3_atoms["O (DOL)"]])
    / len(sys3_atoms["O (DOL)"]),
    label="DOL/FLS",
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys4["r ($\AA$)"].values,
    sum([rdf_sys4[f"g_{sys4_li[0]}-{i}"].values for i in sys4_atoms["O (DOL)"]])
    / len(sys4_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys5["r ($\AA$)"].values,
    sum([rdf_sys5[f"g_{sys5_li[0]}-{i}"].values for i in sys5_atoms["O (DOL)"]])
    / len(sys5_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys6["r ($\AA$)"].values,
    sum([rdf_sys6[f"g_{sys6_li[0]}-{i}"].values for i in sys6_atoms["O (DOL)"]])
    / len(sys6_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys7["r ($\AA$)"].values,
    sum([rdf_sys7[f"g_{sys7_li[0]}-{i}"].values for i in sys7_atoms["O (DOL)"]])
    / len(sys7_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys8["r ($\AA$)"].values,
    sum([rdf_sys8[f"g_{sys8_li[0]}-{i}"].values for i in sys8_atoms["O (DOL)"]])
    / len(sys8_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys9["r ($\AA$)"].values,
    sum([rdf_sys9[f"g_{sys9_li[0]}-{i}"].values for i in sys9_atoms["O (DOL)"]])
    / len(sys9_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys10["r ($\AA$)"].values,
    sum([rdf_sys10[f"g_{sys10_li[0]}-{i}"].values for i in sys10_atoms["O (DOL)"]])
    / len(sys10_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys11["r ($\AA$)"].values,
    sum([rdf_sys11[f"g_{sys11_li[0]}-{i}"].values for i in sys11_atoms["O (DOL)"]])
    / len(sys11_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys12["r ($\AA$)"].values,
    sum([rdf_sys12[f"g_{sys12_li[0]}-{i}"].values for i in sys12_atoms["O (DOL)"]])
    / len(sys12_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)
ax.plot(
    rdf_sys13["r ($\AA$)"].values,
    sum([rdf_sys13[f"g_{sys13_li[0]}-{i}"].values for i in sys13_atoms["O (DOL)"]])
    / len(sys13_atoms["O (DOL)"]),
    color=colors[-1],
    linewidth=2,
)

ax.set_xlim(0, 8)
ax.set_xlabel("r ($\AA$)", fontsize=18)
ax.set_ylabel("g (r)", fontsize=18)
ax.legend(loc="upper right", frameon=False, fontsize=18)
save_fig(fig, working_dir=results_dir, name="rdf_dol.png")


# PS-PS RDF


def get_ps_rdf(file):
    s_rdf = pd.read_csv(file)
    x = np.array(s_rdf["r ($\AA$)"])
    x_smooth = np.linspace(x.min(), x.max(), 100)
    rdf = pd.DataFrame(x_smooth, columns=["r ($\AA$)"])
    y = np.array(s_rdf["g_5-5"])
    spl = splrep(x, y)
    y_smooth = splev(x_smooth, spl)
    rdf["g_5-5"] = y_smooth
    rdf = rdf.melt("r ($\AA$)", var_name="cols", value_name="g(r)")
    rdf = rdf[rdf["g(r)"] >= 0]
    return rdf


rdf_ps_sys2 = get_ps_rdf(f"{sys2_dir}/RDF/rdf_mol.csv")
rdf_ps_sys3 = get_ps_rdf(f"{sys3_dir}/RDF/rdf_mol.csv")
rdf_ps_sys4 = get_ps_rdf(f"{sys4_dir}/RDF/rdf_mol.csv")
rdf_ps_sys5 = get_ps_rdf(f"{sys5_dir}/RDF/rdf_mol.csv")
rdf_ps_sys6 = get_ps_rdf(f"{sys6_dir}/RDF/rdf_mol.csv")
rdf_ps_sys7 = get_ps_rdf(f"{sys7_dir}/RDF/rdf_mol.csv")
rdf_ps_sys8 = get_ps_rdf(f"{sys8_dir}/RDF/rdf_mol.csv")
rdf_ps_sys9 = get_ps_rdf(f"{sys9_dir}/RDF/rdf_mol.csv")
rdf_ps_sys10 = get_ps_rdf(f"{sys10_dir}/RDF/rdf_mol.csv")
rdf_ps_sys11 = get_ps_rdf(f"{sys11_dir}/RDF/rdf_mol.csv")
rdf_ps_sys12 = get_ps_rdf(f"{sys12_dir}/RDF/rdf_mol.csv")
rdf_ps_sys13 = get_ps_rdf(f"{sys13_dir}/RDF/rdf_mol.csv")

fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.plot(
    rdf_ps_sys2["r ($\AA$)"].values,
    rdf_ps_sys2["g(r)"].values,
    label="DOL/DME (PS)",
    color=colors[1],
    linewidth=2,
)
ax.plot(
    rdf_ps_sys4["r ($\AA$)"].values,
    rdf_ps_sys4["g(r)"].values,
    label="DOL/FLS",
    color=colors[-1],
    linewidth=2,
)

# for rdf in [rdf_ps_sys4,
#             rdf_ps_sys5]:
#     ax.plot(
#         rdf["r ($\AA$)"].values,
#         rdf["g(r)"].values,
#         label="DOL/FLS",
#         color=colors[-1],
#         linewidth=2,
#     )

# ax.set_xlim(0, 8)
ax.set_xticks(np.arange(0, 20, 3.0))
ax.set_xlabel("r ($\AA$)", fontsize=18)
ax.set_ylabel("g (r)", fontsize=18)
ax.legend(loc="upper right", frameon=False, fontsize=18)
save_fig(fig, working_dir=results_dir, name="rdf_ps_1.png")

# Bridging effect
sys2_clusters = pd.read_csv(f"{sys2_dir}/Clusters/Li/df.csv")
sys3_clusters = pd.read_csv(f"{sys3_dir}/Clusters/Li/df.csv")
sys4_clusters = pd.read_csv(f"{sys4_dir}/Clusters/Li/df.csv")
sys5_clusters = pd.read_csv(f"{sys5_dir}/Clusters/Li/df.csv")
sys6_clusters = pd.read_csv(f"{sys6_dir}/Clusters/Li/df.csv")
sys7_clusters = pd.read_csv(f"{sys7_dir}/Clusters/Li/df.csv")
sys8_clusters = pd.read_csv(f"{sys8_dir}/Clusters/Li/df.csv")
sys9_clusters = pd.read_csv(f"{sys9_dir}/Clusters/Li/df.csv")
sys10_clusters = pd.read_csv(f"{sys10_dir}/Clusters/Li/df.csv")
sys11_clusters = pd.read_csv(f"{sys11_dir}/Clusters/Li/df.csv")
sys12_clusters = pd.read_csv(f"{sys12_dir}/Clusters/Li/df.csv")
sys13_clusters = pd.read_csv(f"{sys13_dir}/Clusters/Li/df.csv")


x = [0, 70, 10, 70, 66.67, 40, 50, 60, 41.67, 33.33, 41.67, 41.67]
y = []
for i in [
    sys2_clusters,
    sys3_clusters,
    sys4_clusters,
    sys5_clusters,
    sys6_clusters,
    sys7_clusters,
    sys8_clusters,
    sys9_clusters,
    sys10_clusters,
    sys11_clusters,
    sys12_clusters,
    sys13_clusters,
]:
    y.append(sum(i[i["S4"] > 1]["count"].values) * 100 / sum(i["count"].values))

fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.bar(x, y, width=2.5, color="#DDDDDD", edgecolor="k", linewidth=2)
# ax.set_xticks(np.arange(0, 90, 20))
# ax.set_yticks(np.arange(0, 10, 2.0))
ax.set_xlabel("Fluorination degree (%)", fontsize=24)
ax.set_ylabel("Clustering (%)", fontsize=24)
save_fig(fig, working_dir=results_dir, name="bridging_effect.png")


# clusters
def custom_cmap(
    numcolors=100,
    name="custom_cmap",
    mincol="#F9FF95",
    midcol1="#C0C0C0",
    midcol2="#CEDEF2",
    maxcol="#7462A3",
):
    cmap = LinearSegmentedColormap.from_list(
        name=name, colors=[mincol, midcol1, midcol2, maxcol], N=numcolors
    )
    return cmap


sys2_dist = pd.read_csv(f"{sys2_dir}/Clusters/PS/df.csv")
sys3_dist = pd.read_csv(f"{sys3_dir}/Clusters/PS/df.csv")
sys4_dist = pd.read_csv(f"{sys4_dir}/Clusters/PS/df.csv")
sys5_dist = pd.read_csv(f"{sys5_dir}/Clusters/PS/df.csv")
sys6_dist = pd.read_csv(f"{sys6_dir}/Clusters/PS/df.csv")
sys7_dist = pd.read_csv(f"{sys7_dir}/Clusters/PS/df.csv")
sys8_dist = pd.read_csv(f"{sys8_dir}/Clusters/PS/df.csv")
sys9_dist = pd.read_csv(f"{sys9_dir}/Clusters/PS/df.csv")
sys10_dist = pd.read_csv(f"{sys10_dir}/Clusters/PS/df.csv")
sys11_dist = pd.read_csv(f"{sys11_dir}/Clusters/PS/df.csv")
sys12_dist = pd.read_csv(f"{sys12_dir}/Clusters/PS/df.csv")
sys13_dist = pd.read_csv(f"{sys13_dir}/Clusters/PS/df.csv")
size = [i for i in range(1, 8)]
p2 = []
p3 = []
p4 = []
p5 = []
p6 = []
p7 = []
p8 = []
p9 = []
p10 = []
p11 = []
p12 = []
p13 = []

for i in size:
    p2.append(sum(sys2_dist[sys2_dist["S4"] == i]["%"]) / 100)
for i in size:
    p3.append(sum(sys3_dist[sys3_dist["S4"] == i]["%"]) / 100)
for i in size:
    p4.append(sum(sys4_dist[sys4_dist["S4"] == i]["%"]) / 100)
for i in size:
    p5.append(sum(sys5_dist[sys5_dist["S4"] == i]["%"]) / 100)
for i in size:
    p6.append(sum(sys6_dist[sys6_dist["S4"] == i]["%"]) / 100)
for i in size:
    p7.append(sum(sys7_dist[sys7_dist["S4"] == i]["%"]) / 100)
for i in size:
    p8.append(sum(sys8_dist[sys8_dist["S4"] == i]["%"]) / 100)
for i in size:
    p9.append(sum(sys9_dist[sys9_dist["S4"] == i]["%"]) / 100)
for i in size:
    p10.append(sum(sys10_dist[sys10_dist["S4"] == i]["%"]) / 100)
for i in size:
    p11.append(sum(sys11_dist[sys11_dist["S4"] == i]["%"]) / 100)
for i in size:
    p12.append(sum(sys12_dist[sys12_dist["S4"] == i]["%"]) / 100)
for i in size:
    p13.append(sum(sys13_dist[sys13_dist["S4"] == i]["%"]) / 100)

cmap = custom_cmap()
new_df = pd.DataFrame(
    columns=size, data=[p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13]
)
new_df["System"] = [
    "DOL/DME",
    "DOL/FLS1",
    "DOL/FLS2",
    "DOL/FLS3",
    "DOL/FLS4",
    "DOL/FLS5",
    "DOL/FLS6",
    "DOL/FLS7",
    "DOL/FLS8",
    "DOL/FLS9",
    "DOL/FLS10",
    "DOL/FLS11",
]
new_df.set_index("System", inplace=True)

fig, ax = create_fig(figsize=(12, 12), set_axis=True, axis="both")
p = sns.heatmap(
    new_df,
    cmap=cmap,
    annot=True,
    fmt=".2f",
    annot_kws={"size": 18},
    cbar_kws={"label": "Probability (size)"},
    linewidths=0.2,
    square=True,
    ax=ax,
)
ax.figure.axes[-1].yaxis.label.set_size(20)
ax.figure.axes[-1].yaxis.set_tick_params(labelsize=20)
ax.figure.axes[-1].yaxis.get_major_ticks()[0].set_visible(False)
ax.tick_params(axis="x", which="both", bottom=False, top=False, length=0)
ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20)
ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20)
for text in ax.texts:
    if text.get_text() == "0.00":
        text.set_visible(False)
for _, spine in ax.spines.items():
    spine.set_visible(True)
ax.set_ylabel("")
ax.set_title("{}".format(r"$\mathtt{Li_2S_4}$"), fontsize=20)
ax.set_xlabel("Cluster size", fontsize=20)
save_fig(fig, working_dir=results_dir, name="clusters.png")

cmap = custom_cmap()
new_df = pd.DataFrame(columns=size, data=[p9, p10, p11, p12, p13])
new_df["System"] = ["DOL/FLS7", "DOL/FLS8", "DOL/FLS9", "DOL/FLS10", "DOL/FLS11"]
new_df.set_index("System", inplace=True)

fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
p = sns.heatmap(
    new_df,
    cmap=cmap,
    annot=True,
    fmt=".2f",
    annot_kws={"size": 16},
    cbar_kws={"label": "Probability (size)"},
    linewidths=0.2,
    square=True,
    ax=ax,
)
ax.figure.axes[-1].yaxis.label.set_size(20)
ax.figure.axes[-1].yaxis.set_tick_params(labelsize=20)
ax.figure.axes[-1].yaxis.get_major_ticks()[0].set_visible(False)
ax.tick_params(axis="x", which="both", bottom=False, top=False, length=0)
ax.set_xticklabels(ax.get_xmajorticklabels(), fontsize=20)
ax.set_yticklabels(ax.get_ymajorticklabels(), fontsize=20)
for text in ax.texts:
    if text.get_text() == "0.00":
        text.set_visible(False)
for _, spine in ax.spines.items():
    spine.set_visible(True)
ax.set_ylabel("")
ax.set_title("{}".format(r"$\mathtt{Li_2S_4}$"), fontsize=20)
ax.set_xlabel("Cluster size", fontsize=20)
save_fig(fig, working_dir=results_dir, name="clusters_2.png")

# Diffusion
msd_int = pd.read_csv(f"{sys2_dir}/Diffusion/r1/msd_int.csv").reset_index()
msd_int["diff"] = msd_int["msd"] / (6 * (10 ** -15) * 50000 * 4)
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
dist_plot(ax, msd_int[msd_int["type"] == 2]["diff"] * 10 ** 9, 20, "DME", colors[1])
dist_plot(ax, msd_int[msd_int["type"] == 6]["diff"] * 10 ** 9, 40, "Li", colors[-1])
# ax.set_yticks(np.arange(0, 3.0, 0.5))
ax.set_xlabel(
    "Diffusivity ({} {})".format(r"$\mathtt{1O^{-9}}$", r"$\mathtt{m^2/s}$"),
    fontsize=24,
)
ax.set_ylabel("Frequency", fontsize=24)
save_fig(fig, working_dir=results_dir, name="diffusion_dist_dme.png")

msd_int = pd.read_csv(f"{sys5_dir}/Diffusion/r1/msd_int.csv").reset_index()
msd_int["diff"] = msd_int["msd"] / (6 * (10 ** -15) * 50000 * 4)
fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
dist_plot(ax, msd_int[msd_int["type"] == 2]["diff"] * 10 ** 9, 20, "FLS", colors[1])
dist_plot(ax, msd_int[msd_int["type"] == 6]["diff"] * 10 ** 9, 40, "Li", colors[-1])
# ax.set_yticks(np.arange(0, 3.0, 0.5))
ax.set_xlabel(
    "Diffusivity ({} {})".format(r"$\mathtt{1O^{-9}}$", r"$\mathtt{m^2/s}$"),
    fontsize=24,
)
ax.set_ylabel("Frequency", fontsize=24)
save_fig(fig, working_dir=results_dir, name="diffusion_dist_5.png")

# plotting diffusion
diff_sys1 = pd.read_csv(f"{sys1_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys2 = pd.read_csv(f"{sys2_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys3 = pd.read_csv(f"{sys3_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys4 = pd.read_csv(f"{sys4_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys5 = pd.read_csv(f"{sys5_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys6 = pd.read_csv(f"{sys6_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys7 = pd.read_csv(f"{sys7_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys8 = pd.read_csv(f"{sys8_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys9 = pd.read_csv(f"{sys9_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys10 = pd.read_csv(f"{sys10_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys11 = pd.read_csv(f"{sys11_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys12 = pd.read_csv(f"{sys12_dir}/Diffusion/diffusion.csv", index_col=0)
diff_sys13 = pd.read_csv(f"{sys13_dir}/Diffusion/diffusion.csv", index_col=0)

diff_list = [
    diff_sys2,
    diff_sys3,
    diff_sys4,
    diff_sys5,
    diff_sys6,
    diff_sys7,
    diff_sys8,
    diff_sys9,
    diff_sys10,
    diff_sys11,
    diff_sys12,
    diff_sys13,
]
diff = np.zeros((12, 6))
std = np.zeros((12, 6))

fig, ax = create_fig(figsize=(10, 8), set_axis=True, axis="both")
for i, d in enumerate(diff_list):
    diff[i] = d["diffusion (m2/s)"]
    std[i] = d["std"]
diff = diff.T
std = std.T
x_label = [
    "DOL/DME",
    "DOL/FLS1",
    "DOL/FLS2",
    "DOL/FLS3",
    "DOL/FLS4",
    "DOL/FLS5",
    "DOL/FLS6",
    "DOL/FLS7",
    "DOL/FLS8",
    "DOL/FLS9",
    "DOL/FLS10",
    "DOL/FLS11",
]
for i, j, c, l in zip(
    diff, std, colors, ["DOL", "Co-Solvent", "Additive", "Salt", "PS", "Li"]
):
    ax.errorbar(
        x_label,
        i * 10 ** 9,
        yerr=j * 10 ** 9,
        fmt="o",
        ecolor="k",
        capthick=2,
        elinewidth=2,
        mfc=c,
        ms=15,
        mew=2,
        mec="k",
        label=l,
    )

for i in np.arange(0.5, 11.0, 1):
    ax.axvline(x=i, ls="--", color="black")

ax.set_ylabel(
    "Diffusivity ({} {})".format(r"$\mathtt{1O^{-9}}$", r"$\mathtt{m^2/s}$"),
    fontsize=18,
)
ax.legend(
    loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=3, frameon=False, fontsize=18
)
ax.set_xticks(range(len(x_label)), x_label)
ax.set_xticklabels(x_label, rotation="vertical", fontsize=18)
save_fig(fig, working_dir=results_dir, name="diffusion.png")


# Conductivity
fig, ax = create_fig(figsize=(10, 8), set_axis=True, axis="both")
x_label = [
    "Base",
    "DOL/DME",
    "DOL/FLS1",
    "DOL/FLS2",
    "DOL/FLS3",
    "DOL/FLS4",
    "DOL/FLS5",
    "DOL/FLS6",
    "DOL/FLS7",
    "DOL/FLS8",
    "DOL/FLS9",
    "DOL/FLS10",
    "DOL/FLS11",
]

cond = [
    1.36,
    1.055,
    0.1324,
    0.2676,
    0.1096,
    0.1488,
    0.1758,
    0.255,
    0.237,
    0.1965,
    0.2183,
    0.2168,
    0.2673,
]

std = [
    0.3,
    0.42,
    0.047,
    0.053,
    0.0175,
    0.0125,
    0.0884,
    0.0565,
    0.0807,
    0.0453,
    0.0084,
    0.10699,
    0.0223,
]

ax.errorbar(
    x_label,
    np.array(cond) * 10,
    yerr=np.array(std) * 10,
    fmt="o",
    ecolor="k",
    capthick=2,
    elinewidth=2,
    mfc=colors[1],
    ms=15,
    mew=2,
    mec="k",
)

for i in np.arange(0.5, 12.0, 1):
    ax.axvline(x=i, ls="--", color="black")
ax.axhline(5, color=colors[1])
ax.set_ylabel(
    "Conductivity ({})".format(r"$\mathtt{mS/cm}$"), fontsize=18, color=colors[1]
)
ax.legend(
    loc="upper center", bbox_to_anchor=(0.5, 1.2), ncol=3, frameon=False, fontsize=18
)
ax.set_xticks(range(len(x_label)), x_label)
ax.set_xticklabels(x_label, rotation="vertical", fontsize=18)
ax2 = ax.twinx()
from ebatt.visualization.structural import _set_axis
_set_axis(ax2)
visc = [
    1.669,
    2.3117,
    3.996,
    2.8353,
    3.193,
    5.229,
    3.738,
    4.07,
    4.822,
    3.3533,
    3.394,
    6.751,
    4.1,
]

std = [
    0.0655,
    0.16464,
    0.1798,
    0.7884,
    0.1655,
    0.680,
    0.148,
    0.3902,
    0.0171,
    1.1139,
    0.1517,
    1.316,
    0.3512,
]
ax2.errorbar(
    x_label,
    visc,
    yerr=std,
    fmt="o",
    ecolor="k",
    capthick=2,
    elinewidth=2,
    mfc=colors[2],
    ms=15,
    mew=2,
    mec="k",
)
ax2.set_ylabel(
    "Viscosity ({})".format(r"$\mathtt{cP}$"), fontsize=18, color=colors[2]
)
ax2.axhline(2, color=colors[2])
ax.set_xticks(range(len(x_label)), x_label)
ax.set_xticklabels(x_label, rotation="vertical", fontsize=18)
save_fig(fig, working_dir=results_dir, name="cond.png")


def cond(num_mol, charges, diffusion, volume):
    T = 298.15
    e = 1.60217657e-19
    k = 1.3806488e-23
    cond = 0
    for n, q, d in zip(num_mol, charges, diffusion):
        cond += n * q ** 2 * d
    cond *= e ** 2 / k / T / volume
    return cond


dme_d = pd.read_csv(f"{analysis_dir}/DOL-DME/NVT/Diffusion/diffusion.csv", index_col=0)
dme_d = dme_d.loc["diffusion ($m^2$/s)"].values
dme_q = [
    0,
    0,
    -0.7198157507486944,
    -0.7198157507486946,
    -1.439631501497389,
    0.7198157507486945,
]
dme_n = [710, 478, 20, 130, 33, 216]
dme_v = 199181.05 * 10 ** -30
dme_cond = cond(dme_n, dme_q, dme_d, dme_v)

edf_d = pd.read_csv(f"{analysis_dir}/DOL-EDF/NVT/Diffusion/diffusion.csv", index_col=0)
edf_d = edf_d.loc["diffusion ($m^2$/s)"].values
edf_q = [
    0,
    0,
    -0.7198157507486944,
    -0.7198157507486946,
    -1.439631501497389,
    0.7198157507486945,
]
edf_n = [710, 478, 21, 130, 33, 217]
edf_v = 198488.83 * 10 ** -30
edf_cond = cond(edf_n, edf_q, edf_d, edf_v)

tte_d = pd.read_csv(f"{analysis_dir}/DOL-TTE/NVT/Diffusion/diffusion.csv", index_col=0)
tte_d = tte_d.loc["diffusion ($m^2$/s)"].values
tte_q = [
    0,
    0,
    -0.7474350927519358,
    -0.7474350927519361,
    -1.494870185503872,
    0.7474350927519359,
]
tte_n = [710, 328, 25, 130, 33, 221]
tte_v = 200970.9 * 10 ** -30
tte_cond = cond(tte_n, tte_q, tte_d, tte_v)

otp_d = pd.read_csv(f"{analysis_dir}/DOL-OTP/NVT/Diffusion/diffusion.csv", index_col=0)
otp_d = otp_d.loc["diffusion ($m^2$/s)"].values
otp_q = [
    0,
    0,
    -0.7474350927519358,
    -0.7474350927519361,
    -1.494870185503872,
    0.7474350927519359,
]
otp_n = [710, 269, 25, 130, 33, 221]
otp_v = 196768.54 * 10 ** -30
otp_cond = cond(otp_n, otp_q, otp_d, otp_v)

fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="y")
x = ["DOL:DME", "DOL:EDF", "DOL:TTE", "DOL:OTP"]
y = np.array([dme_cond, edf_cond, tte_cond, otp_cond]) / dme_cond
ax.bar(x, y, color="#DDDDDD", edgecolor="black", linewidth=2, align="center", width=0.5)
ax.set_yticks(np.arange(0, 1.2, 0.2))
ax.set_xticklabels(x, ma="center", ha="center")
ax.tick_params(axis="x", direction="in", labelsize=24)
ax.set_ylabel("Relative conductivity", fontsize=24)
ax.xaxis.set_tick_params(pad=5)
save_fig(fig, working_dir=results_dir, name="cond.png")

# Residence Time
dme_r = pd.read_csv(
    f"{analysis_dir}/DOL-DME/NVT/Residence_Time/" f"auto_correlation.csv"
)
edf_r = pd.read_csv(
    f"{analysis_dir}/DOL-EDF/NVT/Residence_Time/" f"auto_correlation.csv"
)
tte_r = pd.read_csv(
    f"{analysis_dir}/DOL-TTE/NVT/Residence_Time/" f"auto_correlation.csv"
)
otp_r = pd.read_csv(
    f"{analysis_dir}/DOL-OTP/NVT/Residence_Time/" f"auto_correlation.csv"
)

fig, ax = create_fig(figsize=(8, 6), set_axis=True, axis="both")
ax.plot(
    dme_r["Time (ps)"],
    dme_r["15-13"],
    label="DME",
    linewidth=2,
    color=palette.as_hex()[0],
)
ax.plot(
    edf_r["Time (ps)"],
    edf_r["15-13"],
    label="EDF",
    linewidth=2,
    color=palette.as_hex()[1],
)
ax.plot(
    tte_r["Time (ps)"],
    tte_r["15-13"],
    label="TTE",
    linewidth=2,
    color=palette.as_hex()[2],
)
ax.plot(
    otp_r["Time (ps)"],
    otp_r["15-13"],
    label="OTP",
    linewidth=2,
    color=palette.as_hex()[3],
)
ax.legend(fontsize=24, frameon=False)
ax.xaxis.set_tick_params(pad=5)
ax.set_xlabel("Time (ps)", fontsize=24)
ax.set_ylabel("C(t)", fontsize=24)
save_fig(fig, working_dir=results_dir, name="autocorrelation.png")
