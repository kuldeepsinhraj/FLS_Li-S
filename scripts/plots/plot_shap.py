import shap

import numpy as np
import pandas as pd
import xgboost as xgb
import matplotlib.pyplot as plt

from ebatt.descriptors.structural import rdkit_descriptors, get_Mordred_descriptors


def run_shap(
    data,
    x_cols,
    y_col,
    n_estimators=100,
    max_depth=2,
    reg_alpha=0.00001,
    reg_lambda=0.00001,
    max_display=15,
    plot_feature_dependence=False,
    cmap="RdYlBu",
    filename="shap",
):
    # get data
    x_data = data[x_cols]
    y_data = data[y_col]

    # train XGBoost model
    model = xgb.XGBRegressor(
        n_estimators=n_estimators,
        max_depth=max_depth,
        reg_alpha=reg_alpha,
        reg_lambda=reg_lambda,
    )
    model.fit(x_data, y_data)

    # explain model predictions using SHAP values
    explainer = shap.TreeExplainer(model)
    shap_values = explainer.shap_values(x_data)

    # plot the summary of feature importance
    plt.figure()
    shap.summary_plot(
        shap_values,
        x_data,
        max_display=max_display,
        show=False,
        plot_size=0.2,
        color_bar=False,
    )
    for fc in plt.gcf().get_children():
        for fcc in fc.get_children():
            if hasattr(fcc, "set_cmap"):
                fcc.set_cmap(cmap)
    ax = plt.gca()
    ax.tick_params(
        axis="x", which="both", direction="in", length=10, labelsize=14, pad=10
    )
    ax.tick_params(axis="y", which="both", direction="in", length=10, pad=10)
    ax.set_xlabel("SHAP Value", fontsize=14, fontweight="bold")

    ytick_labels = [label.get_text() for label in ax.get_yticklabels()]

    label_mapping = {
        "R": "Radius (Å)",
        "B1": "B1 (Å)",
        "B5": "B5 (Å)",
        "Polarizability": "Polarizability (Å³)",
        "FDeg": "Fluorination Degree (%)",
        "ProximityO": "Proximity of F to O (Å)",
        "StericEffect": "Steric Effect",
        "ShapeInd": "Shape Index",
        "vdWV": "Van der Waals Volume (Å³)",
        "MolecularWeight": "Molecular Weight (g/mol)",
        "Dipole": "Dipole Moment (Debye)",
        "MinPartialCharge": "Min Partial Charge",
        "MaxPartialCharge": "Max Partial Charge",
        "MinAbsPartialCharge": "Min Abs Partial Charge",
        "MaxAbsPartialCharge": "Max Abs Partial Charge",
    }
    ytick_labels_custom = [
        label_mapping[label] if label in label_mapping else label
        for label in ytick_labels
    ]
    ax.set_yticklabels(ytick_labels_custom)
    x_limits = ax.get_xlim()
    limits = max(abs(x_limits[0]), abs(x_limits[1]))
    # ax.set_xticks([-limits, 0.0, limits])
    xticks = np.arange(-limits, limits, limits / 2)
    if limits not in xticks:
        xticks = np.append(np.arange(-limits, limits, limits / 2), limits)
    ax.set_xticks(xticks)
    # ax.set_xticklabels([f"{-limits:.1f}", "0.0", f"{limits:.1f}"], fontsize=14)
    ax.set_xticklabels([f"{i:.1f}" for i in xticks], fontsize=14)
    cbar = plt.colorbar()
    cbar.ax.tick_params(labelsize=16)
    cbar.ax.tick_params(direction="out", length=10, grid_alpha=0.5)
    y_limits = cbar.ax.get_ylim()
    cbar.ax.set_yticks([y_limits[0], y_limits[1]])
    cbar.ax.set_yticklabels(["Low", "High"], fontsize=14)
    cbar.ax.set_ylabel("Feature Value", fontsize=14, fontweight="bold", labelpad=-5)
    for line in ax.lines:
        if line.get_label() == "_child0":
            line.set_linewidth(2.5)
            line.set_color("black")
    ax.spines["bottom"].set_linewidth(2)
    plt.savefig(
        f"../../files/results/{filename}.pdf",
        dpi=1200,
        bbox_inches="tight",
        pad_inches=0.1,
    )

    # plot the dependence of each feature on the model output
    if plot_feature_dependence:
        for i in range(x_data.shape[1]):
            shap.dependence_plot(f"rank({i})", shap_values, x_data)


df = pd.read_csv("../../files/results/solvents.csv")
df = df[~df.isna().any(axis=1)]

estimators = 20
depth = 3
alpha = 0.01
lambda_ = 0.01
display = 10
plot_dependence = False
color_map = "RdYlBu"

# basic features
x = [
    "MolecularWeight",
    "Complexity",
    # "RotatableBondCount",
    # "Volume3D",
    # "XStericQuadrupole3D",
    # "YStericQuadrupole3D",
    # "ZStericQuadrupole3D",
    "R",
    # "GyrationR",
    "vdWV",
    "ShapeInd",
    "StericEffect",
    "ProximityO",
    "FDeg",
    "ElectroCount",
    # "L",
    "B1",
    "B5",
    "Dipole",
    "Polarizability",
    "MinPartialCharge",
    "MaxPartialCharge",
    "MinAbsPartialCharge",
    "MaxAbsPartialCharge",
]

for y in ["BE_Salt", "BE_PS", "IP", "EA", "EW"]:
    run_shap(
        df,
        x,
        y,
        n_estimators=estimators,
        max_depth=depth,
        reg_alpha=alpha,
        reg_lambda=lambda_,
        max_display=display,
        plot_feature_dependence=plot_dependence,
        cmap=color_map,
        filename=f"shap_{y}",
    )

# rdkit features
mol_descriptors, desc_names = rdkit_descriptors(df["CanonicalSMILES"])
rdkit_df = pd.DataFrame(mol_descriptors, columns=desc_names)
# fr_cols = [col for col in rdkit_df.columns if col.startswith("fr_")]  # mostly zeros
# rdkit_df.drop(fr_cols, axis=1, inplace=True)
rdkit_df.drop(
    [
        "MaxPartialCharge",
        "MinPartialCharge",
        "MaxAbsPartialCharge",
        "MinAbsPartialCharge",
    ],
    axis=1,
    inplace=True,
)  # remove charges that are already included in the basic features
del rdkit_df["MolWt"]  # repeated
del rdkit_df["ExactMolWt"]  # repeated?
# del rdkit_df["NumRotatableBonds"]  # repeated
# del rdkit_df["NumRadicalElectrons"]  # mostly zero
# del rdkit_df["Ipc"]  # mostly zero
final_df = df.merge(rdkit_df, left_index=True, right_index=True)
final_df = final_df.dropna().reset_index(drop=True)

x_rdkit = x + list(rdkit_df.columns)

for y in ["BE_Salt", "BE_PS", "IP", "EA", "EW"]:
    run_shap(
        final_df,
        x_rdkit,
        y,
        n_estimators=estimators,
        max_depth=depth,
        reg_alpha=alpha,
        reg_lambda=lambda_,
        max_display=display,
        plot_feature_dependence=plot_dependence,
        cmap=color_map,
        filename=f"shap_rdkit_{y}",
    )

# mordred features
mordred_df = get_Mordred_descriptors(df["CanonicalSMILES"])
mordred_df = mordred_df.select_dtypes(include=[np.number])
del mordred_df["MW"]
final_df = df.merge(mordred_df, left_index=True, right_index=True)

x_mordred = x + list(mordred_df.columns)

for y in ["BE_Salt", "BE_PS", "IP", "EA", "EW"]:
    run_shap(
        final_df,
        x_mordred,
        y,
        n_estimators=estimators,
        max_depth=depth,
        reg_alpha=alpha,
        reg_lambda=lambda_,
        max_display=display,
        plot_feature_dependence=plot_dependence,
        cmap=color_map,
        filename=f"shap_mordred_{y}",
    )


# # lasso regression
# from sklearn.linear_model import Lasso
#
# ejer = final_df[x]
# ejer = (ejer - ejer.mean()) / np.linalg.norm(ejer - ejer.mean(), axis=0)
# model2 = Lasso(alpha=0.00001)
# model2.fit(
#     ejer.fillna(0),
#     (final_df["IP"] - final_df["IP"].mean()) / np.linalg.norm(final_df["IP"]),
# )
# pd.Series(model2.coef_, index=ejer.columns).abs().sort_values().tail(10)
