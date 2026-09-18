#!/usr/bin/env python3
"""

 FLS property prediction - model training, evaluation, and figure generation

"""
import argparse, os, warnings
warnings.filterwarnings("ignore")

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from sklearn.model_selection import train_test_split, KFold, cross_val_predict
from sklearn.metrics import r2_score, mean_absolute_error, root_mean_squared_error
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
import shap

RANDOM_STATE = 42
DATA_PATH = "files/results/solvents.csv"
CACHE_DIR = "descriptor_cache"
DROP_EA_OUTLIERS = False

RADICAL_CIDS = [58537745, 58537760, 58700661, 59344306,
                58544205, 58544217, 101084693, 101938351]

FEATURES = ["MolecularWeight", "Complexity", "R", "vdWV", "ShapeInd", "StericEffect",
            "ProximityO", "FDeg", "ElectroCount", "B1", "B5", "Dipole",
            "Polarizability", "MinPartialCharge", "MaxPartialCharge",
            "MinAbsPartialCharge", "MaxAbsPartialCharge"]

TARGETS = {
    "IP":      dict(label="Ionization Potential",    axis="IP (eV)", units="eV"),
    "EA":      dict(label="Electron Affinity",       axis="EA (eV)", units="eV"),
    "BE_Salt": dict(label="Binding Energy - LiTFSI", axis="BE-LiTFSI (kcal/mol)",
                    units="kcal/mol", drop_above=0.0),
    "BE_PS":   dict(label="Binding Energy - Li2S4",  axis="BE-Li2S4 (kcal/mol)",
                    units="kcal/mol"),
}


PERF_SET = {"IP": "mordred", "EA": "rdkit", "BE_Salt": "mordred", "BE_PS": "mordred"}
SHAP_SET = {t: "custom" for t in ("IP", "EA", "BE_Salt", "BE_PS")}

# Backwards-compatibility alias; not used directly by any function below.
BEST_CONFIG = PERF_SET

XGB_PARAMS = dict(n_estimators=500, max_depth=5, learning_rate=0.04, subsample=0.8,
                  colsample_bytree=0.6, reg_lambda=2.0, min_child_weight=2)
RF_PARAMS = dict(n_estimators=300, max_depth=None, min_samples_leaf=1, max_features=0.3)

COLOR_TRAIN_FACE, COLOR_TRAIN_EDGE = "#BFD3F2", "#5B8FD9"
COLOR_TEST_FACE, COLOR_TEST_EDGE = "#E8743B", "#9C3D14"
SHAP_CMAP = "RdYlBu"

LABEL_MAP = {
    "R": "Radius (A)", "B1": "B1 (A)", "B5": "B5 (A)",
    "Polarizability": "Polarizability (A3)", "FDeg": "Fluorination Degree (%)",
    "ProximityO": "Proximity of F to O (A)", "StericEffect": "Steric Effect",
    "ShapeInd": "Shape Index", "vdWV": "Van der Waals Volume (A3)",
    "MolecularWeight": "Molecular Weight (g/mol)", "Dipole": "Dipole Moment (Debye)",
    "MinPartialCharge": "Min Partial Charge", "MaxPartialCharge": "Max Partial Charge",
    "MinAbsPartialCharge": "Min Abs Partial Charge",
    "MaxAbsPartialCharge": "Max Abs Partial Charge",
    "ElectroCount": "Electronegative Atom Count", "Complexity": "Complexity",
}

SET_LABEL = {"custom": "Custom (17)", "rdkit": "Custom + RDKit",
             "mordred": "Custom + Mordred"}


def _clean(d):
    d = d.replace([np.inf, -np.inf], np.nan)
    d = d.loc[:, d.isna().mean() < 0.05]
    return d.loc[:, d.std() > 1e-8]


def _prune_correlated(d, thresh=0.95):
    filled = d.fillna(d.median())
    corr = filled.corr().abs()
    upper = corr.where(np.triu(np.ones(corr.shape), 1).astype(bool))
    return d.drop(columns=[c for c in upper.columns if any(upper[c] > thresh)])


def build_descriptors():
    from rdkit import Chem
    from rdkit.Chem import Descriptors
    from rdkit.ML.Descriptors import MoleculeDescriptors
    from mordred import Calculator, descriptors as mordred_desc

    os.makedirs(CACHE_DIR, exist_ok=True)
    df = pd.read_csv(DATA_PATH)
    df = df[~df.CID.isin(RADICAL_CIDS)].reset_index(drop=True)
    mols = [Chem.MolFromSmiles(s) for s in df.CanonicalSMILES]

    calc = MoleculeDescriptors.MolecularDescriptorCalculator(
        [d[0] for d in Descriptors._descList])
    names = calc.GetDescriptorNames()
    rows = [calc.CalcDescriptors(Chem.AddHs(m)) if m else [np.nan] * len(names)
            for m in mols]
    rd = pd.DataFrame(rows, columns=names)
    for c in ["MaxPartialCharge", "MinPartialCharge", "MaxAbsPartialCharge",
              "MinAbsPartialCharge", "MolWt", "ExactMolWt"]:
        if c in rd.columns:
            del rd[c]
    rd = _clean(rd)
    rd.to_csv(f"{CACHE_DIR}/rdkit.csv", index=False)
    print(f"  RDKit descriptors: {rd.shape[1]}")

    mo = Calculator(mordred_desc, ignore_3D=False).pandas(mols)
    mo = mo.select_dtypes(include=[np.number])
    if "MW" in mo.columns:
        del mo["MW"]
    mo = _clean(mo)
    mo.to_csv(f"{CACHE_DIR}/mordred.csv", index=False)
    df.to_csv(f"{CACHE_DIR}/base.csv", index=False)
    print(f"  Mordred descriptors: {mo.shape[1]}")


def load_sets():
    if not os.path.exists(f"{CACHE_DIR}/mordred.csv"):
        print("Descriptor cache not found - building it now.")
        build_descriptors()
    base = pd.read_csv(f"{CACHE_DIR}/base.csv")
    rdk = pd.read_csv(f"{CACHE_DIR}/rdkit.csv")
    mrd = _prune_correlated(pd.read_csv(f"{CACHE_DIR}/mordred.csv"))
    return {
        "custom": (base, FEATURES),
        "rdkit": (pd.concat([base, rdk], axis=1), FEATURES + list(rdk.columns)),
        "mordred": (pd.concat([base, mrd], axis=1), FEATURES + list(mrd.columns)),
    }


def get_xy(target, sets, which):
    data, cols = sets[which]
    cfg = TARGETS[target]
    d = data[cols + [target]].dropna()
    if "drop_above" in cfg:
        d = d[d[target] <= cfg["drop_above"]]
    if target == "EA" and DROP_EA_OUTLIERS:
        d = d[d[target] <= 3.3]
    return d[cols], d[target]


def make_model(name, seed=RANDOM_STATE):
    if name == "XGBoost":
        return xgb.XGBRegressor(random_state=seed, n_jobs=-1,
                                tree_method="hist", **XGB_PARAMS)
    if name == "Random forest":
        return RandomForestRegressor(random_state=seed, n_jobs=-1, **RF_PARAMS)
    return Pipeline([("scale", StandardScaler()),
                     ("model", RidgeCV(alphas=np.logspace(-3, 3, 25)))])


def split(X, y):
    strata = pd.qcut(y, 5, labels=False, duplicates="drop")
    return train_test_split(X, y, test_size=0.20,
                            random_state=RANDOM_STATE, stratify=strata)


def compare(targets, outdir):
    sets = load_sets()
    rows = []
    for t in targets:
        for which in ["custom", "rdkit", "mordred"]:
            X, y = get_xy(t, sets, which)
            for mname in ["XGBoost", "Random forest", "Ridge"]:
                cvp = cross_val_predict(make_model(mname), X, y,
                                        cv=KFold(5, shuffle=True,
                                                 random_state=RANDOM_STATE))
                X_tr, X_te, y_tr, y_te = split(X, y)
                p = make_model(mname).fit(X_tr, y_tr).predict(X_te)
                rows.append(dict(
                    target=t, property=TARGETS[t]["label"], units=TARGETS[t]["units"],
                    descriptor_set=SET_LABEL[which], model=mname,
                    n_features=len(X.columns), N=len(X),
                    cv_R2=r2_score(y, cvp), cv_MAE=mean_absolute_error(y, cvp),
                    cv_RMSE=root_mean_squared_error(y, cvp),
                    test_R2=r2_score(y_te, p), test_MAE=mean_absolute_error(y_te, p),
                    test_RMSE=root_mean_squared_error(y_te, p)))
                print(f"  {t:<9}{SET_LABEL[which]:<18}{mname:<15}"
                      f"cv={rows[-1]['cv_R2']:.3f} test={rows[-1]['test_R2']:.3f}")
    pd.DataFrame(rows).to_csv(f"{outdir}/model_comparison.csv", index=False)


def _style_shap_axes(ax):
    """Apply the established SHAP plot styling (shared by all SHAP plots)."""
    ax.tick_params(axis="x", which="both", direction="in", length=10,
                   labelsize=14, pad=10)
    ax.tick_params(axis="y", which="both", direction="in", length=10, pad=10)
    ax.set_xlabel("SHAP Value", fontsize=14, fontweight="bold")
    ax.set_yticklabels([LABEL_MAP.get(l.get_text(), l.get_text())
                        for l in ax.get_yticklabels()])
    xl = ax.get_xlim()
    lim = max(abs(xl[0]), abs(xl[1]))
    xt = np.arange(-lim, lim, lim / 2)
    if lim not in xt:
        xt = np.append(xt, lim)
    ax.set_xticks(xt)
    ax.set_xticklabels([f"{i:.1f}" for i in xt], fontsize=14)
    cb = plt.colorbar()
    cb.ax.tick_params(labelsize=16)
    cb.ax.tick_params(direction="out", length=10, grid_alpha=0.5)
    yl = cb.ax.get_ylim()
    cb.ax.set_yticks([yl[0], yl[1]])
    cb.ax.set_yticklabels(["Low", "High"], fontsize=14)
    cb.ax.set_ylabel("Feature Value", fontsize=14, fontweight="bold", labelpad=-5)
    for line in ax.lines:
        if line.get_label() == "_child0":
            line.set_linewidth(2.5)
            line.set_color("black")
    ax.spines["bottom"].set_linewidth(2)


def parity(targets, outdir):
    """One parity plot per property, using each property's best configuration
    (PERF_SET). R2 shown on these plots is held-out TEST R2."""
    sets = load_sets()
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 1.6,
                         "xtick.direction": "in", "ytick.direction": "in",
                         "xtick.major.size": 6, "ytick.major.size": 6,
                         "xtick.major.width": 1.4, "ytick.major.width": 1.4})
    for t in targets:
        cfg = TARGETS[t]
        X, y = get_xy(t, sets, PERF_SET[t])
        X_tr, X_te, y_tr, y_te = split(X, y)
        m = make_model("XGBoost").fit(X_tr, y_tr)
        p_tr, p_te = m.predict(X_tr), m.predict(X_te)

        fig, ax = plt.subplots(figsize=(4.4, 4.4))
        lo = min(y.min(), p_tr.min(), p_te.min())
        hi = max(y.max(), p_tr.max(), p_te.max())
        pad = 0.06 * (hi - lo)
        lo, hi = lo - pad, hi + pad
        ax.plot([lo, hi], [lo, hi], color="0.35", lw=1.4, ls="--", zorder=1)
        ax.scatter(y_tr, p_tr, s=26, facecolor=COLOR_TRAIN_FACE,
                   edgecolor=COLOR_TRAIN_EDGE, linewidth=0.5, alpha=0.75,
                   label="Train", zorder=2)
        ax.scatter(y_te, p_te, s=34, facecolor=COLOR_TEST_FACE,
                   edgecolor=COLOR_TEST_EDGE, linewidth=0.6, alpha=0.9,
                   label="Test", zorder=3)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_xlabel(f"DFT-computed {cfg['axis']}", fontsize=11,
                      fontweight="bold", labelpad=7)
        ax.set_ylabel(f"Predicted {cfg['axis']}", fontsize=11,
                      fontweight="bold", labelpad=7)
        ax.xaxis.set_major_locator(MaxNLocator(5))
        ax.yaxis.set_major_locator(MaxNLocator(5))
        ax.tick_params(labelsize=10, pad=6)
        ax.text(0.04, 0.96,
                f"$R^2$ = {r2_score(y_te, p_te):.3f}\n"
                f"MAE = {mean_absolute_error(y_te, p_te):.3f}\n"
                f"RMSE = {root_mean_squared_error(y_te, p_te):.3f}",
                transform=ax.transAxes, fontsize=9.5, va="top",
                bbox=dict(boxstyle="round,pad=0.45", fc="white", ec="0.7", lw=0.9))
        ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.03), fontsize=9,
                  frameon=False, handletextpad=0.4)
        ax.set_title(cfg["label"], fontsize=12, fontweight="bold", pad=10)
        plt.tight_layout()
        for ext in ("png", "pdf"):
            plt.savefig(f"{outdir}/parity_{t}.{ext}",
                        dpi=600 if ext == "png" else None, bbox_inches="tight")
        plt.close()
        print(f"  {t}: parity written ({PERF_SET[t]} descriptors)")


def shap_plots(targets, outdir, max_display=10):

    sets = load_sets()
    plt.rcParams.update({"font.family": "DejaVu Sans", "axes.linewidth": 1.6,
                         "xtick.direction": "in", "ytick.direction": "in",
                         "xtick.major.size": 6, "ytick.major.size": 6,
                         "xtick.major.width": 1.4, "ytick.major.width": 1.4})
    for t in targets:
        for which in ["custom", "rdkit", "mordred"]:
            X, y = get_xy(t, sets, which)
            m = make_model("XGBoost").fit(X, y)
            sv = shap.TreeExplainer(m).shap_values(X)

            plt.figure()
            shap.summary_plot(sv, X, max_display=max_display, show=False,
                              plot_size=0.2, color_bar=False)
            for fc in plt.gcf().get_children():
                for fcc in fc.get_children():
                    if hasattr(fcc, "set_cmap"):
                        fcc.set_cmap(SHAP_CMAP)
            _style_shap_axes(plt.gca())
            for ext in ("png", "pdf"):
                plt.savefig(f"{outdir}/shap_{t}_{which}.{ext}", dpi=1200,
                            bbox_inches="tight", pad_inches=0.1)
            plt.close()
            print(f"  {t} / {SET_LABEL[which]}: shap_{t}_{which} written "
                  f"({len(X.columns)} features)")


def figures(targets, outdir):
    """Convenience wrapper: runs parity() then shap_plots()."""
    parity(targets, outdir)
    shap_plots(targets, outdir)


def tables(targets, outdir, n_seeds=10):
    sets = load_sets()
    s6, s8 = [], []
    for t in targets:
        X, y = get_xy(t, sets, PERF_SET[t])
        r2s = [r2_score(y, cross_val_predict(
            make_model("XGBoost", seed=s), X, y,
            cv=KFold(5, shuffle=True, random_state=s))) for s in range(n_seeds)]
        p0 = cross_val_predict(make_model("XGBoost", seed=0), X, y,
                               cv=KFold(5, shuffle=True, random_state=0))
        s6.append(dict(Property=TARGETS[t]["label"], N=len(X),
                       Units=TARGETS[t]["units"],
                       Descriptors=SET_LABEL[PERF_SET[t]],
                       CV_R2=f"{np.mean(r2s):.3f} +/- {np.std(r2s):.3f}",
                       CV_MAE=f"{mean_absolute_error(y, p0):.2f}",
                       CV_RMSE=f"{root_mean_squared_error(y, p0):.2f}"))
        print(f"  S6 {t}: {np.mean(r2s):.3f} +/- {np.std(r2s):.3f}")

        ranks = {f: [] for f in FEATURES}
        # S8 uses SHAP_SET, not PERF_SET: this table validates the SHAP
        # figures, so it must describe the same model those figures come from.
        Xs, ys = get_xy(t, sets, SHAP_SET[t])
        idx = [Xs.columns.get_loc(f) for f in FEATURES]
        for s in range(n_seeds):
            m = make_model("XGBoost", seed=s).fit(Xs, ys)
            imp = np.abs(shap.TreeExplainer(m).shap_values(Xs)[:, idx]).mean(0)
            order = [FEATURES[i] for i in np.argsort(-imp)]
            for f in FEATURES:
                ranks[f].append(order.index(f) + 1)
        for f, r in sorted(ranks.items(), key=lambda kv: np.mean(kv[1]))[:6]:
            s8.append(dict(Property=TARGETS[t]["label"],
                           Descriptor_set=SET_LABEL[SHAP_SET[t]],
                           Descriptor=LABEL_MAP.get(f, f),
                           Mean_rank=f"{np.mean(r):.1f}",
                           Rank_range=f"{min(r)} - {max(r)}"))
        print(f"  S8 {t} done ({SHAP_SET[t]} descriptors)")
    pd.DataFrame(s6).to_csv(f"{outdir}/TableS6_cv_performance.csv", index=False)
    pd.DataFrame(s8).to_csv(f"{outdir}/TableS8_shap_rank_stability.csv", index=False)


def depth_sensitivity(targets, outdir, depths=(2, 3, 4, 5, 6, 8, 10)):

    sets = load_sets()
    rows = []
    for t in targets:
        X, y = get_xy(t, sets, SHAP_SET[t])
        for d in depths:
            params = dict(XGB_PARAMS)
            params["max_depth"] = d
            model = xgb.XGBRegressor(random_state=RANDOM_STATE, n_jobs=-1,
                                     tree_method="hist", **params)
            cvp = cross_val_predict(model, X, y,
                                    cv=KFold(5, shuffle=True,
                                             random_state=RANDOM_STATE))
            model.fit(X, y)
            rows.append(dict(
                Property=TARGETS[t]["label"],
                Descriptor_set=SET_LABEL[SHAP_SET[t]],
                Max_tree_depth=d,
                CV_R2=f"{r2_score(y, cvp):.3f}",
                InSample_R2=f"{r2_score(y, model.predict(X)):.3f}"))
        print(f"  S7 {t} done")
    pd.DataFrame(rows).to_csv(f"{outdir}/TableS7_depth_sensitivity.csv", index=False)


def cross_model(targets, outdir):

    from scipy.stats import spearmanr
    sets = load_sets()
    detail_rows, summary_rows = [], []
    for t in targets:
        X, y = get_xy(t, sets, SHAP_SET[t])
        idx = [X.columns.get_loc(f) for f in FEATURES]
        Xs = X.sample(min(400, len(X)), random_state=0)

        ranks = {}
        for name in ["XGBoost", "Random forest"]:
            m = make_model(name).fit(X, y)
            sv = shap.TreeExplainer(m).shap_values(Xs)
            imp = np.abs(np.asarray(sv)[:, idx]).mean(0)
            order = [FEATURES[i] for i in np.argsort(-imp)]
            ranks[name] = {f: order.index(f) + 1 for f in FEATURES}

        for f in FEATURES:
            detail_rows.append(dict(
                Property=TARGETS[t]["label"], Descriptor=LABEL_MAP.get(f, f),
                XGBoost_rank=ranks["XGBoost"][f],
                RandomForest_rank=ranks["Random forest"][f]))

        xgb_r = [ranks["XGBoost"][f] for f in FEATURES]
        rf_r = [ranks["Random forest"][f] for f in FEATURES]
        rho, _ = spearmanr(xgb_r, rf_r)
        top3_x = set(sorted(FEATURES, key=lambda f: ranks["XGBoost"][f])[:3])
        top3_r = set(sorted(FEATURES, key=lambda f: ranks["Random forest"][f])[:3])
        summary_rows.append(dict(Property=TARGETS[t]["label"],
                                 Spearman_rho=f"{rho:.2f}",
                                 Top3_overlap=f"{len(top3_x & top3_r)} of 3"))
        print(f"  {t}: Spearman rho={rho:.2f}, top-3 overlap="
              f"{len(top3_x & top3_r)}/3")

    pd.DataFrame(detail_rows).to_csv(
        f"{outdir}/TableS9_cross_model_detail.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(
        f"{outdir}/TableS9_cross_model_summary.csv", index=False)


def summary_table(targets, outdir, top_n=6):

    sets = load_sets()
    cols = {}
    for t in targets:
        X, y = get_xy(t, sets, SHAP_SET[t])
        idx = [X.columns.get_loc(f) for f in FEATURES]
        m = make_model("XGBoost").fit(X, y)
        sv = shap.TreeExplainer(m).shap_values(X)
        imp = np.abs(np.asarray(sv)[:, idx]).mean(0)
        order = [FEATURES[i] for i in np.argsort(-imp)][:top_n]
        cols[TARGETS[t]["label"]] = [LABEL_MAP.get(f, f) for f in order]
        print(f"  {t}: top {top_n} computed")
    pd.DataFrame(cols).to_csv(f"{outdir}/Table_summary_top_descriptors.csv",
                              index=False)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--descriptors", action="store_true")
    ap.add_argument("--compare", action="store_true")
    ap.add_argument("--figures", action="store_true",
                    help="parity plots + all 12 SHAP plots (parity + shap)")
    ap.add_argument("--parity", action="store_true", help="parity plots only")
    ap.add_argument("--shap", action="store_true",
                    help="all 12 SHAP plots only (3 per property x 4 properties)")
    ap.add_argument("--tables", action="store_true",
                    help="Tables S6 and S8 (seed-repeated CV performance and "
                         "SHAP rank stability, custom-17 descriptors only)")
    ap.add_argument("--depth", action="store_true",
                    help="Table S7 (tree-depth sensitivity)")
    ap.add_argument("--crossmodel", action="store_true",
                    help="Table S9 (XGBoost vs Random Forest SHAP agreement)")
    ap.add_argument("--summary", action="store_true",
                    help="simple cross-property top-descriptor table, no "
                         "numeric ranks shown (optional; your call whether "
                         "to keep it)")
    ap.add_argument("--target", default=None)
    ap.add_argument("--outdir", default="ml_output")
    a = ap.parse_args()

    if not os.path.exists(DATA_PATH):
        raise SystemExit(f"Cannot find {DATA_PATH}. Run from the repository root.")
    os.makedirs(a.outdir, exist_ok=True)
    targets = [a.target] if a.target else list(TARGETS)

    if a.descriptors:
        print("\nComputing descriptors")
        build_descriptors()
    if a.all or a.compare:
        print("\nModel and descriptor-set comparison")
        compare(targets, a.outdir)
    if a.all or a.figures or a.parity:
        print("\nParity plots")
        parity(targets, a.outdir)
    if a.all or a.figures or a.shap:
        print("\nSHAP plots (3 per property)")
        shap_plots(targets, a.outdir)
    if a.all or a.tables:
        print("\nTables S6 and S8")
        tables(targets, a.outdir)
    if a.all or a.depth:
        print("\nTable S7: depth sensitivity")
        depth_sensitivity(targets, a.outdir)
    if a.all or a.crossmodel:
        print("\nTable S9: cross-model agreement")
        cross_model(targets, a.outdir)
    if a.all or a.summary:
        print("\nSummary table (optional)")
        summary_table(targets, a.outdir)
    ran = [a.all, a.compare, a.figures, a.parity, a.shap, a.tables, a.depth,
           a.crossmodel, a.summary, a.descriptors]
    if not any(ran):
        ap.print_help()
    else:
        print(f"\nDone. Outputs in {a.outdir}/")


if __name__ == "__main__":
    main()
