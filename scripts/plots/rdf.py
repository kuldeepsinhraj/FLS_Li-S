import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

from ebatt.visualization.structural import custom_cmap, save_fig


rdf_dir = "/Users/rashatwi/Desktop/fls/files/analysis/rdf"
rdf_files = ["dme.csv"] + ["tte.csv"] + [f"fls{i}.csv" for i in range(0, 19)]
cmap = custom_cmap(
    ["#C7CAFA", "#FFFFFF", "#F9554A"],
    numcolors=50,
)


def format_rdf_ax(axs, xticks, labels):
    for ax in axs:
        ax.figure.axes[-1].xaxis.label.set_size(56)
        ax.figure.axes[-2].xaxis.label.set_size(56)
        ax.figure.axes[-3].xaxis.label.set_size(56)
        ax.figure.axes[-3].tick_params(labelsize=46, direction="in", length=10, pad=10)
        ax.figure.axes[-2].tick_params(labelsize=46, direction="in", length=10, pad=10)
        ax.figure.axes[-1].tick_params(labelsize=46, direction="in", length=10, pad=10)
        for _, spine in ax.spines.items():
            spine.set_visible(True)
        for axis in ["top", "bottom", "left", "right"]:
            ax.spines[axis].set_linewidth(5)
        ax.tick_params(axis="x", which="both", direction="in", length=20, pad=10)
        ax.tick_params(
            axis="y",
            which="both",
            bottom=False,
            top=False,
            direction="in",
            length=20,
            pad=10,
        )
        # ax.set_xlim([2, 8])
        ax.set_xticks(xticks)
        ax.set_xticklabels(labels, fontsize=46, rotation="horizontal")
        ax.set_yticklabels(
            ax.get_ymajorticklabels(), fontsize=36, rotation="horizontal"
        )
        ax.set_xlabel(r"$\mathbf{r\ (\AA)}$", weight="bold", fontsize=56)
        ax.set_ylabel("")
        cbar = ax.collections[0].colorbar
        cbar.outline.set_linewidth(5)
        cbar.set_label("g(r)", weight="bold", fontsize=56, labelpad=40)
    axs[1].set_yticklabels("")
    axs[2].set_yticklabels("")


def x_axis_ticks_labels(x_values, num_of_ticks):
    val_range = [int(min(x_values)) + 1, int(max(x_values)) + 1]
    step_size = (val_range[1] - val_range[0]) / num_of_ticks
    if step_size < 0.5:
        ticks_step_size = 0.5
    else:
        ticks_step_size = int(step_size) or 1
    print(val_range, ticks_step_size)
    x_step_size = x_values[1] - x_values[0]
    labels_ = list(np.arange(*val_range, ticks_step_size))
    x_ticks = []
    for label in labels_:
        for x_ind, x_value in enumerate(x_values):
            if label > x_value:
                x_ticks.append(x_ind + (label - x_value) / x_step_size)
                break
    return x_ticks, labels_


# St (PS) - St (PS)
fig, axs = plt.subplots(1, 3, figsize=(30, 15))
for idx, sys in enumerate(["0.25_li2s2", "0.25_li2s4", "0.25_li2s8"]):
    rdf = {}
    for file in rdf_files:
        name = file[0:-4]
        rdf_df = pd.read_csv(f"{rdf_dir}/ps-ps/{sys}/{file}")
        new = rdf_df.iloc[:, 2:]
        rdf[name] = rdf_df[new.columns[new.max().argmax()]].to_list()
    df = pd.DataFrame(rdf)
    df["r ($\AA$)"] = rdf_df["r ($\AA$)"]
    df.columns = [i.upper() for i in df.columns[0:-1]] + ["r ($\AA$)"]
    min_ = np.ceil(df.iloc[:, 0:-1].min().min())
    max_ = round(df.iloc[:, 0:-1].max().max(), 0)
    mid_ = np.mean([min_, max_]).round(0)

    rdf_m = df[(df["r ($\AA$)"] <= 7) & (df["r ($\AA$)"] >= 2)]
    rdf_m = rdf_m.melt("r ($\AA$)", var_name="cols", value_name="g(r)")
    rdf_m = rdf_m.pivot(index="cols", columns="r ($\AA$)", values="g(r)")
    rdf_m = rdf_m.reindex(["DME"] + ["TTE"] + [f"FLS{i}" for i in range(0, 19)])

    xticks, labels = x_axis_ticks_labels(list(rdf_m.columns), 6)

    g = sns.heatmap(
        rdf_m,
        vmin=min_,
        vmax=max_,
        cbar_kws={
            "label": "g(r)",
            "use_gridspec": True,
            "location": "top",
            "shrink": 0.8,
            "ticks": [min_, mid_, max_],
        },
        cmap=cmap,
        ax=axs[idx],
    )
format_rdf_ax(axs, xticks, labels)
save_fig(fig, working_dir="../../files/results", name="rdf_ps_ps1.pdf")


# St (PS) - O (Solvent)
fig, axs = plt.subplots(1, 3, figsize=(30, 15))
for idx, sys in enumerate(["0.25_li2s2", "0.25_li2s4", "0.25_li2s8"]):
    rdf = {}
    for file in rdf_files:
        name = file[0:-4]
        rdf_df = pd.read_csv(f"{rdf_dir}/ps-sol/{sys}/{file}")
        new = rdf_df.iloc[:, 2:]
        rdf[name] = rdf_df[new.columns[new.max().argmax()]].to_list()
    df = pd.DataFrame(rdf)
    df["r ($\AA$)"] = rdf_df["r ($\AA$)"]
    df.columns = [i.upper() for i in df.columns[0:-1]] + ["r ($\AA$)"]
    min_ = np.ceil(df.iloc[:, 0:-1].min().min())
    max_ = round(df.iloc[:, 0:-1].max().max(), 0)
    mid_ = np.mean([min_, max_]).round(0)
    rdf_m = df[(df["r ($\AA$)"] <= 5) & (df["r ($\AA$)"] >= 2)]
    rdf_m = rdf_m.melt("r ($\AA$)", var_name="cols", value_name="g(r)")
    rdf_m = rdf_m.pivot(index="cols", columns="r ($\AA$)", values="g(r)")
    rdf_m = rdf_m.reindex(["DME"] + ["TTE"] + [f"FLS{i}" for i in range(0, 19)])

    xticks, labels = x_axis_ticks_labels(list(rdf_m.columns), 6)

    g = sns.heatmap(
        rdf_m,
        vmin=min_,
        vmax=max_,
        cbar_kws={
            "label": "g(r)",
            "use_gridspec": True,
            "location": "top",
            "shrink": 0.8,
            "ticks": [min_, mid_, max_],
        },
        cmap=cmap,
        ax=axs[idx],
    )
format_rdf_ax(axs, xticks, labels)
save_fig(fig, working_dir="../../files/results", name="rdf_ps_sol.pdf")

# Li - O (Solvent)
fig, axs = plt.subplots(1, 3, figsize=(30, 15))
for idx, sys in enumerate(["0.25_li2s2", "0.25_li2s4", "0.25_li2s8"]):
    rdf = {}
    for file in rdf_files:
        name = file[0:-4]
        rdf_df = pd.read_csv(f"{rdf_dir}/li-sol/{sys}/{file}")
        new = rdf_df.iloc[:, 2:]
        rdf[name] = rdf_df[new.columns[new.max().argmax()]].to_list()
    df = pd.DataFrame(rdf)
    df["r ($\AA$)"] = rdf_df["r ($\AA$)"]
    df.columns = [i.upper() for i in df.columns[0:-1]] + ["r ($\AA$)"]
    min_ = np.ceil(df.iloc[:, 0:-1].min().min())
    max_ = round(df.iloc[:, 0:-1].max().max(), 0)
    mid_ = np.mean([min_, max_]).round(0)
    rdf_m = df[(df["r ($\AA$)"] <= 3) & (df["r ($\AA$)"] >= 0.5)]

    rdf_m = rdf_m.melt("r ($\AA$)", var_name="cols", value_name="g(r)")
    rdf_m = rdf_m.pivot(index="cols", columns="r ($\AA$)", values="g(r)")
    rdf_m = rdf_m.reindex(["DME"] + ["TTE"] + [f"FLS{i}" for i in range(0, 19)])

    xticks, labels = x_axis_ticks_labels(list(rdf_m.columns), 6)

    g = sns.heatmap(
        rdf_m,
        vmin=min_,
        vmax=max_,
        cbar_kws={
            "label": "g(r)",
            "use_gridspec": True,
            "location": "top",
            "shrink": 0.8,
            "ticks": [min_, mid_, max_],
        },
        cmap=cmap,
        ax=axs[idx],
    )
format_rdf_ax(axs, xticks, labels)
save_fig(fig, working_dir="../../files/results", name="rdf_li_sol.pdf")

# Li - O (TFSI)
fig, axs = plt.subplots(1, 3, figsize=(30, 15))
for idx, sys in enumerate(["0.25_li2s2", "0.25_li2s4", "0.25_li2s8"]):
    rdf = {}
    for file in rdf_files:
        name = file[0:-4]
        rdf_df = pd.read_csv(f"{rdf_dir}/li-tfsi/{sys}/{file}")
        new = rdf_df.iloc[:, 2:]
        rdf[name] = rdf_df[new.columns[new.max().argmax()]].to_list()
    df = pd.DataFrame(rdf)
    df["r ($\AA$)"] = rdf_df["r ($\AA$)"]
    df.columns = [i.upper() for i in df.columns[0:-1]] + ["r ($\AA$)"]
    min_ = np.ceil(df.iloc[:, 0:-1].min().min())
    max_ = round(df.iloc[:, 0:-1].max().max(), 0)
    mid_ = np.mean([min_, max_]).round(0)
    rdf_m = df[(df["r ($\AA$)"] <= 3) & (df["r ($\AA$)"] >= 0.5)]

    rdf_m = rdf_m.melt("r ($\AA$)", var_name="cols", value_name="g(r)")
    rdf_m = rdf_m.pivot(index="cols", columns="r ($\AA$)", values="g(r)")
    rdf_m = rdf_m.reindex(["DME"] + ["TTE"] + [f"FLS{i}" for i in range(0, 19)])

    xticks, labels = x_axis_ticks_labels(list(rdf_m.columns), 6)

    g = sns.heatmap(
        rdf_m,
        vmin=min_,
        vmax=max_,
        cbar_kws={
            "label": "g(r)",
            "use_gridspec": True,
            "location": "top",
            "shrink": 0.8,
            "ticks": [min_, mid_, max_],
        },
        cmap=cmap,
        ax=axs[idx],
    )
format_rdf_ax(axs, xticks, labels)
save_fig(fig, working_dir="../../files/results", name="rdf_li_tfsi.pdf")

# Li - O (DOL)
fig, axs = plt.subplots(1, 3, figsize=(30, 15))
for idx, sys in enumerate(["0.25_li2s2", "0.25_li2s4", "0.25_li2s8"]):
    rdf = {}
    for file in rdf_files:
        name = file[0:-4]
        rdf_df = pd.read_csv(f"{rdf_dir}/li-dol/{sys}/{file}")
        new = rdf_df.iloc[:, 2:]
        rdf[name] = rdf_df[new.columns[new.max().argmax()]].to_list()
    df = pd.DataFrame(rdf)
    df["r ($\AA$)"] = rdf_df["r ($\AA$)"]
    df.columns = [i.upper() for i in df.columns[0:-1]] + ["r ($\AA$)"]
    min_ = np.ceil(df.iloc[:, 0:-1].min().min())
    max_ = round(df.iloc[:, 0:-1].max().max(), 0)
    mid_ = np.mean([min_, max_]).round(0)
    rdf_m = df[(df["r ($\AA$)"] <= 3) & (df["r ($\AA$)"] >= 0.5)]

    rdf_m = rdf_m.melt("r ($\AA$)", var_name="cols", value_name="g(r)")
    rdf_m = rdf_m.pivot(index="cols", columns="r ($\AA$)", values="g(r)")
    rdf_m = rdf_m.reindex(["DME"] + ["TTE"] + [f"FLS{i}" for i in range(0, 19)])

    xticks, labels = x_axis_ticks_labels(list(rdf_m.columns), 6)

    g = sns.heatmap(
        rdf_m,
        vmin=min_,
        vmax=max_,
        cbar_kws={
            "label": "g(r)",
            "use_gridspec": True,
            "location": "top",
            "shrink": 0.8,
            "ticks": [min_, mid_, max_],
        },
        cmap=cmap,
        ax=axs[idx],
    )
format_rdf_ax(axs, xticks, labels)
save_fig(fig, working_dir="../../files/results", name="rdf_li_dol.pdf")
