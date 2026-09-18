import os

import pylab

import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from mpl_toolkits.mplot3d import art3d
from matplotlib.ticker import AutoMinorLocator
from matplotlib.colors import LinearSegmentedColormap

sns.set_palette("colorblind")


def _set_axis(ax, axis="both", labelsize=18):
    """
    Formats the axes on a plot.

    Parameters
    ----------
    ax: matplotlib.axes._subplots.AxesSubplot
    axis: str
        Axis to format: x, y, or both

    Returns
    -------
    None
    """
    if axis == "both":
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: "{:g}".format(x))
        )
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda y, _: "{:g}".format(y))
        )
        ax.tick_params(which="major", length=10, pad=10)
        ax.tick_params(which="minor", length=5, pad=10)
        ax.tick_params(
            axis="both", which="both", direction="in", labelsize=labelsize, pad=10
        )
    elif axis == "x":
        ax.xaxis.set_minor_locator(AutoMinorLocator(2))
        ax.xaxis.set_major_formatter(
            ticker.FuncFormatter(lambda x, _: "{:g}".format(x))
        )
        ax.tick_params(which="major", length=10, pad=10)
        ax.tick_params(which="minor", length=5, pad=10)
        ax.tick_params(
            axis="x", which="both", direction="in", labelsize=labelsize, pad=10
        )
    elif axis == "y":
        ax.yaxis.set_minor_locator(AutoMinorLocator(2))
        ax.yaxis.set_major_formatter(
            ticker.FuncFormatter(lambda y, _: "{:g}".format(y))
        )
        ax.tick_params(which="major", length=10, pad=10)
        ax.tick_params(which="minor", length=5, pad=10)
        ax.tick_params(
            axis="y", which="both", direction="in", labelsize=labelsize, pad=10
        )


def lims(mplotlims):
    scale = 1.021
    offset = (mplotlims[1] - mplotlims[0]) * scale
    return mplotlims[1] - offset, mplotlims[0] + offset


def format_3d_axes(ax):
    xlims, ylims, zlims = lims(ax.get_xlim()), lims(ax.get_ylim()), lims(ax.get_zlim())
    axes = [
        [(0, 1), (1, 1), (0, 0)],
        [(0, 1), (0, 0), (1, 1)],
        [(1, 1), (1, 1), (0, 1)],
        [(1, 1), (0, 1), (1, 1)],
        [(0, 0), (0, 0), (0, 1)],
        [(0, 0), (0, 1), (0, 0)],
    ]
    for x, y, z in axes:
        i = pylab.array([xlims[x[0]], ylims[y[0]], zlims[z[0]]])
        f = pylab.array([xlims[x[1]], ylims[y[1]], zlims[z[1]]])
        p = art3d.Poly3DCollection(pylab.array([[i, f]]))
        p.set_color("black")
        p.set_linewidth(2.5)
        ax.add_collection3d(p)

    ax.xaxis.pane.fill = False
    ax.yaxis.pane.fill = False
    ax.zaxis.pane.fill = False

    ax.xaxis._axinfo["grid"].update(
        {"linewidth": 0.5, "linestyle": "--", "color": "lightgray"}
    )
    ax.yaxis._axinfo["grid"].update(
        {"linewidth": 0.5, "linestyle": "--", "color": "lightgray"}
    )
    ax.zaxis._axinfo["grid"].update(
        {"linewidth": 0.5, "linestyle": "--", "color": "lightgray"}
    )

    ax.tick_params(which="major", direction="in", labelsize=18, length=10, pad=1)
    ax.set_box_aspect((2, 3.5, 2))
    # ax.view_init(azim=-165, elev=5)
    ax.view_init(azim=136, elev=19)


def custom_cmap(
    colors,
    numcolors=256,
):
    cmap = LinearSegmentedColormap.from_list(
        name="custom_cmap",
        colors=colors,
        N=numcolors,
    )
    return cmap


def create_fig(
    xlabel=None,
    ylabel=None,
    title=None,
    figsize=(8, 6),
    set_axis=True,
    axis="both",
    fontsize=18,
    labelsize=18,
    labelpad=10,
):
    fig, ax = plt.subplots(figsize=figsize)
    if set_axis:
        _set_axis(ax, axis=axis, labelsize=labelsize)
    if xlabel:
        ax.set_xlabel(xlabel, fontsize=fontsize, weight="bold", labelpad=labelpad)
    if ylabel:
        ax.set_ylabel(ylabel, fontsize=fontsize, weight="bold", labelpad=labelpad)
    if title:
        ax.set_title(title, fontsize=fontsize, weight="bold")
    return fig, ax


def dist_plot(ax, data, bins, label, color=None, alpha=0.5):
    ax = sns.distplot(
        data,
        hist=True,
        kde=True,
        bins=bins,
        color=color,
        hist_kws={"edgecolor": "black", "alpha": alpha},
        kde_kws={"linewidth": 2},
        label=label,
        ax=ax,
    )
    return ax


def ecdf_plot(ax, data, color):
    ax = sns.ecdfplot(data, ax=ax, color=color, linewidth=4)
    return ax


def scatter_plot(ax, x, y, color, size=100, marker="x", label=None):
    ax.scatter(x, y, color=color, s=size, marker=marker, edgecolor="k", label=label)


def line_plot(
    ax,
    x,
    y,
    fmt,
    markerfacecolor,
    markeredgecolor,
    color,
    markersize=10,
    linewidth=2.5,
    label=None,
    yerr=None,
    xerr=None,
    show_errorbars=False,
):
    if show_errorbars:
        ax.errorbar(
            x,
            y,
            yerr=yerr,
            xerr=xerr,
            fmt=fmt,
            markerfacecolor=markerfacecolor,
            markeredgecolor=markeredgecolor,
            color=color,
            markersize=markersize,
            label=label,
            linewidth=linewidth,
        )
    else:
        ax.plot(
            x,
            y,
            fmt,
            markerfacecolor=markerfacecolor,
            markeredgecolor=markeredgecolor,
            color=color,
            markersize=markersize,
            label=label,
            linewidth=linewidth,
        )

    return ax


def edit_plot(ax, legend=None, loc="best", linewidth=2):
    if legend:
        ax.legend(legend, frameon=False, fontsize=18, loc=loc)
    for axis in ["top", "bottom", "left", "right"]:
        ax.spines[axis].set_linewidth(linewidth)
    return ax


def save_fig(fig, working_dir=None, name="fig.png"):
    if not working_dir:
        working_dir = os.getcwd()
    if not os.path.isdir(working_dir):
        os.mkdir(working_dir)
    fig.tight_layout(pad=5)
    fig.savefig(f"{working_dir}/{name}", bbox_inches="tight", pad_inches=0.1)
