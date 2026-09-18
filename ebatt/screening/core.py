from functools import reduce

import numpy as np
import pandas as pd

from ebatt.config.solvents_config import S1_RULES


def _cond_det(prop, prop_stats, epsilon):
    """
    Defines conditions used in the deterministic screening approach.

    Parameters
    ----------
    prop: pd.DataFrame
        Dataframe of molecules with their properties.
    prop_stats:
        Dataframe of statistics on molecule properties to be included in the
        sampling process. Names of properties should match those in the prop
        df. Should contain the mean and std of each property to be included in
        the sampling process. Does not need to have all properties included in
        the prop df.
    epsilon: float, optional
        Threshold to be used in the det method.

    Returns
    -------
    conditions: list
        List of conditions.
    """
    conditions = []
    for i in prop_stats:
        if S1_RULES[i] == "low":
            cond = (
                prop[i]
                <= prop_stats.loc["mean", i] - (1 - epsilon) * prop_stats.loc["std", i]
            )
        elif S1_RULES[i] == "high":
            cond = (
                prop[i]
                >= prop_stats.loc["mean", i] + (1 - epsilon) * prop_stats.loc["std", i]
            )
        else:
            cond = (
                prop[i]
                >= prop_stats.loc["mean", i] - (1 + epsilon) * prop_stats.loc["std", i]
            ) & (
                prop[i]
                <= prop_stats.loc["mean", i] + (1 + epsilon) * prop_stats.loc["std", i]
            )
        conditions.append(cond)
    return conditions


def _prop_scores(prop, prop_stats):
    """
    Calculates a score for each molecule based on a set of desired properties.

    Parameters
    ----------
    prop: pd.DataFrame
        Dataframe of molecules with their properties.
    prop_stats:
        Dataframe of statistics on molecule properties to be included in the
        sampling process. Names of properties should match those in the prop
        df. Should contain the mean and std of each property to be included in
        the sampling process. Does not need to have all properties included in
        the prop df.

    Returns
    -------
    scores = np.array
        Scores of molecules.
    """
    scores = np.zeros(shape=(prop.shape[0],))
    for i in prop_stats:
        if S1_RULES[i] == "low":
            cond = (
                (prop[i] <= prop_stats.loc["mean", i] - prop_stats.loc["std", i])
                .astype("int")
                .values
            )
        elif S1_RULES[i] == "high":
            cond = (
                (prop[i] >= prop_stats.loc["mean", i] + prop_stats.loc["std", i])
                .astype("int")
                .values
            )
        else:
            cond = (
                (
                    (prop[i] >= prop_stats.loc["mean", i] - prop_stats.loc["std", i])
                    & (prop[i] <= prop_stats.loc["mean", i] + prop_stats.loc["std", i])
                )
                .astype("int")
                .values
            )
        scores += cond
    return scores


def s1_mol(
    prop,
    prop_stats,
    method="det",
    eps=1.5,
    step_size=0.01,
    desired_ratio=0.75,
    lib_size=430,
):
    """
    Samples a library of molecules for those meeting desired criteria. Uses
    either a deterministic or a probabilistic approach. In either approach, the
    majority of the returned dataset includes the desired molecules,
    but will still contain those not meeting the criteria on each property to
    ensure they are still represented.

    Parameters
    ----------
    prop: pd.DataFrame
        Dataframe of molecules with their properties.
    prop_stats:
        Dataframe of statistics on molecule properties to be included in the
        sampling process. Names of properties should match those in the prop
        df. Should contain the mean and std of each property to be included in
        the sampling process. Does not need to have all properties included in
        the prop df.
    method: str, optional
        Method to use for sampling. Supports two methods:
        det: deterministic approach; looks for molecules meeting the
        following criteria; if the desired prop is:
        1. Low: (<= mean_prop - (1 - eps) * std_prop)
        2. Intermediate: (>= mean - (1 + eps) * std) & (<= mean + (1 + eps) * std)
        3. High: (>= mean + (1 - eps) * std)
        Stops when it finds N = desired_ratio * lib_size molecules, and
        finds the rest by randomly sampling from the remaining molecules
        prob: probabilistic approach; samples the library by giving scores to
        each molecule. Molecules meeting most of the desired criteria have
        higher scores.
    eps: float, optional
        Maximum threshold to be used in the det method.
    step_size: float, optional
        Step size to be used in the det method.
    desired_ratio: float, optional
        Ratio of number of molecules with desired properties to total number
        of molecules in the sampled library. Only used in the det approach.
    lib_size: int, optional
        Total number of sampled molecules.

    Returns
    -------
    df: pd.DataFrame
        Dataframe with sampled molecules and their properties.
    summary_stats: pd.DataFrame
        Dataframe with statistics on sampled molecules using properties
        included in the sampling process.
    """
    if method == "det":
        final_eps = 0
        det_df = None
        for e in np.arange(0, eps, step_size):
            conditions = _cond_det(prop, prop_stats, e)
            det_df = prop[reduce(lambda x, y: x & y, conditions)]
            if len(det_df) > desired_ratio * lib_size:
                final_eps = e
                break
        conditions = _cond_det(prop, prop_stats, final_eps)
        df_comp = prop[~reduce(lambda x, y: x & y, conditions)]
        df = pd.concat([det_df, df_comp.sample(lib_size - det_df.shape[0])])

    elif method == "prop":
        scores = _prop_scores(prop, prop_stats)
        print(np.where(scores == np.amax(scores)))
        print(np.where(scores == np.amin(scores)))
        df = prop.sample(lib_size, weights=pd.Series(scores))

    else:
        raise ValueError(
            f"{method} not supported. Choose between det and prop " f"approaches."
        )

    summary_stats = df[prop_stats.columns].describe()
    summary_stats = np.round(summary_stats, decimals=2)
    return df, summary_stats


def s2_dft():
    pass


def s3_md():
    pass
