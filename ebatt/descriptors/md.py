import numpy as np
import pandas as pd

AVOGADRO = 6.023 * (10**23)
CONVERSION_CONSTANT = 10 ** (-24)

rho = {
    "LiTFSI": 1330,
    "DOL": 1060,
    "DME": 868,
    "Li2S2": 1660,
    "Li2S4": 1660,
    "Li2S6": 1660,
    "Li2S8": 1660,
    "LiNO3": 2380,
}

mw = {
    "LiTFSI": 287.085518,
    "DOL": 74.08,
    "DME": 90.1192,
    "Li2S2": 78.012,
    "Li2S4": 142.182,
    "Li2S6": 206.272,
    "Li2S8": 270.402,
    "LiNO3": 68.946,
}


def number_molecules(cosolvent, prop_df, ps="Li2S4", conc_ps=0.25, cosolvent_ratio=0.5, length=6):
    """
    Calculates the number of species of each type for a solution of 1 M LiTFSI,
    Li2S4, 0.3 M LiNO3 in 1:1 DOL:cosolvent mixture.
    Parameters
    ----------
    cosolvent: int
        CID of the cosolvent.
    prop_df: pd.DataFrame
        Dataframe of solvent properties including CID, MW, and density.
    ps: str
        Type of the polysulfide.
    conc_ps: float
        Concentration of PS to use. Defaults to 0.25 M.
    cosolvent_ratio: float
        volumetric ratio of FLS co-solvent.
    length: int
        Simulation box in nm.

    Returns
    -------
    num_mols: dictionary
        Number of molecules of each type
    """
    cosolvent_prop = prop_df[prop_df["CID"] == cosolvent]
    volume = (length**3) * CONVERSION_CONSTANT  # nm^3 to liters
    number_salt = np.round(1 * AVOGADRO * volume)  # conc_salt = 1
    volume_salt = number_salt * mw["LiTFSI"] / (AVOGADRO * rho["LiTFSI"])
    number_additive = np.round(0.3 * AVOGADRO * volume)  # conc_additive = 0.3
    volume_additive = number_additive * mw["LiNO3"] / (AVOGADRO * rho["LiNO3"])
    number_ps = np.round(conc_ps * AVOGADRO * volume, 0)
    volume_ps = number_ps * mw[ps] / (AVOGADRO * rho[ps])

    remaining_volume = volume - volume_salt - volume_additive - volume_ps

    volume_cosolvent = remaining_volume * cosolvent_ratio
    volume_dol = remaining_volume * (1 - cosolvent_ratio)
    number_dol = np.round(volume_dol * AVOGADRO * rho["DOL"] / mw["DOL"])

    # use more accurate props of DME than the COSMO-RS ones
    if cosolvent == 8071:
        number_cosolvent = np.round(
            volume_cosolvent * AVOGADRO * rho["DME"] / mw["DME"]
        )
    else:
        number_cosolvent = np.round(
            volume_cosolvent
            * AVOGADRO
            * cosolvent_prop["Density"].values[0]
            / cosolvent_prop["MolecularWeight"].values[0]
        )

    number_li = number_ps * 2 + number_salt + number_additive

    num_mols = {
        "li": number_li,
        "salt": number_salt,
        "ps": number_ps,
        "additive": number_additive,
        cosolvent: number_cosolvent,
        "dol": number_dol,
    }
    return num_mols


def gibbs_mixing():
    """
    Based on Flory–Huggins solution theory
    Returns
    -------

    """
    pass
