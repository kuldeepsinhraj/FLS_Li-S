import json

import numpy as np
import pandas as pd

from IPython.display import HTML

from ebatt.data.molmd import electronic, binding_energy, redox_potential
from ebatt.config.db_config import DB_FILE
from ebatt.config.solvents_config import PROPS

from mispr.gaussian.database import GaussianCalcDb

# provide database credentials
db = GaussianCalcDb.from_db_file(DB_FILE)

# read dataframe containing pubchem, structural, and cosmors properties
df = pd.read_csv("../files/results/solvents.csv")

elec = {
    "CID": [],
    "Dipole": [],
    "Polarizability": [],
    "MinPartialCharge": [],
    "MaxPartialCharge": [],
    "MinAbsPartialCharge": [],
    "MaxAbsPartialCharge": [],
}
for cid in df["CID"]:
    try:
        dipole, polar, charges = electronic(
            functional="B3LYP",
            basis_set="6-31+G*",
            tag="lis",
            db=db,
            **{"fluorinated-solvent": str(cid)},
        )
        elec["CID"].append(cid)
        elec["Dipole"].append(dipole)
        elec["Polarizability"].append(polar)
        elec["MinPartialCharge"].append(min(charges))
        elec["MaxPartialCharge"].append(max(charges))
        elec["MinAbsPartialCharge"].append(min([abs(i) for i in charges]))
        elec["MaxAbsPartialCharge"].append(max([abs(i) for i in charges]))
    except Exception as e:
        print("esp calculation not complete")
elec_df = pd.DataFrame.from_dict(elec)
df = df.merge(elec_df, on="CID", how="left")

be = {"CID": [], "BE_Salt": []}
for cid in df["CID"]:
    try:
        value = binding_energy(
            functional="B3LYP",
            basis_set="6-31+G*",
            tag="lis",
            db=db,
            **{"fluorinated-solvent": str(cid), "ligand": "salt"},
        )
        be["CID"].append(cid)
        be["BE_Salt"].append(value)
    except Exception as e:
        print("binding energy calculation not complete")
be_df = pd.DataFrame.from_dict(be)
df = df.merge(be_df, on="CID", how="left")

be = {"CID": [], "BE_PS": []}
for cid in df["CID"]:
    try:
        value = binding_energy(
            functional="B3LYP",
            basis_set="6-31+G*",
            tag="lis",
            db=db,
            **{"fluorinated-solvent": str(cid), "ligand": "ps"},
        )
        be["CID"].append(cid)
        be["BE_PS"].append(value)
    except Exception as e:
        print("binding energy calculation not complete")
be_df = pd.DataFrame.from_dict(be)
df = df.merge(be_df, on="CID", how="left")

redox = {"CID": [], "IP": [], "EA": []}
for cid in df["CID"]:
    try:
        oxd, red = redox_potential(
            functional="wB97XD",
            basis_set="Def2TZVP",
            tag="lis",
            db=db,
            **{"fluorinated-solvent": str(cid)},
        )
        redox["CID"].append(cid)
        redox["IP"].append(oxd)
        redox["EA"].append(red)
    except Exception as e:
        print("redox calculation not complete")
redox_df = pd.DataFrame.from_dict(redox)
df = df.merge(redox_df, on="CID", how="left")
df["EW"] = df["IP"] - df["EA"]

df.to_csv("../files/results/solvents.csv", index=False)
df_html = df.to_html(open("../files/results/solvents.html", "w"), escape=False)
HTML(df.to_html(open("../files/results/solvents.html", "w"), escape=False))

# get summary statistics of properties
df_prop = df[PROPS]
summary = df_prop.describe()
summary = np.round(summary, decimals=2)
summary.to_csv("../files/results/stats.csv")
