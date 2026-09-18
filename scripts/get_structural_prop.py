import json

import numpy as np
import pandas as pd

from pymatgen.core.structure import Molecule
from pymatgen.io.gaussian import GaussianInput

from ebatt.descriptors.structural import all_props, sterimol_params, type_connected_carbon


# read cids file
with open("../files/cids.txt", "r") as f:
    cids = json.loads(f.read())

# read pubchem dataframe from file (to avoid parsing pubchem again) and
# compute structural descriptors
df = pd.read_csv("../files/results/solvents.csv")
for ind, cid in enumerate(cids):
    print(f"Processing mol number {ind+1} out of {len(cids)}")
    mol = Molecule.from_file(f"../files/xyz/optimized/{cid}.xyz")
    result = all_props(mol)
    g = GaussianInput.from_file(f"../files/com/{cid}.com")
    ster_param = sterimol_params(g.molecule, f"{cid}.com", f"../files/com")
    type_carbon = type_connected_carbon(mol)
    # keep only the minimum steric effect
    result["steric effect"] = result["steric effect"][2]
    for key, value in result.items():
        df.loc[df["CID"] == cid, key] = value
    for key, value in ster_param.items():
        df.loc[df["CID"] == cid, key] = value
    for key, value in type_carbon[0].items():
        df.loc[df["CID"] == cid, key] = value
df = np.round(df, decimals=2)
df = df.rename(
    columns={
        "radius": "R",
        "gyration radius": "GyrationR",
        "vdW volume": "vdWV",
        "steric effect": "StericEffect",
        "proximity to oxygen": "ProximityO",
        "fluorination degree": "FDeg",
        "shape index": "ShapeInd",
        "# electronegative atoms": "ElectroCount",
    }
)
df.to_csv("../files/results/solvents.csv", index=False)
