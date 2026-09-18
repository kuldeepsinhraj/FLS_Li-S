import json

import pandas as pd

from ebatt.data.pubchem import get_pubchem_computed_prop

# read cids file
with open("../files/cids.txt", "r") as f:
    cids = json.loads(f.read())

# parse properties from pubchem
pubchem_props = []
for cid in cids:
    output = get_pubchem_computed_prop(
        cid,
        [
            "CanonicalSMILES",
            "Complexity",
            "IUPACName",
            "InChI",
            "MolecularFormula",
            "MolecularWeight",
            "RotatableBondCount",
            "Volume3D",
            "XStericQuadrupole3D",
            "YStericQuadrupole3D",
            "ZStericQuadrupole3D",
        ],
    )
    output["2D"] = f'<img src="../files/png/{cid}.png" />'
    pubchem_props.append(output)

df = pd.DataFrame(pubchem_props)
df.to_csv("../files/results/solvents.csv", index=False)
