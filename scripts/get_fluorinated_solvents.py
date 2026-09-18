import json

from ebatt.data.pubchem import (
    get_pubchem_cid,
    download_pubchem_file,
    get_pubchem_similar_mol,
    filter_pubchem_molecules,
)

# read inchis of fluorinated solvents used in li-s battery lit.
with open("../files/lit_inchis.txt") as f:
    inchis = f.read().splitlines()

# get unique molecules sharing similar structural features to solvents in lit.
cids = []
lit_cids = []
for num, inchi in enumerate(inchis):
    cid = get_pubchem_cid(inchi)
    cids += [cid]
    lit_cids += [cid]
    cids += get_pubchem_similar_mol(cid)
flat_cids = list(set(cids))
print(f"number of similar training cids: {len(flat_cids)}")

# filter for neutral, single-molecule, non-aromatic, and C-O-H-F containing mols
cids = filter_pubchem_molecules(
    flat_cids,
    quant_criteria={
        "Charge": [0, 0],
        "CovalentUnitCount": [1, 1],
        "FeatureRingCount3D": [0, 0],
        "DefinedAtomStereoCount": [0, 0],
        "IsotopeAtomCount": [0, 0],
    },
    atom_types={"C", "O", "H", "F"},
)
lit_cids = [lit_cid for lit_cid in lit_cids if lit_cid in cids]

# add DME and DOL cids to cids
cids += [8071, 12586]

# write cids to file
with open("../files/lit_cids.txt", "w") as f:
    f.write(json.dumps(lit_cids))
with open("../files/cids.txt", "w") as f:
    f.write(json.dumps(cids))

# download sdf and png files
for cid in cids:
    download_pubchem_file(cid, fmt="sdf", working_dir="../files")
    download_pubchem_file(cid, fmt="png", working_dir="../files")
