# abbreviations for literature fluorinated solvents
ABBREVIATIONS = [
    "BTFE",
    "DFEP",
    "ETE",
    "ETTFP",
    "FDE",
    "FEP",
    "OTE",
    "PTTFEP",
    "TFEE",
    "TFEP",
    "TFTFE",
    "TTE",
    "TTFEP",
    "FDMB",
]

# radical cids
radical_cids = [
    58537745,
    58537760,
    58700661,
    59344306,
    58544205,
    58544217,
    101084693,
    101938351,
]

# screening rules
S1_RULES = {
    "MolecularWeight": "low",
    "Complexity": "low",
    "R": "low",
    "GyrationR": "low",
    "vdWV": "low",
    "StericEffect": "high",
    "ProximityO": "intermediate",
    "FDeg": "intermediate",
    "ShapeInd": "high",
    "ElectroCount": "low",
    "BP": "high",
    "FP": "high",
    "Solubility_LiTFSI": "high",
    "Solubility_LiFSI": "high",
    "BE_Salt": "low",
    "BE_PS": "high",
    "EW": "high",
    "Polarizability": "low",
    "Dipole": "low"
}

# properties to be used in screening
PROPS = S1_RULES.keys()

# quantitative limits
RULES = {
    "BP": 320,
    "FP": 250,
    "Solubility_LiTFSI": 1,
    "Solubility_LiFSI": 1,
    "EW": 2.7,
}
