import os
from collections import OrderedDict

import numpy as np
import pandas as pd

from fireworks import LaunchPad
from mispr.lammps.workflows.base import lammps_workflow
from ebatt.descriptors.md import number_molecules
from pymatgen.core.structure import Molecule


lpad = LaunchPad.auto_load()

working_dir = os.getcwd()
df = pd.read_csv(f"{working_dir}/solvents.csv")

# define dictionary of cosolvents to be modeled using MD simulations
cosolvents_dict = {
    "FLS1": 91607507,
    "FLS2": 89393880,
    "FLS3": 139775149,
    "FLS4": 54465191,
    "FLS5": 88850722,
    "FLS6": 21083128,
    "FLS7": 21024143,
    "FLS8": 140985949,
    "FLS9": 150733090,
    "FLS10": 54128414,
    "FLS11": 54212056,
    "FLS12": 54489129,
    "FLS13": 57271144,
    "FLS14": 132186867,
    "FLS15": 59650563,
    "FLS16": 53916352,
    "FLS17": 54183389,
    "FLS18": 18335957,
    "DME": 8071,
    "TTE": 2776662,
    "FLS0": 92551,
}

# define cation force field parameters - Dang
li = Molecule.from_file(f"{working_dir}/pdb/li.pdb")
li.set_charge_and_spin(1, 1)
li_param = {
    "Molecule": li,
    "Labels": ["li"],
    "Masses": OrderedDict({"li": 6.941}),
    "Nonbond": [[0.165, 1.506]],
    "Bonds": [],
    "Angles": [],
    "Dihedrals": [],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray([1]),
}

# define polysulfide force field parameters - Dreiding & Rajput
ps = Molecule.from_file(f"{working_dir}/pdb/s4.pdb")
ps.set_charge_and_spin(-2, 1)
ps_param = {
    "Molecule": ps,
    "Labels": ["sp", "ss", "ss", "sp"],
    "Masses": OrderedDict({"sp": 32.06, "ss": 32.06}),
    "Nonbond": [
        [0.344001, 3.592],
        [0.344001, 3.592],
    ],
    "Bonds": [
        {"coeffs": [157.87, 2.133], "types": [("sp", "ss")]},
        {"coeffs": [154.36, 2.117], "types": [("ss", "ss")]},
    ],
    "Angles": [
        {"coeffs": [73.56, 112.74], "types": [("sp", "ss", "ss")]},
    ],
    "Dihedrals": [
        {
            "coeffs": [
                "fourier",
                4,
                5.4317,
                1,
                0,
                -3.59,
                2,
                180,
                1.0556,
                3,
                0,
                0.05832,
                4,
                180,
            ],
            "types": [("sp", "ss", "ss", "sp")],
        },
    ],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(
        [
            -0.7702,
            -0.2298,
            -0.2298,
            -0.7702,
        ]
    ),
}

# define additive force field parameters - CL&P
no3 = Molecule.from_file(f"{working_dir}/pdb/no3.pdb")
no3.set_charge_and_spin(-1, 1)
no3_param = {
    "Molecule": no3,
    "Labels": ["oa", "oa", "oa", "na"],
    "Masses": OrderedDict({"oa": 15.999400, "na": 14.006700}),
    "Nonbond": [
        [0.1457935, 2.77],
        [0.08078394, 3.06],
    ],
    "Bonds": [
        {"coeffs": [634.2, 1.226], "types": [("oa", "na")]},
    ],
    "Angles": [
        {"coeffs": [120.815, 120.0], "types": [("oa", "na", "oa")]},
    ],
    "Dihedrals": [],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(
        [
            -0.65,
            -0.65,
            -0.65,
            0.95,
        ]
    ),
}

# define anion force field parameters - CL&P
tfsi = Molecule.from_file(f"{working_dir}/pdb/tfsi.pdb")
tfsi.set_charge_and_spin(-1, 1)
tfsi_param = {
    "Molecule": tfsi,
    "Labels": [
        "nt",
        "st",
        "ot",
        "ot",
        "ct",
        "ft",
        "ft",
        "ft",
        "st",
        "ot",
        "ot",
        "ct",
        "ft",
        "ft",
        "ft",
    ],
    "Masses": OrderedDict(
        {
            "nt": 14.006700,
            "st": 32.06,
            "ot": 15.999400,
            "ct": 12.010700,
            "ft": 18.998403,
        }
    ),
    "Nonbond": [
        [0.17, 3.25],
        [0.25, 3.55],
        [0.21, 2.96],
        [0.066, 3.5],
        [0.053, 2.95],
    ],
    "Bonds": [
        {"coeffs": [441.8, 1.323], "types": [("ct", "ft")]},
        {"coeffs": [235.4, 1.818], "types": [("ct", "st")]},
        {"coeffs": [637.1, 1.442], "types": [("st", "ot")]},
        {"coeffs": [372.0, 1.570], "types": [("st", "nt")]},
    ],
    "Angles": [
        {"coeffs": [93.3, 107.1], "types": [("ft", "ct", "ft")]},
        {"coeffs": [82.90, 111.8], "types": [("ft", "ct", "st")]},
        {"coeffs": [104.0, 102.6], "types": [("ct", "st", "ot")]},
        {"coeffs": [115.8, 118.5], "types": [("ot", "st", "ot")]},
        {"coeffs": [94.30, 113.6], "types": [("ot", "st", "nt")]},
        {"coeffs": [97.50, 100.2], "types": [("ct", "st", "nt")]},
        {"coeffs": [80.20, 125.6], "types": [("st", "nt", "st")]},
    ],
    "Dihedrals": [
        {
            "coeffs": ["fourier", 1, 0.1734, 3, 0],
            "types": [("ot", "st", "ct", "ft")],
        },
        {
            "coeffs": ["fourier", 1, 0.158, 3, 0],
            "types": [("nt", "st", "ct", "ft")],
        },
        {
            "coeffs": ["fourier", 1, -0.0018, 3, 0],
            "types": [("ot", "st", "nt", "st")],
        },
        {
            "coeffs": ["fourier", 3, 3.91646, 1, 0, -1.245, 2, 180, -0.3818, 3, 0],
            "types": [("ct", "st", "nt", "st")],
        },
    ],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(
        [
            -0.66,
            1.02,
            -0.53,
            -0.53,
            0.35,
            -0.16,
            -0.16,
            -0.16,
            1.02,
            -0.53,
            -0.53,
            0.35,
            -0.16,
            -0.16,
            -0.16,
        ]
    ),
}

# define DOL force field parameters
dol = Molecule.from_file(f"{working_dir}/pdb/dol.pdb")
dol.set_charge_and_spin(0, 1)
dol_param = {
    "Molecule": dol,
    "Labels": ["os", "os", "ct", "ct", "co", "hc", "hc", "hc", "hc", "hc", "hc"],
    "Masses": OrderedDict(
        {"os": 15.999400, "ct": 12.010700, "co": 12.010700, "hc": 1.007940}
    ),
    "Nonbond": [[0.14, 2.9], [0.066, 3.5], [0.066, 3.5], [0.03, 2.5]],
    "Bonds": [
        {"coeffs": [320.0, 1.41], "types": [("ct", "os")]},
        {"coeffs": [320.0, 1.38], "types": [("co", "os")]},
        {"coeffs": [268.0, 1.529], "types": [("ct", "ct")]},
        {"coeffs": [340.0, 1.09], "types": [("ct", "hc")]},
        {"coeffs": [340.0, 1.09], "types": [("co", "hc")]},
    ],
    "Angles": [
        {"coeffs": [60.0, 109.5], "types": [("co", "os", "ct")]},
        {"coeffs": [50.0, 109.5], "types": [("ct", "ct", "os")]},
        {"coeffs": [35.0, 109.5], "types": [("hc", "ct", "os")]},
        {"coeffs": [37.5, 110.7], "types": [("ct", "ct", "hc")]},
        {"coeffs": [33.0, 107.8], "types": [("hc", "ct", "hc")]},
        {"coeffs": [92.6, 111.55], "types": [("os", "co", "os")]},
        {"coeffs": [35.0, 109.5], "types": [("hc", "co", "os")]},
        {"coeffs": [33.0, 109.5], "types": [("hc", "co", "hc")]},
    ],
    "Dihedrals": [
        {
            "coeffs": ["fourier", 3, 0.0625, 1, 0, -0.011, 2, 180, 0.2805, 3, 0],
            "types": [("ct", "os", "ct", "ct")],
        },
        {
            "coeffs": ["fourier", 3, 0.0625, 1, 0, -0.011, 2, 180, 0.2805, 3, 0],
            "types": [("co", "os", "ct", "ct")],
        },
        {
            "coeffs": ["fourier", 1, 0.3705, 3, 0],
            "types": [("ct", "os", "ct", "hc")],
        },
        {
            "coeffs": ["fourier", 1, 0.3705, 3, 0],
            "types": [("co", "os", "ct", "hc")],
        },
        {
            "coeffs": ["fourier", 3, -0.5335, 1, 0, -0.7835, 2, 180, 0.3375, 3, 0],
            "types": [("ct", "os", "co", "os")],
        },
        {
            "coeffs": ["fourier", 1, 0.3705, 3, 0],
            "types": [("ct", "os", "co", "hc")],
        },
        {
            "coeffs": ["fourier", 3, 1.119, 1, 0, -1.1635, 2, 180, -0.3415, 3, 0],
            "types": [("os", "ct", "ct", "os")],
        },
        {
            "coeffs": ["fourier", 1, 0.234, 3, 0],
            "types": [("hc", "ct", "ct", "os")],
        },
        {"coeffs": ["fourier", 1, 0.15, 3, 0], "types": [("hc", "ct", "ct", "hc")]},
    ],
    "Impropers": [],
    "Improper Topologies": None,
    "Charges": np.asarray(
        [
            -0.4000,
            -0.4000,
            0.1400,
            0.1400,
            0.2000,
            0.0300,
            0.0300,
            0.0300,
            0.0300,
            0.1000,
            0.1000,
        ]
    ),
}

for k, v in cosolvents_dict.items():
    num_mols = number_molecules(
        cosolvent=v, prop_df=df, conc_ps=0.25, cosolvent_ratio=0.5, length=6
    )
    system_species_data = {
        "dol": {
            "molecule": dol,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_dict",
            "ff_param_data": dol_param,
            "mol_mixture_type": "Solvents",
            "mixture_data": num_mols["dol"],
            "save_ff_to_file": True,
        },
        k.lower(): {
            "molecule": f"{working_dir}/pdb/{v}.pdb",
            "molecule_operation_type": "get_from_file",
            "ff_param_method": "get_from_opls",
            "ff_param_data": f"{working_dir}/pdb/{v}.pdb",
            "mol_mixture_type": "Solvents",
            "mixture_data": num_mols[v],
            "save_ff_to_file": True,
        },
        "tfsi": {
            "molecule": tfsi,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_dict",
            "ff_param_data": tfsi_param,
            "mol_mixture_type": "Solutes",
            "mixture_data": num_mols["salt"],
            "save_ff_to_file": True,
        },
        "no3": {
            "molecule": no3,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_dict",
            "ff_param_data": no3_param,
            "mol_mixture_type": "Solutes",
            "mixture_data": num_mols["additive"],
            "save_ff_to_file": True,
        },
        "ps": {
            "molecule": ps,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_dict",
            "ff_param_data": ps_param,
            "mol_mixture_type": "Solutes",
            "mixture_data": num_mols["ps"],
            "save_ff_to_file": True,
        },
        "li": {
            "molecule": li,
            "molecule_operation_type": "get_from_mol",
            "ff_param_method": "get_from_dict",
            "ff_param_data": li_param,
            "mol_mixture_type": "Solutes",
            "mixture_data": num_mols["li"],
            "save_ff_to_file": True,
        },
    }

    wf = lammps_workflow(
        system_species_data=system_species_data,
        system_mixture_type="number of molecules",
        box_data=60.0,
        box_data_type="cubic",
        scale_charges=True,
        charge_scaling_factor=1 / np.sqrt(1.93),
        working_dir=f"{working_dir}/{k.lower()}",
        recipe=[
            ["min", ["template_filename", "min.in"]],
            ["npt", ["template_filename", "npt.in"]],
            ["melt_quench", ["template_filename", "melt_quench.in"]],
            ["nvt", ["template_filename", "nvt.in"]],
        ],
        template_dir=f"{working_dir}/md_templates",
    )
    lpad.add_wf(wf)
