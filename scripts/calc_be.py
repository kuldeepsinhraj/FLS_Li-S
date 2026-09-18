import os

from fireworks import LaunchPad

from pymatgen.core.structure import Molecule

from mispr.gaussian.utilities.fw_utilities import (
    control_worker,
    modify_queue_parameters,
)
from mispr.gaussian.workflows.base.binding_energy import get_binding_energies

lpad = LaunchPad.auto_load()

opt_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "Opt": "(calcfc, tight)",
        "SCF": "(tight, xqc)",
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
        "EmpiricalDispersion": "GD3",
    },
    "link0_parameters": {"%chk": "opt.chk", "%mem": "45GB", "%NProcShared": "28"},
}

freq_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "Freq": None,
        "iop(7/33=1)": None,
        "NoSymmetry": None,
        "test": None,
        "EmpiricalDispersion": "GD3",
    },
    "link0_parameters": {"%chk": "freq.chk", "%mem": "45GB", "%NProcShared": "28"},
}

working_dir = os.getcwd()
for file in os.listdir(working_dir):
    if file.endswith(".xyz"):
        cid = file.strip(".xyz")
        mol = Molecule.from_file(f"{working_dir}/{file}")
        oxygen_sites = []
        for ind, site in enumerate(mol.sites):
            if str(site.specie) == "O":
                oxygen_sites.append(ind)
        wf = get_binding_energies(
            mol_operation_type=["get_from_gout_file", "get_from_mol"],
            mol=[f"{working_dir}/LiTFSI.out", mol],
            index=[15, oxygen_sites[0]],
            process_mol_func=False,
            mol_name=["salt", cid],
            working_dir=f"{working_dir}/{cid}",
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            save_to_db=True,
            save_to_file=True,
            skips=[["opt", "freq"], None],
            additional_prop_doc_fields={"fluorinated-solvent": cid, "ligand": "salt"},
            wall_time=172800,
            tag="lis",
        )
        wf = modify_queue_parameters(
            wf,
            walltime="00-48:00:00",
            queue="rajput-28core",
            other_parameters={
                "nodes": 1,
                "ntasks_per_node": 28,
                "ntasks_per_core": 1,
                "ntasks": None,
                "mem": "48GB",
            },
        )
        wf = control_worker(wf, fworker="Seawulf")
        lpad.add_wf(wf)
