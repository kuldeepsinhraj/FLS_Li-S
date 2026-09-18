import os

from fireworks import LaunchPad

from pymatgen.core.structure import Molecule

from mispr.gaussian.workflows.base.ip_ea import get_ip_ea
from mispr.gaussian.utilities.fw_utilities import (
    control_worker,
    modify_queue_parameters,
)

lpad = LaunchPad.auto_load()

opt_gaussian_inputs = {
    "functional": "wB97XD",
    "basis_set": "Def2TZVP",
    "route_parameters": {
        "Opt": "(calcfc, tight)",
        "SCF": "(tight, xqc)",
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
        "Extrabasis": None,
    },
    "link0_parameters": {"%chk": "opt.chk", "%mem": "45GB", "%NProcShared": "28"},
    "input_parameters": {"@/gpfs/home/ratwi/files/def2-tzvppd1gbs": None},
}

freq_gaussian_inputs = {
    "functional": "wB97XD",
    "basis_set": "Def2TZVP",
    "route_parameters": {
        "Freq": None,
        "iop(7/33=1)": None,
        "NoSymmetry": None,
        "test": None,
        "Extrabasis": None,
    },
    "link0_parameters": {"%chk": "freq.chk", "%mem": "45GB", "%NProcShared": "28"},
    "input_parameters": {"@/gpfs/home/ratwi/files/def2-tzvppd1gbs": None},
}

working_dir = os.getcwd()
for file in os.listdir(working_dir):
    if file.endswith(".xyz"):
        cid = file.strip(".xyz")
        mol = Molecule.from_file(f"{working_dir}/{file}")
        wf = get_ip_ea(
            mol_operation_type="get_from_file",
            mol=file,
            ref_charge=0,
            process_mol_func=False,
            mol_name=cid,
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            solvent_gaussian_inputs=f"(SMD, Solvent=Acetone)",
            num_electrons=1,
            vertical=False,
            states=["cation", "anion"],
            phases=["gas", "solution"],
            gibbs_elec=0,
            job_type="better_guess",
            lower_functional="B3LYP",
            lower_basis_set="6-31+G*",
            save_to_db=True,
            save_to_file=True,
            additional_prop_doc_fields={"fluorinated-solvent": cid},
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
