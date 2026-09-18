import os

from fireworks import LaunchPad

from mispr.gaussian.workflows.base.esp import get_esp_charges
from mispr.gaussian.utilities.fw_utilities import (
    control_worker,
    modify_queue_parameters,
)

lpad = LaunchPad.auto_load()

opt_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "Opt": "(calcfc, tight)",
        "SCF": "Tight",
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
    },
    "link0_parameters": {"%chk": "opt.chk", "%mem": "45GB", "%NProcShared": "24"},
}

freq_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "Freq": None,
        "iop(7/33=1)": None,
        "NoSymmetry": None,
        "Polar": None,
        "test": None,
    },
    "link0_parameters": {"%chk": "freq.chk", "%mem": "45GB", "%NProcShared": "24"},
}

esp_gaussian_inputs = {
    "functional": "B3LYP",
    "basis_set": "6-31+G*",
    "route_parameters": {
        "pop": "MK",
        "iop(6/50=1)": None,
        "NoSymmetry": None,
        "Density": None,
        "test": None,
    },
    "link0_parameters": {"%chk": "esp.chk", "%mem": "45GB", "%NProcShared": "24"},
}

working_dir = os.getcwd()
for file in os.listdir(working_dir):
    if file.endswith(".out"):
        cid = file.strip(".out")
        wf = get_esp_charges(
            mol_operation_type="get_from_gout_file",
            mol=file,
            process_mol_func=False,
            mol_name=cid,
            skips=["opt", "freq"],
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            esp_gaussian_inputs=esp_gaussian_inputs,
            format_chk=True,
            save_to_db=True,
            save_to_file=True,
            additional_prop_doc_fields={"fluorinated-solvent": cid},
            tag="lis",
        )
        wf = modify_queue_parameters(
            wf,
            walltime="00-04:00:00",
            queue="rajput-24core",
            other_parameters={
                "nodes": 1,
                "ntasks_per_node": 24,
                "ntasks_per_core": 1,
                "ntasks": None,
                "mem": "48GB",
            },
        )

        wf = control_worker(wf, fworker="Seawulf")
        lpad.add_wf(wf)
