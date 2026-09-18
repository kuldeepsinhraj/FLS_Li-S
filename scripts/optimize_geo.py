import os

from fireworks import Workflow, LaunchPad

from mispr.gaussian.workflows.base.core import common_fw
from mispr.gaussian.utilities.fw_utilities import (
    replace_runtask,
    modify_queue_parameters,
)

lpad = LaunchPad.auto_load()

opt_gaussian_inputs = {
    "functional": "wB97XD",
    "basis_set": "Def2TZVPP",
    "route_parameters": {
        "Opt": "(calcfc, tight)",
        "SCF": "(tight, xqc)",
        "int": "ultrafine",
        "NoSymmetry": None,
        "test": None,
        "ExtraBasis": None,
    },
    "link0_parameters": {"%chk": "opt.chk", "%mem": "45GB", "%NProcShared": "24"},
    "input_parameters": {"@/gpfs/home/ratwi/files/def2-tzvppd1gbs": None},
}

freq_gaussian_inputs = {
    "functional": "wB97XD",
    "basis_set": "Def2TZVPP",
    "route_parameters": {
        "Freq": None,
        "iop(7/33=1)": None,
        "NoSymmetry": None,
        "test": None,
        "ExtraBasis": None,
    },
    "link0_parameters": {"%chk": "freq.chk", "%mem": "45GB", "%NProcShared": "24"},
    "input_parameters": {"@/gpfs/home/ratwi/files/def2-tzvppd1gbs": None},
}

working_dir = os.getcwd()
for file in os.listdir(working_dir):
    if file.endswith(".xyz"):
        name = file.strip(".xyz")
        _, _, fws = common_fw(
            mol_operation_type="get_from_file",
            mol=f"{working_dir}/{file}",
            working_dir=working_dir,
            opt_gaussian_inputs=opt_gaussian_inputs,
            freq_gaussian_inputs=freq_gaussian_inputs,
            cart_coords=True,
            oxidation_states={},
            process_mol_func=False,
            format_chk=True,
            mol_name=name,
            tag=f"lis-{name}",
        )
        wf = Workflow(fws, name=name)
        wf = replace_runtask(wf)
        wf = modify_queue_parameters(
            wf,
            walltime="00-48:00:00",
            queue="rajput-24core",
            other_parameters={
                "nodes": 1,
                "ntasks_per_node": 24,
                "ntasks_per_core": 1,
                "ntasks": None,
                "mem": "48GB",
            },
        )
        lpad.add_wf(wf)
