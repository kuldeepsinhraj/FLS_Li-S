import numpy as np


def charge(functional, basis_set, phase, job, tag, db):
    runs = db.retrieve_run(
        functional=functional, basis=basis_set, phase=phase, job_type=job, tag=tag
    )
    oxygen_charges = []
    esp_dict = runs[0]["output"]["output"]["ESP_charges"]
    for key, value in esp_dict.items():
        if value[0] == "O":
            oxygen_charges.append(value[1])
    return np.average(oxygen_charges)


def electronic(functional, basis_set, tag, db, **kwargs):
    doc = db.retrieve_doc(
        collection_name="esp",
        functional=functional,
        basis=basis_set,
        tag=tag,
        **kwargs,
    )
    dipole_moment = doc[0]["dipole_moment"]["Tot"]
    polar = doc[0]["polarizability"] * 0.529177 ** 3  # from Bohr^3 to Angstrom^3
    charges = [i[1] for i in doc[0]["esp"].values()]
    return dipole_moment, polar, charges


def binding_energy(functional, basis_set, tag, db, **kwargs):
    doc = db.retrieve_doc(
        collection_name="binding_energy",
        functional=functional,
        basis=basis_set,
        tag=tag,
        **kwargs,
    )
    be = doc[0]["be_eV"]["value"] / 0.043363  # from eV to kcal/mol
    return be


def redox_potential(functional, basis_set, tag, db, **kwargs):
    doc = db.retrieve_doc(
        collection_name="ip_ea",
        functional=functional,
        basis=basis_set,
        tag=tag,
        **kwargs,
    )
    ip = doc[0]["potentials"]["oxidation"]["solution"]["0e_-1e"]["lithium"]
    ea = doc[0]["potentials"]["reduction"]["solution"]["0e_1e"]["lithium"]
    return ip, ea


def enthalpy():
    pass
