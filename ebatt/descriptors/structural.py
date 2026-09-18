import os

import collections

from itertools import combinations

import numpy as np
import networkx as nx

from openbabel import pybel as pb

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors

from sterimol import sterimoltools

from mordred import Calculator, descriptors

from pymatgen.io.babel import BabelMolAdaptor
from pymatgen.core.bonds import CovalentBond
from pymatgen.analysis.graphs import MoleculeGraph
from pymatgen.analysis.local_env import OpenBabelNN

VDW_VOLUME = {"C": 20.58, "O": 14.71, "H": 7.24, "N": 15.60, "F": 13.31}


def radius(mol):
    """
    Calculates radius of molecule based on topology. If r_max is the largest
    distance entry in row i of the distance matrix D, then the radius is
    defined as the smallest of r_max. Based on 10.1021/ci00008a012.

    Parameters
    ----------
    mol: Molecule object
        pymatgen Molecule object

    Returns
    -------
    radius: float
        Radius of molecule
    """
    distance = mol.distance_matrix
    max_dist = [max(i) for i in distance]
    return min(max_dist)


def gyration_radius(mol):
    """
    Calculates radius of gyration based on the distribution of atomic masses in
    a molecule.

    Parameters
    ----------
    mol: Molecule object
        pymatgen molecule object.

    Returns
    -------
        radius: float
            Radius of gyration (angstroms).
    """
    com = mol.center_of_mass
    tot_mass = mol.composition.weight
    diff = 0
    for site in mol:
        diff += site.species.weight * (np.linalg.norm(site.coords - com)) ** 2
    return np.sqrt(diff / tot_mass)


def vdw_volume(mol):
    """
    Calculates the van der Waals volume of a molecule, as the space occupied by
    the molecule and not penetrable by other molecules. Uses atomic volume
    contributions and number of bonds. Assumes no aromatic and non-aromatic
    rings in the molecules. The true volume will be slightly different since
    the atomic contribution is not constant but depends on the nature of the
    surrounding atoms. Based on 10.1021/jo034808o.

    Parameters
    ----------
    mol: Molecule object
        pymatgen Molecule object.

    Returns
    -------
    volume: float
        van der Waals volume of molecule in angstroms^3.
    """
    mol_graph = MoleculeGraph.with_local_env_strategy(mol, OpenBabelNN())
    num_bonds = len(mol_graph.graph.edges)
    mol_comp = mol.composition.as_dict()
    atomic_cont = sum(mol_comp[key] * VDW_VOLUME[key] for key in mol_comp)
    volume = atomic_cont - 5.92 * num_bonds
    return volume


def shape_index(mol):
    """
    Calculates the shape index of a molecule (spherical < 0.5 < linear). First,
    finds the shortest path between each non-hydrogen atom pair in the
    molecule. Then, returns the shape index as the longest of these
    paths/number of non-hydrogen atoms. Based on
    https://openmolecules.org/forum/index.php?t=msg&goto=197&.

    Parameters
    ----------
    mol: Molecule object
        pymatgen molecule object.

    Returns
    -------
    shape_ind: float
        Shape index of the molecule.
    """
    removed_h = mol.copy()
    removed_h.remove_species("H")
    mol_graph = MoleculeGraph.with_local_env_strategy(removed_h, OpenBabelNN())
    mol_graph_und = mol_graph.graph.to_undirected()
    shortest_dist = []
    node_pairs = [pair for pair in combinations(mol_graph_und.nodes, 2)]
    for pair in node_pairs:
        num_nodes = []
        for path in nx.all_simple_paths(mol_graph_und, source=pair[0], target=pair[1]):
            num_nodes.append(len(path))
        # there can be multiple possible paths connecting 2 nodes (e.g. cyclic),
        # so take shortest
        shortest_dist.append(min(num_nodes))
    shape_ind = max(shortest_dist) / (mol.num_sites - mol.composition.as_dict()["H"])
    return shape_ind


def steric_effect(mol):
    """
    Calculates steric effect of an oxygen atom in a fluorinated ether
    solvent using covalent radii values and geometrical distances. Based on
    10.1021/ci034266b (Eq. 8b). This is only an approximation and not an exact
    value of steric strains.

    Parameters
    ----------
    mol: Molecule object
        pymatgen molecule object.

    Returns
    -------
    effect: dictionary
        Steric effect at each oxygen site.

    avg_effect: float
        Average steric effect at oxygen sites in the molecule.
    """
    mol_graph = MoleculeGraph.with_local_env_strategy(mol, OpenBabelNN())
    mol_graph_und = mol_graph.graph.to_undirected()
    oxygen_sites = []
    for ind, site in enumerate(mol.sites):
        if str(site.specie) == "O":
            oxygen_sites.append(ind)

    node_pairs = [pair for pair in combinations(mol_graph_und.nodes, 2)]

    # get pairs involving an oxygen site only
    oxygen_node_pairs = []
    for oxygen in oxygen_sites:
        oxygen_pairs = []
        for pair in node_pairs:
            if oxygen in pair:
                # put oxygen site index first
                oxygen_pairs.append((oxygen, [i for i in pair if i != oxygen][0]))
        oxygen_node_pairs.append(oxygen_pairs)

    effect = {}
    for oxygen_site, oxygen in zip(oxygen_sites, oxygen_node_pairs):
        distance = 0
        for pair in oxygen:
            site_radius = mol[pair[1]].specie.atomic_radius
            site_radius = site_radius / 0.7
            num_nodes = {}
            for path in nx.all_simple_paths(
                mol_graph_und, source=pair[0], target=pair[1]
            ):
                num_nodes[len(path)] = path
            shortest_dist = min(list(num_nodes.keys()))
            shortest_path = num_nodes[shortest_dist]
            dist = 0
            for n1, n2 in zip(shortest_path, shortest_path[1:]):
                dist += mol.get_distance(n1, n2)
            dist = dist / (0.7 * 2)
            distance += site_radius**3 / dist**3
        effect[oxygen_site] = distance
    avg_effect = np.average(list(effect.values()))
    min_effect = np.min(list(effect.values()))
    return effect, avg_effect, min_effect


def sterimol_params(mol, g_file, path):
    """
    Calculates Verloop's Sterimol parameters (10.1021/acscatal.8b04043) of a molecule
    along the direction of the coordinated lithium cation and the oxygen atom of the
    solvent molecule using DFT optimized structures.

    Parameters
    ----------
    mol: Molecule object
        pymatgen molecule object.
    g_file: str
        Gaussian input file name.
    path: str
        Path to the Gaussian input file.

    Returns
    -------
    ster_param: dictionary
        Sterimol parameters of a molecule.
    """
    working_dir = os.getcwd()
    oxygen_sites = []
    for ind, site in enumerate(mol.sites):
        if str(site.specie) == "O":
            oxygen_sites.append(ind)

    os.chdir(path)
    ster_param = sterimoltools.calcSterimol(
        g_file, "bondi", oxygen_sites[0] + 1, 1, True
    )
    ster_param = {"L": ster_param.lval, "B1": ster_param.B1, "B5": ster_param.newB5}
    os.chdir(working_dir)
    return ster_param


def proximity_to_oxygen(mol):
    """
    Calculates the proximity of fluorine atoms to oxygen atoms in a fluorinated
    solvent molecule. Correlated to the steric hindrance effect which can
    affect the solubility of polysulfides in a Li-S electrolyte. Assumes the
    presence of at least one oxygen and one fluorine atom. If more than one
    pair exists, returns the minimum distance. Uses direct distances.

    Parameters
    ----------
    mol: Molecule object
        pymatgen molecule object.

    Returns
    -------
    distance: float
        Minimum distance between a fluorine and an oxygen atom in a fluorinated
        solvent molecule (angstroms).
    """
    # get oxygen and fluorine site indexes
    oxygen_sites = []
    fluorine_sites = []
    for ind, site in enumerate(mol.sites):
        if str(site.specie) == "O":
            oxygen_sites.append(ind)
        elif str(site.specie) == "F":
            fluorine_sites.append(ind)
    # calculate distance between each fluorine atom and each oxygen atom
    distances = []
    for fluorine_site in fluorine_sites:
        for oxygen_site in oxygen_sites:
            distances.append(mol.get_distance(fluorine_site, oxygen_site))
    return np.min(distances) if distances else np.nan


def fluorination_degree(mol):
    """
    Calculates the degree of fluorination of a solvent molecule as a %.

    Parameters
    ----------
    mol: Molecule object
        pymatgen molecule object.

    Returns
    -------
    Fluorination degree: float
        Degree of fluorination of a solvent molecule.

    """
    comp_dict = mol.composition.as_dict()
    return comp_dict["F"] * 100 / (comp_dict["F"] + comp_dict["H"])


def perc_fluorinated_carbon(mol):
    """
    Calculates the % of Carbon atoms that are fluorinated.

    Parameters
    ----------
    mol: Molecule object
        pymatgen molecule object.

    Returns
    -------
    percent: float
        Percentage of fluorinated carbon atoms.

    """

    def _get_percent_non_empty(dictionary):
        non_empty_count = 0
        for key, value in dictionary.items():
            if value:  # Check if the value is non-empty
                non_empty_count += 1
        return non_empty_count / len(dictionary) * 100

    carbon_sites = []
    fluorine_sites = []
    for ind, site in enumerate(mol.sites):
        if str(site.specie) == "C":
            carbon_sites.append(ind)
        elif str(site.specie) == "F":
            fluorine_sites.append(ind)

    bonds = {}
    for carbon_site in carbon_sites:
        bonds[carbon_site] = []
        for fluorine_site in fluorine_sites:
            if CovalentBond.is_bonded(mol[carbon_site], mol[fluorine_site]):
                bonds[carbon_site].append(fluorine_site)
    return _get_percent_non_empty(bonds)


def num_ether_oxygen(mol):
    """
    Calculates the number of ether oxygens in the fluorinated solvent.

    Parameters
    ----------
    mol: Molecule object
        pymatgen molecule object.

    Returns
    -------
    num: int
        Number of oxygen atoms.

    """

    oxygen_sites = []
    for ind, site in enumerate(mol.sites):
        if str(site.specie) == "O":
            oxygen_sites.append(ind)
    return len(oxygen_sites)


def type_connected_carbon(mol):
    def convert_list(lst):
        result = []
        for d in lst:
            s = ""
            for k, v in d.items():
                s += k + str(v)
            result.append(s)
        return result

    oxygen_sites = []
    connections = []
    types_carbon = ["F3", "F2", "F2H1", "F1H1", "H3", "H2", "F1H2", "H1"]
    for ind, site in enumerate(mol.sites):
        if str(site.specie) == "O":
            oxygen_sites.append(ind)
    connected_atom_indices = []
    for oxygen_site in oxygen_sites:
        connected_carbons = mol.get_neighbors(mol[oxygen_site], 1.5)
        if connected_carbons:
            for carbon in connected_carbons:
                if str(carbon.specie) == "C":
                    connected_atoms = mol.get_neighbors(carbon, 1.5)
                    connected_atom_indices.append(
                        [i.index + 1 for i in connected_atoms if str(i.specie) != "O"]
                        + [carbon.index + 1]
                    )
                    connected_atoms = [
                        str(i.specie) for i in connected_atoms if str(i.specie) != "O"
                    ]
                    count_dict = dict(collections.Counter(connected_atoms))
                    connections.append({k: v for k, v in sorted(count_dict.items())})
    connections = convert_list(connections)
    result = {}

    for type_carbon in types_carbon:
        result[type_carbon] = 0

    for connection in connections:
        if connection in result:
            result[connection] += 1
    return result, connected_atom_indices


def electronegative_atoms(mol):
    """
    Counts number of electronegative atoms in a molecule.

    Parameters
    ----------
    mol: Moleucle object
        pymatgen molecule object.

    Returns
    -------
    nm_elec: int
        Number of electronegative atoms.
    """
    electo_atoms = ["N", "O", "P", "S", "F", "Cl", "Br", "I", "As", "Se"]
    mol_comp = mol.composition.as_dict()
    num_elec = int(sum([mol_comp[i] for i in electo_atoms]))
    return num_elec


def similarity_measure(pattern_mol, target_mol, fptype="ECFP4"):
    """
    Calculates the similarity measure between two molecules using the Tanimoto
    index as implemented in openbabel. Defined as number of bits set in
    (patternFP & targetFP) / Number of bits in (patternFP | targetFP).

    Parameters
    ----------
    fptype: str
        Type of fingerprint. Supported onces in open babel are:
            ECFP0
            ECFP10
            ECFP2
            ECFP4
            ECFP6
            ECFP8
            FP2
            FP3
            FP4
            MACCS
    pattern_mol: Molecule object
        pymatgen Molecule object.
    target_mol: Molecule object
        pymatgen Molecule object.

    Returns
    -------
    similarity measure: float
        Level of similarity between two molecules (value between 0 and 1).
    """
    principle_mol_pb = pb.Molecule(BabelMolAdaptor(pattern_mol).openbabel_mol)
    other_mol_pb = pb.Molecule(BabelMolAdaptor(target_mol).openbabel_mol)
    mols = [principle_mol_pb, other_mol_pb]
    fps = [x.calcfp(fptype) for x in mols]
    return fps[0] | fps[1]


def get_Mordred_descriptors(data):
    calc = Calculator(descriptors, ignore_3D=False)
    mols = [Chem.MolFromSmiles(smi) for smi in data]
    df = calc.pandas(mols)
    return df


def rdkit_descriptors(smiles):
    mols = [Chem.MolFromSmiles(i) for i in smiles]
    calc = MoleculeDescriptors.MolecularDescriptorCalculator(
        [x[0] for x in Descriptors._descList]
    )
    desc_names = calc.GetDescriptorNames()

    mol_descriptors = []
    for mol in mols:
        # add hydrogens to molecules
        mol = Chem.AddHs(mol)
        # Calculate all 200 descriptors for each molecule
        descriptors = calc.CalcDescriptors(mol)
        mol_descriptors.append(descriptors)
    return mol_descriptors, desc_names


def all_props(mol):
    """
    Get the dictionary of all structural descriptors of a given molecule.

    Parameters
    ----------
    mol: Molecule object
        pymatgen Molecule object.

    Returns
    -------
    result: dictionary
        Eight calculated properties of a given molecule:
            radius
            gyration radius
            vdW volume
            shape index
            steric hindrance
            proximity to oxygen
            fluorination degree
            number of electronegative atoms

    """
    result = {}
    properties = {
        "radius": radius,
        "gyration radius": gyration_radius,
        "vdW volume": vdw_volume,
        "shape index": shape_index,
        "steric effect": steric_effect,
        "proximity to oxygen": proximity_to_oxygen,
        "fluorination degree": fluorination_degree,
        "# electronegative atoms": electronegative_atoms,
    }

    for prop in properties.keys():
        result[prop] = properties[prop](mol)
    return result
