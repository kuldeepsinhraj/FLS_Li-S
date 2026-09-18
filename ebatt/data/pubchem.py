import os
import ssl
import json
import urllib.request

import pubchempy as pcp

ssl._create_default_https_context = ssl._create_unverified_context


def get_pubchem_cid(inchi):
    """
    Uses pubchempy to get the cid of a compound using its inchi representation.
    Returns 0 if molecule is not found.

    Parameters
    ----------
    inchi : str
        The InChI representation of the molecule.

    Returns
    -------
    cid : int
        PubChem Compound ID.
    """
    comps = pcp.get_compounds("InChI={}".format(inchi), "inchi")
    if comps:
        return comps[0].cid
    else:
        return 0


def get_pubchem_similar_mol(cid):
    """
    Uses the CID of a molecule to find similar structures in the PubChem
    database.

    Parameters
    ----------
    cid: int
        PubChem Compound ID.

    Returns
    -------
    list_similar_cids: list of int
        List of CIDs of molecules with similar structural features to the input.
    """
    pubchem_pug_rest_api_link = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/"
    pubchem_pug_rest_api_link += "compound/fastsimilarity_2d/cid/%s/cids/txt" % cid
    results = urllib.request.urlopen(pubchem_pug_rest_api_link).read()
    results = results.decode("utf-8").split("\n")
    list_similar_cids = [int(cid) for cid in results if cid]
    return list_similar_cids


def download_pubchem_file(cid, fmt="sdf", working_dir=None):
    """
    Downloads files from PubChem database using a molecule cid.

    Parameters
    ----------
    cid: int
        PubChem Compound CID.
    fmt: int, optional
         Format of file to download. Can be  XML, ASNT/B, JSON, SDF, CSV, PNG,
         TXT. Default is SDF.
    working_dir: str, optional
        Path where to create a folder for the specified format and save the
        file. Default is current directory.

    Returns
    -------
    None
    """
    if not working_dir:
        working_dir = os.getcwd()
    full_dir = f"{working_dir}/{fmt.lower()}"
    if not os.path.isdir(full_dir):
        os.makedirs(full_dir)
    try:
        pcp.download(fmt, f"{full_dir}/{cid}.{fmt.lower()}", cid, "cid")
    except Exception as e:
        print(e)


def parse_pubchem(url):
    """
    Takes PubChem API link, parses it, and converts it to a dictionary.

    Parameters
    ----------
    url: str
        URL link to use for parsing.

    Returns
    -------
    result_dict: dict
        Dictionary of parsed PubChem output.
    """
    request = urllib.request.Request(url)
    result = urllib.request.urlopen(request).read()
    result_dict = json.loads(result.decode("latin-1"))
    return result_dict


def get_pubchem_computed_prop(cid, properties):
    """
    Parses PubChem database using molecule cid and outputs the requested
    properties (works for selected properties only) in a dictionary format.

    Parameters
    ----------
    cid: int
        PubChem Compound CID.
    properties: list of str
        Properties to parse from PubChem database. Supported properties include:
            AtomStereoCount
            BondStereoCount
            CanonicalSMILES
            Charge
            Complexity
            ConformerCount3D
            ConformerModelRMSD3D
            CovalentUnitCount
            DefinedAtomStereoCount
            DefinedBondStereoCount
            EffectiveRotorCount3D
            ExactMass
            FeatureAcceptorCount3D
            FeatureAnionCount3D
            FeatureCationCount3D
            FeatureCount3D
            FeatureDonorCount3D
            FeatureHydrophobeCount3D
            FeatureRingCount3D
            Fingerprint2D
            HBondAcceptorCount
            HBondDonorCount
            HeavyAtomCount
            IUPACName
            InChI
            InChIKey
            IsomericSMILES
            IsotopeAtomCount
            MolecularFormula
            MolecularWeight
            MonoisotopicMass
            RotatableBondCount
            TPSA
            UndefinedAtomStereoCount
            UndefinedBondStereoCount
            Volume3D
            XLogP
            XStericQuadrupole3D
            YStericQuadrupole3D
            ZStericQuadrupole3D

    Returns
    -------
    computed_properties: dict
        Dictionary of parsed properties.
    """
    requested_properties = ",".join(properties)
    rest_api_link = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/"
        f"property/{requested_properties}/JSON"
    )
    computed_properties = parse_pubchem(rest_api_link)["PropertyTable"]["Properties"][0]
    return computed_properties


def get_pubchem_other_prop(cid, properties):
    """
    Parses PubChem database using molecule cid and outputs the requested
    properties - if available - (works for selected properties only) in a
    dictionary format.

    Parameters
    ----------
    cid: int
        PubChem Compound CID.
    properties: list of str
        Properties to parse from PubChem database. Supported properties include:
             Auto-Ignition
             Boiling Point
             Density
             Dissociation Constants
             Flash Point
             Heat of Combustion
             Heat of Vaporization
             Ionization Potential
             LogP
             Melting Point
             Solubility
             Stability
             Surface Tension
             Vapor Density
             Vapor Pressure
             Viscosity
             Wikipedia

    Returns
    -------
    available_properties: dict
        Dictionary of parsed properties.
    """
    rest_api_link = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/{cid}/JSON"
    )
    output = parse_pubchem(rest_api_link)["Record"]["Section"]

    output_data = []
    for heading in ["Names and Identifiers", "Chemical and Physical Properties"]:
        output_data.append(
            list(filter(lambda section: section["TOCHeading"] == heading, output))[0]
        )
    identifiers = [
        "Computed Descriptors",
        "Other Identifiers",
        "Experimental Properties",
    ]
    all_data = []
    for dict_ in output_data:
        all_data += list(
            filter(
                lambda section: section["TOCHeading"] in identifiers, dict_["Section"]
            )
        )

    required_data = []
    for data in all_data:
        required_data += list(
            filter(lambda section: section["TOCHeading"] in properties, data["Section"])
        )

    available_properties = {}
    for prop in required_data:
        available_properties[prop["TOCHeading"]] = prop["Information"]
    return available_properties


def filter_pubchem_molecules(cids, quant_criteria=None, atom_types=None):
    """
    Filters a list of molecules represented by PubChem CIDs by parsing their
    properties from PubChem, comparing them to user defined criteria, and
    returning a modified list of molecules with only molecules matching the
    criteria. Works only with quantitative properties defined as int or float
    and with PubChem computed (table) properties.

    Parameters
    ----------
    cids: list of int
        PubChem Compound CIDs.
    quant_criteria: dictionary, optional
        Keys correspond to PubChem computed properties and values correspond to
        a list of lower and upper limit of each computed property,
        e.g. {"Charge": [-1, 0]}
    atom_types: set of str, optional
        Set of atom types to include in each molecule. Any molecule containing
        additional atom types will be excluded.

    Returns
    -------
    filtered_cids: list of int
        Filtered PubChem Compound CIDs.
    """
    filtered_cids = []
    if atom_types:
        properties = list(quant_criteria.keys()) + ["MolecularFormula"]
    else:
        properties = list(quant_criteria.keys())

    for cid in cids:
        computed_prop = get_pubchem_computed_prop(cid, properties)
        match = True
        for key, value in computed_prop.items():
            if key not in ["CID", "MolecularFormula"]:
                if quant_criteria[key][0] <= value <= quant_criteria[key][1]:
                    match = True
                else:
                    match = False
                    break
            if atom_types:
                atoms = set(
                    "".join(
                        [
                            i
                            for i in computed_prop["MolecularFormula"]
                            if not i.isdigit()
                        ]
                    )
                )
                if all(flag in atom_types for flag in atoms):
                    match = True
                else:
                    match = False
        if match:
            filtered_cids.append(cid)
    return filtered_cids
