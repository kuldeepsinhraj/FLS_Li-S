# ebatt

Code and data for the computational screening of fluorinated ether solvents (FLS) for lithium–sulfur (Li–S) batteries, combining PubChem structure search, DFT, COSMO-RS, molecular dynamics (MD), and machine learning (XGBoost/Random Forest/Ridge + SHAP).


## What's here

A library of fluorinated ether solvents is built from PubChem, screened down using COSMO-RS solubility/boiling-point/flash-point cutoffs and DFT-derived redox and binding-energy criteria, and interpreted with XGBoost/SHAP over three descriptor sets (17 custom descriptors, +RDKit, +Mordred). The most promising candidates are then evaluated with MD simulations and validated experimentally.

```
ebatt/            Python package: descriptor calculation, ML models, and plotting utilities
  data/           PubChem queries, DFT/MD result parsing
  descriptors/    Structural descriptors, DFT-based descriptors, and the ML pipeline 
  screening/      Library down-selection logic
  visualization/  Shared plotting helpers
  config/         Solvent lists, screening thresholds

scripts/          All scripts used for the simulations and calculations.
  plots/          Figure-generation scripts

files/            Data: PubChem structures (.sdf/.png), Gaussian inputs (.com), optimized
                  geometries (.xyz), and processed results (files/results/)
```

## Installation

```bash
git clone https://github.com/kuldeepsinhraj/FLS_Li-S.git
cd FLS_Li-S
pip install -e .
```

Requires Python 3.8+. `setup.py` covers the dependencies needed for descriptor generation, screening, and the ML/SHAP analysis (RDKit, Mordred, OpenBabel, pymatgen, sterimol, scikit-learn, XGBoost, SHAP).

Gaussian, COSMORS, and LAMMPS are required to reproduce the raw DFT/COSMO-RS/MD calculations themselves 

