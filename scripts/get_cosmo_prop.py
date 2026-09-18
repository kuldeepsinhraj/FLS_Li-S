import pandas as pd

# read dataframe containing pubchem and computed structural properties
df = pd.read_csv("../files/results/solvents.csv")

# read cosmo dataframe and merge with pubchem and computed properties
df_cosmo = pd.read_csv("../files/results/cosmo.txt")
df = pd.merge(df, df_cosmo, on="CID")

# add density (g/L)
df["Density"] = 1 / (df["MolarVolume"] / df["MolecularWeight"])

df.to_csv("../files/results/solvents.csv", index=False)
