import pandas as pd
from pathlib import Path
import glob

# 1) Find all your per-gene xlsx files:
xlsx_files = ["all_variables_mhci_AHNAK2.xlsx", "all_variables_mhci_ARID1A.xlsx"]

all_dfs = []
for file in xlsx_files:
    path = Path(file)
    stem = path.stem  # e.g. "all_epitopes_mhci_AHNAK2"

    # extract gene = last underscore-delimited segment
    gene = stem.split("_")[-1]  # → "AHNAK2", "ARID1A", etc.
    sheets = pd.read_excel(file, sheet_name=None)
    for mutation_name, df in sheets.items():
        df = df.copy()
        if "Cancers" in df.columns:
            df = df.drop(columns=["Cancers"])
        df["mutation"] = mutation_name
        df["gene"] = gene
        all_dfs.append(df)

# 2) Concatenate into one master table
master = pd.concat(all_dfs, ignore_index=True)

# 4) (Optional) reorder columns so 'gene' and 'mutation' come first
cols = ["gene", "mutation"] + [c for c in master.columns if c not in ("gene", "mutation")]
master = master[cols]

# 5) Save to disk
master.to_excel("master_epitope_table.xlsx", index=False)

# 3) Remove duplicate peptide sequences
master = master.drop_duplicates(subset=["peptide"], keep="first")

# 5) Save to disk
master.to_excel("master_epitope_table_deduped.xlsx", index=False)

with open("all_peptides_anti.fasta", "w") as fh:
    for idx, pep in enumerate(master["peptide"], start=1):
        # you can customize the header as you like
        fh.write(f">peptide_{idx}\n")
        fh.write(f"{pep}\n")
