import pandas as pd


master_df = pd.read_excel("master_epitope_table.xlsx")
antigenicity_df = pd.read_excel("peptide_antigenicity.xlsx")

merged_df = pd.merge(master_df, antigenicity_df, on="peptide", how="left")

merged_df.to_excel("full_master_epitope_table.xlsx", index=False)
