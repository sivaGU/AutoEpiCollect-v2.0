import re
import pandas as pd

# 1) Read the text from a file
with open("vaxijen_test.txt", "r") as file:
    text = file.read()

# Initialize lists to store parsed data
peptides = []
antigenicity_scores = []

# Regular expression patterns
peptide_pattern = r'>peptide_\d+\s+([A-Z]+)'
antigenicity_pattern = r'Overall Prediction for the Protective Antigen = ([-+]?\d*\.?\d+)'

# Find all peptide sequences
peptide_matches = re.findall(peptide_pattern, text)

# Find all antigenicity scores
antigenicity_matches = re.findall(antigenicity_pattern, text)

# Ensure we have the same number of peptides and scores
if len(peptide_matches) != len(antigenicity_matches):
    print(f"Warning: Found {len(peptide_matches)} peptides but {len(antigenicity_matches)} scores")

# Create the DataFrame
df = pd.DataFrame({
    'peptide': peptide_matches,
    'antigenicity': [float(score) for score in antigenicity_matches]
})


# 5) Save to Excel
output_file = "peptide_antigenicity.xlsx"
df.to_excel(output_file, index=False)

print(f"Saved {len(df)} records to '{output_file}'.")


