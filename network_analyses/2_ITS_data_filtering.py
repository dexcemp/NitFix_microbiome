import pandas as pd
import re

# 1. Load the data
metadata = pd.read_csv('sample-metadata.ITS.csv')
table = pd.read_csv('table.hostremoved-l6.ITS.csv')

# 2. Filter metadata for only 'nodule' samples
metadata = metadata[metadata['sample_type'] == 'nodule'].copy()
metadata = metadata.rename(columns={'genus': 'host_genus'})

# 3. Prepare the abundance table
table.set_index(table.columns[0], inplace=True)
table_transposed = table.T
table_transposed.index.name = 'sample-id'

# 4. Extraction logic (Updated to handle "unidentified")
def extract_genus(tax_string):
    # Check for Genus tag 'g__'
    match = re.search(r'g__([^;]*)', tax_string)
    if match:
        genus = match.group(1).strip()
        # Ensure it's not empty, not underscores, and not 'unidentified'
        if genus and genus != '__' and genus != '' and genus.lower() != 'unidentified':
            return genus.replace('[', '').replace(']', '')
    
    # Fallback: get the last defined taxonomic level that isn't 'unidentified'
    parts = [p.split('__')[-1].strip() for p in tax_string.split(';')]
    for p in reversed(parts):
        # We check for '', '__', and 'unidentified'
        if p and p != '' and p != '__' and p.lower() != 'unidentified':
            return f"Unclassified_{p.replace('[', '').replace(']', '')}"
            
    return "unclassified"

# 5. Process Fungal columns
fungal_cols = [c for c in table_transposed.columns if str(c).startswith('k__Fungi')]
tax_data = table_transposed[fungal_cols].copy()
tax_data.columns = [extract_genus(col) for col in tax_data.columns]

# 6. Aggregate duplicate genera
tax_data_grouped = tax_data.T.groupby(level=0).sum().T

# 7. Final Clean (Removes "unclassified" AND "unidentified" prefixes + numbers)
final_tax_cols = [
    col for col in tax_data_grouped.columns 
    if "unclassified" not in col.lower() 
    and "unidentified" not in col.lower()
    and not any(char.isdigit() for char in col)
]
tax_data_filtered = tax_data_grouped[final_tax_cols].copy()

# 8. Relative Abundance (0-1)
row_sums = tax_data_filtered.sum(axis=1)
tax_data_rel = tax_data_filtered.div(row_sums, axis=0).fillna(0)

# 9. Merge and Save
merged_df = pd.merge(metadata, tax_data_rel.reset_index(), on='sample-id', how='inner')
merged_df.to_csv('merged_and_cleaned_ITS_data.csv', index=False)

print(f"Done.")
