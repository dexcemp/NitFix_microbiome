import pandas as pd
import re

# 1. Load the metadata and abundance table
metadata = pd.read_csv('sample-metadata.16S.csv')
table = pd.read_csv('table.hostremoved-l6.16S.csv')

# 2. Filter metadata for only 'nodule' samples
metadata = metadata[metadata['sample_type'] == 'nodule'].copy()

# 3. Prepare the metadata
# Renaming 'genus' to 'host_genus' to avoid overlap with bacterial genus names
metadata = metadata.rename(columns={'genus': 'host_genus'})

# 4. Prepare the abundance table
table.set_index(table.columns[0], inplace=True)
table_transposed = table.T
table_transposed.index.name = 'sample-id'

# 5. Filter for Bacterial taxa and extract Genus names
def extract_genus(tax_string):
    match = re.search(r'g__([^;]*)', tax_string)
    if match:
        genus = match.group(1).strip()
        if genus:
            return genus.replace('[', '').replace(']', '')
    
    parts = [p.split('__')[-1].strip() for p in tax_string.split(';')]
    for p in reversed(parts):
        if p:
            return f"Unclassified_{p.replace('[', '').replace(']', '')}"
    return "unclassified"

# Process only bacterial columns
bacterial_cols = [c for c in table_transposed.columns if str(c).startswith('k__Bacteria')]
tax_data = table_transposed[bacterial_cols].copy()
tax_data.columns = [extract_genus(col) for col in tax_data.columns]

# 6. Aggregate duplicates
tax_data_grouped = tax_data.T.groupby(level=0).sum().T

# 7. Filter: Strict removal of Unclassified and Numbers
# - Case-insensitive check for "unclassified"------ Unclassified_<rank>, such as Unclassified_Rhodospirillaceae, Unclassified_Enterobacteriaceae and Unclassified_Comamonadaceae.
# - Removes any column containing a digit (0-9) ------ 1-68, 4-29, A17, BSV43, DA101, Dok59, Ellin506, FFCH10602, G07, GOUTA19, HB118, JG37-AG-70, KSA1, LCP-6, OR-59, R18-435, SC3-56, WAL_1855D, WCHB1-84, ph2, rc4-4
final_tax_cols = [
    col for col in tax_data_grouped.columns 
    if "unclassified" not in col.lower() and not any(char.isdigit() for char in col)
]
tax_data_filtered = tax_data_grouped[final_tax_cols].copy()

# 8. Convert abundance to 0-1 (Relative Abundance)
row_sums = tax_data_filtered.sum(axis=1)
tax_data_rel = tax_data_filtered.div(row_sums, axis=0).fillna(0)

# 9. Merge with Metadata
merged_df = pd.merge(metadata, tax_data_rel.reset_index(), on='sample-id', how='inner')

# 10. Save the result with the new requested filename
merged_df.to_csv('merged_and_cleaned_16S_data.csv', index=False)

# 11. Output simple completion message
print("Done")
