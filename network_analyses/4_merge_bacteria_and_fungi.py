import pandas as pd

# Load the datasets
its_df = pd.read_csv('merged_and_cleaned_ITS_data.csv')
ssu_df = pd.read_csv('merged_and_cleaned_16S_data.csv')

# Identify overlapping metadata columns 
common_cols = its_df.columns.intersection(ssu_df.columns).tolist()
common_cols.remove('sample-id')

# Merge the dataframes
plant_bacteria_data_merged = pd.merge(
    ssu_df, 
    its_df.drop(columns=common_cols), 
    on='sample-id', 
    how='inner'
)

# Display information about the merged dataset
print(f"Bacterial samples: {ssu_df.shape[0]}")
print(f"Fungal samples: {its_df.shape[0]}")
print(f"Merged samples: {plant_bacteria_data_merged.shape[0]}")
print(f"Total columns in merged file: {plant_bacteria_data_merged.shape[1]}")

# Save the merged file
plant_bacteria_data_merged.to_csv('plant_bacteria_and_fungi_data_merged.csv', index=False)
