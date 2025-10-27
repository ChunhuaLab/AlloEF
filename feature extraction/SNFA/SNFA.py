import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer

# Read all data from Excel file
input_path = r'.\1WQW_Aall.xlsx'
df = pd.read_excel(input_path, sheet_name='Sheet1')

# Assume dataframe contains position coordinates 'x', 'y', 'z' and feature columns to process
feature_columns = [
    'NTE5', 'PRScol', 'PRSlin', 'ISPOCKET', 'Entropy',
    'Conservation Score', 'FrstIndex', 'ACC', 'NTECR_AVE',
    'NTECR_MAX', 'NTECR_MIN', 'NTECR_MID','msf'
]

# Replace all non-numeric entries with NaN
df[feature_columns] = df[feature_columns].apply(pd.to_numeric, errors='coerce')

# Now Imputer should work normally
# Handle missing values
imputer = SimpleImputer(strategy='mean')
df[feature_columns] = imputer.fit_transform(df[feature_columns])

grouped = df.groupby('ProteinID')

# Calculate spatial features for n=7 nearest neighbors only
n = 7
all_new_features = []

for protein_id, group in grouped:
    coordinates = group[['x', 'y', 'z']].values
    features = group[feature_columns].values
    num_residues = len(coordinates)

    new_features = np.zeros_like(features)

    for i in range(num_residues):
        distances = np.linalg.norm(coordinates - coordinates[i], axis=1)  # Calculate distances from current residue to all other residues
        nearest_indices = np.argsort(distances)[1:n + 1]  # Sort and select nearest 7 residues (excluding self)
        new_features[i] = features[nearest_indices].mean(axis=0)  # Calculate mean features of nearest 7 residues

    new_features_df = pd.DataFrame(new_features, columns=feature_columns)
    new_features_df['ProteinID'] = protein_id
    new_features_df.index = group.index  # Maintain original index

    all_new_features.append(new_features_df)

# Combine all new features
new_features_combined = pd.concat(all_new_features).sort_index()
for col in feature_columns:
    df[f'space_{col}_{n}nn'] = new_features_combined[col]

# Save results to new Excel file
output_path = rf'.\1WQW_A_allfeature{n}nn.xlsx'
df.to_excel(output_path, index=False)

print(f"Generated file: {output_path}")

print("Spatial feature engineering completed successfully!")