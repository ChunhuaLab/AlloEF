import pandas as pd
import os
import re
from pathlib import Path


def parse_pdb_to_excel(pdb_file_path):
    """
    Read PDB file and extract protein information, generate Excel file

    Parameters:
    pdb_file_path: PDB file path
    """

    # Check if file exists
    if not os.path.exists(pdb_file_path):
        print(f"Error: File not found {pdb_file_path}")
        return

    # Get protein name from filename
    protein_name = Path(pdb_file_path).stem

    # Store extracted information
    residue_data = []

    try:
        with open(pdb_file_path, 'r', encoding='utf-8') as file:
            lines = file.readlines()

        print(f"Parsing PDB file: {pdb_file_path}")
        print(f"Protein name: {protein_name}")

        # Parse each line
        for line in lines:
            # Only process ATOM records focusing on CA (carbon alpha) atoms
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                try:
                    # Parse PDB format fixed column positions
                    atom_serial = line[6:11].strip()  # Atom serial number
                    atom_name = line[12:16].strip()  # Atom name (CA)
                    residue_name = line[17:20].strip()  # Residue name
                    chain_id = line[21].strip()  # Chain identifier
                    residue_seq = line[22:26].strip()  # Residue sequence number
                    x_coord = float(line[30:38].strip())  # X coordinate
                    y_coord = float(line[38:46].strip())  # Y coordinate
                    z_coord = float(line[46:54].strip())  # Z coordinate

                    # Add to data list
                    residue_data.append({
                        'Atom_Serial': atom_serial,
                        'Residue_Name': residue_name,
                        'Residue_Number': residue_seq,
                        'Chain_Name': chain_id,
                        'X_Coordinate': x_coord,
                        'Y_Coordinate': y_coord,
                        'Z_Coordinate': z_coord
                    })

                except (ValueError, IndexError) as e:
                    print(f"Warning: Error parsing line, skipping: {line.strip()}")
                    continue

        # Create DataFrame
        if residue_data:
            df = pd.DataFrame(residue_data)

            # Generate Excel filename
            excel_filename = f"{protein_name}_protein_info.xlsx"

            # Create Excel writer object
            with pd.ExcelWriter(excel_filename, engine='openpyxl') as writer:
                # Write main data
                df.to_excel(writer, sheet_name='Protein_Residue_Info', index=False)

                # Create statistics table
                stats_data = {
                    'Statistics': ['Total_Residues', 'Chain_Count', 'Protein_Name', 'Source_File'],
                    'Value': [len(df), len(df['Chain_Name'].unique()),
                           protein_name, os.path.basename(pdb_file_path)]
                }
                stats_df = pd.DataFrame(stats_data)
                stats_df.to_excel(writer, sheet_name='Statistics', index=False)

                # Chain-grouped statistics
                chain_stats = df.groupby('Chain_Name').size().reset_index()
                chain_stats.columns = ['Chain_Name', 'Residue_Count']
                chain_stats.to_excel(writer, sheet_name='Chain_Statistics', index=False)

            print(f"Successfully generated Excel file: {excel_filename}")
            print(f"Total extracted {len(df)} residue information")
            print(f"Involving {len(df['Chain_Name'].unique())} protein chains: {', '.join(sorted(df['Chain_Name'].unique()))}")

            # Display data preview
            print("\nData preview:")
            print(df.head())

        else:
            print("Error: No valid CA atom records found")

    except Exception as e:
        print(f"Error occurred while processing file: {e}")


def batch_process_pdb_files(directory_path):
    """
    Batch process all PDB files in directory

    Parameters:
    directory_path: Directory path containing PDB files
    """
    pdb_files = []

    # Find all PDB files
    for ext in ['*.pdb', '*.PDB']:
        pdb_files.extend(Path(directory_path).glob(ext))

    if not pdb_files:
        print(f"No PDB files found in directory {directory_path}")
        return

    print(f"Found {len(pdb_files)} PDB files, starting batch processing...")

    for pdb_file in pdb_files:
        print(f"\n{'='*50}")
        parse_pdb_to_excel(str(pdb_file))


# Usage examples
if __name__ == "__main__":
    parse_pdb_to_excel("1WQW_A.pdb")
    # batch_process_pdb_files("./pdb_files/")