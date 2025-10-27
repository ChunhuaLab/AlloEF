import os
import subprocess
import pandas as pd

# Define input folder containing PDB files
input_folder = r"./Independent_test_set"#1WQW_A.pdb
# Define output folder path
output_folder = r"./Independent_test_set_result"#1WQW_A_dssp_result

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

# Iterate through all PDB files in the input folder
for pdb_file in os.listdir(input_folder):
    if pdb_file.endswith(".pdb"):
        # Get filename without extension
        file_name = os.path.splitext(pdb_file)[0]
        # Define corresponding DSSP file path
        dssp_file = os.path.join(output_folder, f"{file_name}.dssp")

        # Complete input PDB file path
        pdb_file_path = os.path.join(input_folder, pdb_file)

        # Build and run DSSP command
        cmd = ["dssp", "-i", pdb_file_path, "-o", dssp_file]
        result = subprocess.run(cmd, capture_output=True, text=True)

        # Check if command was successful
        if result.returncode == 0:
            print(f"Successfully processed {pdb_file} 🌟")

            all_data = []

            # Read DSSP file content and extract information
            with open(dssp_file) as f:
                lines = f.readlines()
                parse_data = False
                for line in lines:
                    if line.startswith("  #  RESIDUE AA STRUCTURE BP1 BP2  ACC"):
                        parse_data = True
                        continue
                    if parse_data:
                        # Remove leading spaces and split fields at fixed intervals
                        if len(line) >= 120:
                            residue_data = {
                                "RESIDUE": line[5:10].strip(),  # Residue number and chain name/ID
                                "AA": line[13:14].strip(),  # Amino acid type
                                "STRUCTURE": line[16:25].strip(),  # Secondary structure type
                                "BP1": line[25:29].strip(),
                                "BP2": line[29:33].strip(),
                                "ACC": line[34:38].strip(),
                                "N-H-->O": line[39:50].strip(),
                                "O-->H-N": line[50:61].strip(),
                                "N-H-->O (2)": line[61:72].strip(),
                                "O-->H-N (2)": line[72:83].strip(),
                                "TCO": line[84:91].strip(),
                                "KAPPA": line[91:97].strip(),
                                "ALPHA": line[97:103].strip(),
                                "PHI": line[103:109].strip(),
                                "PSI": line[109:115].strip(),
                                "X-CA": line[115:122].strip(),
                                "Y-CA": line[122:129].strip(),
                                "Z-CA": line[129:].strip()
                            }
                            all_data.append(residue_data)

            # Convert data to DataFrame
            df = pd.DataFrame(all_data)

            # Add "seq" column
            df.insert(0, 'seq', range(1, len(df) + 1))

            # Define output Excel file path for each protein
            excel_file = os.path.join(output_folder, f"{file_name}.xlsx")

            # Write DataFrame to Excel file
            df.to_excel(excel_file, index=False)
            print(f"Created Excel file for {file_name} with sequence numbers 🗂️")
        else:
            print(f"Failed to process {pdb_file} 🚫")
            print(result.stderr)

print("Batch processing complete 😊")