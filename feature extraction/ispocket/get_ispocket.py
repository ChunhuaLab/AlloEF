import os
import glob
from openpyxl import Workbook


# Function: Read residue information from PDB file
def read_pdb_residues(pdb_file):
    """
    Extract residue information from PDB file

    Parameters:
    pdb_file: Path to PDB file

    Returns:
    set: Set of tuples containing (residue_name, residue_number)
    """
    residues = set()
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                res_name = line[17:20].strip()  # Residue name
                res_num = line[22:26].strip()  # Residue number
                residues.add((res_name, res_num))
    return residues


# Define main function
def main(source_directory, target_directory):
    """
    Main processing function

    Parameters:
    source_directory: Source directory containing subdirectories with PDB files
    target_directory: Target directory for output Excel files
    """
    # Create target directory if it doesn't exist
    os.makedirs(target_directory, exist_ok=True)

    # Iterate through all subdirectories in the main directory
    for subdir in os.listdir(source_directory):
        subdir_path = os.path.join(source_directory, subdir)
        if os.path.isdir(subdir_path):
            # Find PDB files matching the pattern in the subdirectory
            pdb_files = glob.glob(os.path.join(subdir_path, 'this_cavity_[12].pdb'))

            if pdb_files:
                # Create new workbook
                wb = Workbook()
                ws = wb.active
                ws.append(["ProteinID", "File_Name", "Residue_Name", "Residue_Number"])  # Add header

                collected_residues = set()  # Use set to track already added residues

                # Process each PDB file
                for pdb_file in pdb_files:
                    residues = read_pdb_residues(pdb_file)
                    file_name = os.path.basename(pdb_file)

                    # Add residue information to Excel
                    for res_name, res_num in residues:
                        unique_residue = (os.path.basename(subdir_path), file_name, res_name, res_num)
                        if unique_residue not in collected_residues:
                            ws.append(unique_residue)
                            collected_residues.add(unique_residue)

                    # Add empty row as separator
                    ws.append([])

                # Name output Excel file using subdirectory name
                output_file = os.path.join(target_directory, f'{os.path.basename(subdir_path)}.xlsx')
                wb.save(output_file)
                print(f"Excel file saved: {output_file}")


# Call main function with directory settings
if __name__ == "__main__":
    # Define source and target directories
    source_directory = r".\cavityplus2022\datas"
    target_directory = r".\result"#1WQW_A_result.csv

    # Execute main processing
    main(source_directory, target_directory)
    print("Processing completed successfully!")