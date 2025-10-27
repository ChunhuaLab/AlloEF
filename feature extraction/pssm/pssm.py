import os
import subprocess
import pandas as pd
from pathlib import Path

# Try to import openpyxl, install if not available
try:
    import openpyxl
except ImportError:
    os.system('pip install openpyxl')


class PSSMPipeline:
    def __init__(self, query_dir, database, output_dir, excel_output_dir, num_threads=40, num_iterations=3):
        """
        Initialize PSSM Pipeline

        Parameters:
        query_dir: Directory containing query FASTA files
        database: Path to BLAST database
        output_dir: Directory for PSSM output files
        excel_output_dir: Directory for Excel output files
        num_threads: Number of threads for PSI-BLAST
        num_iterations: Number of PSI-BLAST iterations
        """
        self.query_dir = query_dir
        self.database = database
        self.output_dir = output_dir
        self.excel_output_dir = excel_output_dir
        self.num_threads = num_threads
        self.num_iterations = num_iterations

        # Create output directories
        os.makedirs(self.output_dir, exist_ok=True)
        os.makedirs(self.excel_output_dir, exist_ok=True)

    def load_fasta_sequences(self, fasta_file):
        """
        Extract protein ID and sequences from FASTA file
        """
        sequences = {}
        with open(fasta_file, 'r') as file:
            protein_id = None
            sequence_lines = []
            for line in file:
                if line.startswith('>'):
                    if protein_id:
                        sequences[protein_id] = ''.join(sequence_lines)
                    protein_id = line.split()[0][1:]  # Extract protein ID
                    sequence_lines = []
                else:
                    sequence_lines.append(line.strip())
            if protein_id:
                sequences[protein_id] = ''.join(sequence_lines)
        return sequences

    def run_psiblast(self, protein_id, sequence):
        """
        Run PSI-BLAST for a single protein sequence
        """
        # Create temporary FASTA file
        temp_fasta_path = os.path.join(self.output_dir, f"{protein_id}_temp.fasta")
        with open(temp_fasta_path, 'w') as temp_fasta_file:
            temp_fasta_file.write(f">{protein_id}\n{sequence}\n")

        # Define output file paths
        ascii_pssm_out = os.path.join(self.output_dir, f"{protein_id}.pssm")
        checkpoint_out = os.path.join(self.output_dir, f"{protein_id}.chk")
        blast_out = os.path.join(self.output_dir, f"{protein_id}_results.out")

        # Build PSI-BLAST command
        psiblast_cmd = [
            "psiblast",
            "-query", temp_fasta_path,
            "-db", self.database,
            "-num_iterations", str(self.num_iterations),
            "-out_ascii_pssm", ascii_pssm_out,
            "-out_pssm", checkpoint_out,
            "-out", blast_out,
            "-num_threads", str(self.num_threads)
        ]

        print(f"Running PSI-BLAST for {protein_id}...")
        print(f"Command: {' '.join(psiblast_cmd)}")

        # Run PSI-BLAST
        process = subprocess.Popen(psiblast_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = process.communicate()

        # Clean up temporary file
        os.remove(temp_fasta_path)

        if process.returncode == 0:
            print(f"✅ PSSM for {protein_id} generated successfully")
            return ascii_pssm_out
        else:
            print(f"❌ Error processing {protein_id}: {stderr.decode('utf-8')}")
            return None

    def parse_pssm_to_excel(self, pssm_file_path, protein_id):
        """
        Parse PSSM file and convert to Excel format
        """
        try:
            # Read PSSM file
            with open(pssm_file_path, 'r') as file:
                lines = file.readlines()

            # Extract header (third line)
            header_line = lines[2].strip().split()
            headers = [header for header in header_line if header]  # Remove empty strings

            # Extract data (from 4th line to last 5 lines)
            data_lines = lines[3:-5]

            # Process data
            data = []
            max_cols = len(headers)  # Initial max columns
            for line in data_lines:
                row = line.strip().split()
                max_cols = max(max_cols, len(row))  # Update max columns
                data.append(row)

            # Pad rows to match max columns
            for i in range(len(data)):
                if len(data[i]) < max_cols:
                    data[i] += [''] * (max_cols - len(data[i]))

            # Pad headers if needed
            if len(headers) < max_cols:
                headers += [f'Extra_col_{i}' for i in range(len(headers), max_cols)]

            # Remove empty rows
            data = [row for row in data if any(cell.strip() for cell in row)]

            if data:
                # Create DataFrame
                pssm_df = pd.DataFrame(data, columns=headers)

                # Rename first two columns
                if len(pssm_df.columns) > 0:
                    pssm_df.columns.values[0] = 'seq'
                if len(pssm_df.columns) > 1:
                    pssm_df.columns.values[1] = 'residues_id'

                # Regenerate sequence numbers
                pssm_df['seq'] = range(1, len(pssm_df) + 1)

                # Save to Excel
                excel_file_path = os.path.join(self.excel_output_dir, f"{protein_id}.xlsx")
                pssm_df.to_excel(excel_file_path, index=False)
                print(f"✅ Excel file saved: {protein_id}.xlsx")
                return excel_file_path
            else:
                print(f"❌ No valid data found in PSSM file for {protein_id}")
                return None

        except Exception as e:
            print(f"❌ Error parsing PSSM file for {protein_id}: {e}")
            return None

    def process_batch_pssm_files(self):
        """
        Process existing PSSM files and convert to Excel (for files already generated)
        """
        print("Processing existing PSSM files...")

        # Group PSSM files by prefix (first 16 characters)
        pssm_file_dict = {}

        for pssm_file_name in os.listdir(self.output_dir):
            if pssm_file_name.endswith('.pssm'):
                prefix = pssm_file_name[:16]
                if prefix not in pssm_file_dict:
                    pssm_file_dict[prefix] = []
                pssm_file_dict[prefix].append(pssm_file_name)

        # Process each group
        for prefix, file_list in pssm_file_dict.items():
            combined_data = []
            headers = []

            # Sort files (files with '_1_' get priority)
            file_list.sort(key=lambda x: 1 if '_1_' in x else 2)

            for pssm_file_name in file_list:
                pssm_file_path = os.path.join(self.output_dir, pssm_file_name)

                try:
                    with open(pssm_file_path, 'r') as file:
                        lines = file.readlines()

                    header_line = lines[2].strip().split()
                    current_headers = [header for header in header_line if header]

                    data_lines = lines[3:-5]

                    data = []
                    max_cols = len(current_headers)
                    for line in data_lines:
                        row = line.strip().split()
                        max_cols = max(max_cols, len(row))
                        data.append(row)

                    for i in range(len(data)):
                        if len(data[i]) < max_cols:
                            data[i] += [''] * (max_cols - len(data[i]))

                    if len(current_headers) < max_cols:
                        current_headers += [f'Extra_col_{i}' for i in range(len(current_headers), max_cols)]

                    if not headers:
                        headers = current_headers

                    combined_data.extend(data)

                except Exception as e:
                    print(f"❌ Error reading {pssm_file_name}: {e}")

            # Save combined data
            combined_data = [row for row in combined_data if any(cell.strip() for cell in row)]

            if combined_data:
                pssm_df = pd.DataFrame(combined_data, columns=headers)

                if len(pssm_df.columns) > 0:
                    pssm_df.columns.values[0] = 'seq'
                if len(pssm_df.columns) > 1:
                    pssm_df.columns.values[1] = 'residues_id'

                pssm_df['seq'] = range(1, len(pssm_df) + 1)

                excel_file_path = os.path.join(self.excel_output_dir, f"{prefix}.xlsx")
                pssm_df.to_excel(excel_file_path, index=False)
                print(f"✅ Combined Excel file saved: {prefix}.xlsx")

    def run_complete_pipeline(self):
        """
        Run the complete pipeline: FASTA -> PSI-BLAST -> PSSM -> Excel
        """
        print("🚀 Starting complete PSSM pipeline...")

        processed_count = 0

        # Process all FASTA files
        for root, dirs, files in os.walk(self.query_dir):
            for query_file in files:
                if query_file.endswith(".fasta"):
                    input_path = os.path.join(root, query_file)
                    print(f"\n📁 Processing file: {input_path}")

                    sequences = self.load_fasta_sequences(input_path)

                    for protein_id, sequence in sequences.items():
                        print(f"\n🧬 Processing protein: {protein_id}")

                        # Run PSI-BLAST
                        pssm_file = self.run_psiblast(protein_id, sequence)

                        # Convert to Excel if PSSM was generated successfully
                        if pssm_file and os.path.exists(pssm_file):
                            self.parse_pssm_to_excel(pssm_file, protein_id)
                            processed_count += 1

        print(f"\n🎉 Pipeline completed! Processed {processed_count} proteins.")


# Usage example
if __name__ == "__main__":
    # Configuration
    query_dir = "./Independent_test_set"#1WQW_A.pdb
    database = "./datas/uniref50"
    output_dir = "./Independent_test_set_result"#1WQW_A_result
    excel_output_dir = "./Independent_test_set_excel"#1WQW_A_excel

    # Initialize pipeline
    pipeline = PSSMPipeline(
        query_dir=query_dir,
        database=database,
        output_dir=output_dir,
        excel_output_dir=excel_output_dir,
        num_threads=40,
        num_iterations=3
    )

    # Choose operation mode
    print("PSSM Processing Pipeline")
    print("1. Run complete pipeline (FASTA -> PSI-BLAST -> Excel)")
    print("2. Process existing PSSM files to Excel only")

    choice = input("\nSelect mode (1 or 2): ").strip()

    if choice == "1":
        pipeline.run_complete_pipeline()
    elif choice == "2":
        pipeline.process_batch_pssm_files()
    else:
        print("Invalid choice. Running complete pipeline by default.")
        pipeline.run_complete_pipeline()