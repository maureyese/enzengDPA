from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO
import os
import random
import string
import subprocess
import pandas as pd
from pymsaviz import MsaViz

# ------------------------
# Run MSA
# ------------------------

def perform_msa(
        input_fasta: str,
        rid: str,
        num_sequences: int = 100,
        export_folder: str = "files/"
) -> AlignIO.MultipleSeqAlignment:
    '''
    Run MSA using Clustal-Omega using a specific input FASTA file.
    '''
    # 1. Verify input file exists
    if not os.path.exists(input_fasta):
        print(f"Error: The file {input_fasta} was not found.")
        return None

    if not export_folder.endswith("/"):
        export_folder = f"{export_folder}/"

    # 2. Subset the sequences
    # We read the input_fasta and take only the top N sequences
    subset_fasta = f"{export_folder}{rid}_msa_input.fasta"
    temp_aln_out = f"{export_folder}{rid}_aligned.fasta"
    
    try:
        records = list(SeqIO.parse(input_fasta, "fasta"))
        # Ensure we don't try to take more sequences than available
        actual_num = min(len(records), num_sequences)
        SeqIO.write(records[:actual_num], subset_fasta, "fasta")
        
        print(f"Prepared {actual_num} sequences for alignment.")

        # 3. Run Clustal Omega
        # -i: input file, -o: output file, --force: overwrite, --auto: set parameters automatically
        cmd = [
            "clustalo", 
            "-i", subset_fasta, 
            "-o", temp_aln_out, 
            "--outfmt=fasta", 
            "--force", 
            "--auto"
        ]
        
        print(f"Executing: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Clustal Omega Error: {result.stderr}")
            return None

        # 4. Load the result into Biopython AlignIO
        alignment = AlignIO.read(temp_aln_out, "fasta")
        print(f"Alignment completed successfully with {len(alignment)} sequences.")

        # 5. Optional Cleanup
        if os.path.exists(subset_fasta):
            os.remove(subset_fasta)

        return alignment

    except Exception as e:
        print(f"An error occurred during MSA: {e}")
        return None

# ------------------------
# Testing functions
# ------------------------

if __name__ == "__main__":
    df = pd.read_csv("files/testing_blast_results.csv")
    sequence_id = "SPY12701.1"
    num_sequences = 100
    alignment = perform_msa(
        df=df,
        query_id=sequence_id,
        num_sequences=num_sequences
    )
    # Visualize MSA
    #mv = MsaViz(alignment, color_scheme="Taylor", wrap_length=80, show_grid=True, show_consensus=True)
    #mv.savefig(f"plots/{sequence_id}_msa.png")
    #mv.plotfig()