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
        df:pd.DataFrame,
        query_id:str,
        num_sequences: int = 100) -> AlignIO.MultipleSeqAlignment:
    '''Run MSA using Clustal-Omega. Make sure you have it installed (Ubuntu-based)'''
    # Verify df is not empty
    if df.empty:
        print("Error: Input DataFrame is empty.")
        return None

    # Prepare sequences
    sequences = []

    # Extract original query sequence (it's the same for all hits)
    query_seq_str = df.iloc[0]['Query_seq']
    query_record = SeqRecord(
        Seq(query_seq_str),
        id=query_id,
        name=query_id,
        description=f"Query Sequence: {query_id}"
    )

        # Extract top N hit sequences (Hsp_hseq)
    # Ensure num_sequences is at least 1 (for the query)
    num_hits = max(0, num_sequences - 1)
    top_hits = df.head(num_hits)

    for index, row in top_hits.iterrows():
        # Hit_id is typically the most descriptive identifier
        accession_parts = row['Hit_id'].split('|')
        # Tries to get the accession number if the ID is piped (e.g., ref|ACC|def)
        seq_id = accession_parts[1] if len(accession_parts) > 1 else row['Hit_accession']
        seq_def = row['Hit_def']

        hit_record = SeqRecord(
            Seq(row['Hit_seq']),
            id=seq_id,
            name=seq_id,
            description=f"{seq_id} | {seq_def}"
        )
        sequences.append(hit_record)

    if len(sequences) < 2:
        print(f"Error: Only {len(sequences)} sequence(s) available. Need at least two sequences for alignment.")
        return None

    print(f"Prepared {len(sequences)} sequences (1 Query + {len(sequences) - 1} Top Hits) for alignment.")

    # 2. Save sequences to a temporary FASTA file
    random_suffix = ''.join(random.choices(string.ascii_letters + string.digits, k=12))
    temp_fasta_in = f"temp_input_{random_suffix}.fasta"
    temp_aln_out = f"temp_output_{random_suffix}.fasta"

    with open(temp_fasta_in, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
    print(f"Sequences saved to temporary file: {temp_fasta_in}")

    # 3. Run the Clustal alignment
    try:
        cmd = ["clustalo", "-i", temp_fasta_in, "-o", temp_aln_out, "--outfmt=fasta", "--verbose", "--auto"]
        print(f"Running: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(f"Clustal Omega failed with error: {result.stderr}")
            return None

        print("Alignment completed successfully")

        # 4. Read the alignment
        alignment = AlignIO.read(temp_aln_out, "fasta")

        # 5. Clean up temporary files
        os.remove(temp_fasta_in)
        os.remove(temp_aln_out)
        print("Temporary files cleaned up.")

        return alignment

    except Exception as e:
        print(f"Error during alignment: {e}")
        # Clean up on error
        if os.path.exists(temp_fasta_in):
            os.remove(temp_fasta_in)
        if os.path.exists(temp_aln_out):
            os.remove(temp_aln_out)
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