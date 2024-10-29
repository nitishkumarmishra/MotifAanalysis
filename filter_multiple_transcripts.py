from Bio import SeqIO
import pandas as pd

# Define file paths
#fasta_files = ["m39_transcript_region_fasta/three_prime_utr.fasta", "m39_transcript_region_fasta/CDS.fasta"]
fasta_files = [
    "m39_transcript_region_fasta/5UTR_CDS.fasta",
    "m39_transcript_region_fasta/5UTR_start30.fasta",
    "m39_transcript_region_fasta/CDS.fasta",
    "m39_transcript_region_fasta/five_prime_utr.fasta",
    "m39_transcript_region_fasta/three_prime_utr.fasta"
]

#pc_name_file = "m39_canonicalTranscriptsGeneName.txt"
output_files = [
        "m39_transcript_region_fasta/protein_coding_5UTR_CDS.fasta",
        "m39_transcript_region_fasta/protein_coding_5UTR_start30.fasta",
        "m39_transcript_region_fasta/protein_coding_CDS.fasta", 
        "m39_transcript_region_fasta/protein_coding_five_prime_utr.fasta",
        "m39_transcript_region_fasta/protein_coding_three_prime_utr.fasta",

        ]
pc_name_file = "m39_canonicalTranscriptsGeneName.txt"


# Read the sequence names into a list
pc_gene_names = pd.read_csv(pc_name_file, sep='\t')['transcript_id'].tolist()

# Process each FASTA file
for fasta_file, output_file in zip(fasta_files, output_files):
    with open(output_file, "w") as out_f:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in pc_gene_names:
                SeqIO.write(record, out_f, "fasta")
    print(f"Selected sequences from {fasta_file} have been written to {output_file}")


