import pandas as pd
from gtfparse import read_gtf
import pickle

# Path to your GTF file
gtf_file = '/research_jude/rgs01_jude/groups/blancgrp/projects/rRNA_variation/common/EMT_analysis_Nitish/data/mouse/ensembl/Mus_musculus.GRCm39.104.rdna_rn18s.gtf'

# Read the GTF file
gtf_data = read_gtf(gtf_file)

# Save the GTF data to a pickle file
with open('gtf_data.pkl', 'wb') as f:
    pickle.dump(gtf_data, f)

# Filter for protein-coding transcripts
# Filter for protein-coding transcripts using the query method
protein_coding_transcripts = gtf_data.query("feature == 'transcript' and transcript_biotype == 'protein_coding'")

#protein_coding_transcripts = gtf_data[(gtf_data['feature'] == 'transcript') & (gtf_data['transcript_biotype'] == 'protein_coding')]

# Find the largest transcript for each gene
largest_transcripts = protein_coding_transcripts.loc[protein_coding_transcripts.groupby('gene_id')['end'].idxmax()]

# Function to calculate 3'UTR length
def calculate_3utr_length(transcript_id):
    utr3 = gtf_data[(gtf_data['feature'] == 'UTR') & (gtf_data['transcript_id'] == transcript_id) & (gtf_data['type'] == 'three_prime_UTR')]
    return (utr3['end'] - utr3['start'] + 1).sum()

# Calculate 3'UTR lengths for the largest transcripts
largest_transcripts['3UTR_Length'] = largest_transcripts['transcript_id'].apply(calculate_3utr_length)

# Select relevant columns
result = largest_transcripts[['gene_id', 'transcript_id', '3UTR_Length']]

# Save the result to a CSV file
output_file = '3utr_lengths.csv'
result.to_csv(output_file, index=False)

print(f"Output saved to {output_file} and GTF data saved to gtf_data.pkl")

