import pandas as pd
from gtfparse import read_gtf
import pickle

def get_largest_transcript_and_utr(gtf_file, output_file, pickle_file):
    # Specify the data types for each column
    dtype = {
        'seqname': str,
        'source': str,
        'feature': str,
        'start': int,
        'end': int,
        'score': str,
        'strand': str,
        'frame': str,
        'attribute': str
    }

    # Read the GTF file into a pandas DataFrame
    df = pd.read_csv(gtf_file, sep='\t', comment='#', header=None, dtype=dtype, low_memory=False,
                     names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])

    # Save the GTF data to a pickle file
    with open(pickle_file, 'wb') as f:
        pickle.dump(df, f)

    # Filter for protein-coding transcripts
    transcript_df = df.query("feature == 'transcript' and 'protein_coding' in attribute")

    # Extract gene_id, gene_name, and transcript_id from the attribute column
    transcript_df['gene_id'] = transcript_df['attribute'].str.extract('gene_id "([^"]+)"')
    transcript_df['gene_name'] = transcript_df['attribute'].str.extract('gene_name "([^"]+)"')
    transcript_df['transcript_id'] = transcript_df['attribute'].str.extract('transcript_id "([^"]+)"')

    # Calculate the length of each transcript
    transcript_df['length'] = transcript_df['end'] - transcript_df['start'] + 1

    # Rename seqname to Chr
    transcript_df['Chr'] = transcript_df['seqname']

    # Find the largest transcript for each gene
    largest_transcript_df = transcript_df.loc[transcript_df.groupby('gene_id')['length'].idxmax()]

    # Function to calculate 3'UTR length
    def calculate_3utr_length(transcript_id):
        utr3 = df.query("feature == 'UTR' and transcript_id == @transcript_id and 'three_prime_UTR' in attribute")
        return (utr3['end'] - utr3['start'] + 1).sum()

    # Calculate 3'UTR lengths for the largest transcripts
    largest_transcript_df['3UTR_Length'] = largest_transcript_df['transcript_id'].apply(calculate_3utr_length)

    # Select specific columns to write to the output file
    selected_columns = ['transcript_id', 'gene_name', 'gene_id', 'Chr', 'start', 'end', 'length', '3UTR_Length']
    largest_transcript_df[selected_columns].to_csv(output_file, sep='\t', index=False)

    print(f"Largest transcript and 3'UTR information has been written to {output_file}")
    print(f"GTF data has been saved to {pickle_file}")

# Example usage
gtf_file = '/research_jude/rgs01_jude/groups/blancgrp/projects/rRNA_variation/common/EMT_analysis_Nitish/data/mouse/ensembl/Mus_musculus.GRCm39.104.rdna_rn18s.gtf'
output_file = 'largest_transcript_and_utr.txt'
pickle_file = 'gtf_data.pkl'
get_largest_transcript_and_utr(gtf_file, output_file, pickle_file)

