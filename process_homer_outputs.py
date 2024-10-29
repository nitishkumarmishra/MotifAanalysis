import os
import pandas as pd

def read_homer_output(homer_file, protein_coding_transcripts):
    try:
        # Read HOMER output file without headers
        df = pd.read_csv(homer_file, sep='\t', header=None)
        
        # Assuming the first column contains the transcript IDs
        transcript_column = 0
        
        # Filter for protein-coding transcripts
        filtered_df = df[df[transcript_column].isin(protein_coding_transcripts)]
        
        # Debugging: Print the number of matches found
        print(f"File: {homer_file}, Matches found: {len(filtered_df)}")
        
        # Count the number of unique genes
        gene_count = filtered_df[transcript_column].nunique()
        
        return gene_count
    except pd.errors.EmptyDataError:
        print(f"File is empty or not formatted correctly: {homer_file}")
        return 0

def process_homer_outputs(homer_dir, motif_file, protein_coding_transcripts_file):
    # Read the list of protein-coding transcripts
    protein_coding_transcripts = pd.read_csv(protein_coding_transcripts_file, sep='\t')['transcript_id'].tolist()
    
    # Read the list of motif names and sequences
    motifs_df = pd.read_csv(motif_file, sep='\t')
    motif_names = motifs_df['motiffName'].tolist()
    
    # Initialize a dictionary to store counts
    counts = {motif: {'five_prime_utr': 0, 'three_prime_utr': 0, 'CDS': 0} for motif in motif_names}
    
    # Iterate through each motif name and region
    for motif in motif_names:
        for region in ['five_prime_utr', 'three_prime_utr', 'CDS']:
            homer_file = os.path.join(homer_dir, f"{region}_{motif}")
            if os.path.exists(homer_file):
                counts[motif][region] = read_homer_output(homer_file, protein_coding_transcripts)
            else:
                # Debugging: Print a message if the file does not exist
                print(f"File not found: {homer_file}")
    
    # Convert the dictionary to a DataFrame
    counts_df = pd.DataFrame.from_dict(counts, orient='index').reset_index()
    counts_df.columns = ['Motif', 'five_prime_utr', 'three_prime_utr', 'CDS']
    
    return counts_df

# Define the directory containing HOMER output files and the protein-coding transcripts file
homer_dir = 'homer_output'
protein_coding_transcripts_file = 'm39_canonicalTranscriptsGeneName.txt'
motif_file = 'transite_extra_motif.tsv'

# Process HOMER outputs and get counts
counts_df = process_homer_outputs(homer_dir, motif_file, protein_coding_transcripts_file)
print("Counts of protein-coding genes for each motif and region:")
print(counts_df)

# Save the DataFrame to a CSV file
counts_df.to_csv('motif_counts.csv', index=False)

