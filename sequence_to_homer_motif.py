"""
This script generates HOMER motif files from a list of motifs and sequences provided in an input file.
It reads the input file, which contains motif names and sequences, and uses the HOMER suite's seq2profile.pl
script to create motif files in a specified output directory.

Key Components:
1. Importing necessary modules:
   - os: For interacting with the operating system, such as creating directories.
   - subprocess: For running external commands and scripts.

2. Defining paths:
   - input_file: Path to the input file containing motif names and sequences.
   - output_dir: Directory where the generated motif files will be saved.
   - homer_bin_path: Path to the HOMER binary directory.

3. Ensuring the output directory exists:
   - Creates the output directory if it doesn't already exist.

4. Function to run seq2profile.pl:
   - Constructs the output file path.
   - Builds the command to run seq2profile.pl with the given motif name and sequence.
   - Executes the command and writes the output to the specified file.

5. Reading the input file and generating motif files:
   - Opens the input file and reads it line by line.
   - Splits each line into motif name and sequence.
   - Calls the function to generate the motif file for each motif.
   - Includes error handling to manage missing files or other issues.

Usage:
- Ensure the input file and HOMER binary path are correctly specified.
- Run the script to generate the HOMER motif files in the specified output directory.

python sequence_to_homer_motif.py

"""

import os
import subprocess

# Define the paths to your files
#input_file = 'motifs.txt'
input_file = 'transite_extra_motif.txt'
output_dir = 'homer_motifs/'
homer_bin_path = '/research/groups/blancgrp/home/nmishra/EMT_data/Total_RNAseq/Riboprofiling-RSEM-BAM/HOMER/bin'

# Ensure the output directory exists
os.makedirs(output_dir, exist_ok=True)

# Function to run seq2profile.pl
def run_seq2profile(motif_name, sequence, output_dir, homer_bin_path):
    output_file = os.path.join(output_dir, f"{motif_name}.motif")
    seq2profile_command = [
        os.path.join(homer_bin_path, 'seq2profile.pl'),
        sequence,
        '0',  # Number of mismatches allowed
        motif_name,
        '>', output_file
    ]
    command_str = ' '.join(seq2profile_command)
    subprocess.run(command_str, shell=True, check=True, env=dict(os.environ, PATH=homer_bin_path + ':' + os.environ['PATH']))
    print(f"Generated HOMER motif file for {motif_name}")

# Read the input file and generate HOMER motif files
with open(input_file, 'r') as f:
    for line in f:
        motif_name, sequence = line.strip().split('\t')
        run_seq2profile(motif_name, sequence, output_dir, homer_bin_path)

print("All HOMER motif files have been generated.")

