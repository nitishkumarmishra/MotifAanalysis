"""
This script converts HOMER motif files from a specified directory into a single MEME format file.
It reads each HOMER motif file, extracts the motif name and position weight matrix (PWM), and writes
the information into a MEME format file.

Key Components:
1. Importing necessary modules:
   - os: For interacting with the file system.

2. Functions:
   - read_homer_motif(file_path): Reads a HOMER motif file and extracts the motif name and PWM.
   - write_meme_header(output_file): Writes the header for the MEME format file.
   - write_meme_motif(output_file, motif_name, pwm): Appends each motif's information to the MEME format file.
   - convert_homer_to_meme(homer_dir, output_file): Iterates through all HOMER motif files in the specified directory and converts them to MEME format.

3. Main Execution:
   - Defines the directory containing HOMER motif files and the output MEME file.
   - Calls the convert_homer_to_meme function to perform the conversion.

Usage:
- Ensure the HOMER motif files are in the specified directory.
- Run the script to generate a single MEME format file containing all the motifs.
"""

import os

def read_homer_motif(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    #motif_name = lines.strip().split()
    first_line = lines[0]
    #print("First line:", first_line)

    split_line = first_line.split()
    #print("Split line:", split_line)

    motif_name = split_line[1]
    print("Motif name:", motif_name)

    pwm = []
    for line in lines[1:]:
        if line.strip():
            pwm.append([float(x) for x in line.strip().split()])
    
    return motif_name, pwm

def write_meme_header(output_file):
    with open(output_file, 'w') as f:
        f.write("MEME version 4\n\n")
        f.write("ALPHABET= ACGT\n\n")
        f.write("strands: + -\n\n")
        f.write("Background letter frequencies:\n")
        f.write("A 0.25 C 0.25 G 0.25 T 0.25\n\n")

def write_meme_motif(output_file, motif_name, pwm):
    with open(output_file, 'a') as f:
        f.write(f"MOTIF {motif_name} {motif_name}\n\n")
        f.write("letter-probability matrix: alength= 4 w= {} nsites= 20 E= 0\n".format(len(pwm)))
        for row in pwm:
            f.write("  ".join(f"{x:.6f}\t " for x in row) + "\n")
        f.write("\n")

def convert_homer_to_meme(homer_dir, output_file):
    write_meme_header(output_file)
    
    for file_name in os.listdir(homer_dir):
        if file_name.endswith(".motif"):
            file_path = os.path.join(homer_dir, file_name)
            motif_name, pwm = read_homer_motif(file_path)
            write_meme_motif(output_file, motif_name, pwm)
    
    print(f"All HOMER motifs have been converted to MEME format in {output_file}")

# Define the directory containing HOMER motif files and the output MEME file
homer_dir = 'homer_motifs/'
output_file = 'combined_motifs.meme'

# Convert HOMER motif files to MEME format
convert_homer_to_meme(homer_dir, output_file)


