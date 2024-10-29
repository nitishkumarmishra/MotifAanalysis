"""
This script processes Position Probability Matrices (PPMs) to generate consensus sequences with ambiguous nucleotides.
It uses numpy for matrix operations and defines a function to convert PPMs to consensus sequences based on a given threshold.

Key Components:
1. Importing necessary modules:
   - numpy: For handling numerical operations and matrix manipulations.

2. Example PPMs:
   - ppm_PRTE: A PPM for the PRTE motif, sourced from a scientific article.
   - ppm_CERT: A PPM for the CERT motif, sourced from a scientific article.

3. Nucleotides and IUPAC codes:
   - nucleotides: List of standard nucleotides (A, C, G, T).
   - iupac_codes: Dictionary mapping combinations of nucleotides to their IUPAC codes for ambiguous bases.

4. Function to get consensus sequence with ambiguous nucleotides:
   - ppm_to_consensus: Converts a PPM to a consensus sequence using a specified probability threshold.

5. Generating consensus sequences:
   - consensus_sequence_prte: Consensus sequence for the PRTE motif.
   - consensus_sequence_cert: Consensus sequence for the CERT motif.

Usage:
- Ensure numpy is installed in your Python environment.
- Run the script to print the consensus sequences for the provided PPMs.
"""

import numpy as np

# Example PPMs
# PRTE/CERT Letter Probability Matrix : https://www.nature.com/articles/s41586-024-07781-7
# Polypurine motif
# https://www.nature.com/articles/s41467-023-36290-w 
# https://github.com/sherkinglee/RocA/blob/main/data/homer-RocA03-up-CDS-motif.txt
# https://github.com/sherkinglee/RocA/tree/main
ppm_polypurine = np.array([
    [0.686, 0.032, 0.248, 0.034],
    [0.447, 0.240, 0.138, 0.175],
    [0.037, 0.026, 0.884, 0.053],
    [0.709, 0.049, 0.192, 0.050],
    [0.028, 0.215, 0.550, 0.207],
    [0.733, 0.160, 0.037, 0.070],
    [0.040, 0.689, 0.207, 0.064],
    [0.049, 0.631, 0.192, 0.128],
    [0.677, 0.157, 0.045, 0.121],
    [0.149, 0.155, 0.218, 0.477]
])


# Nucleotides and IUPAC codes for ambiguous bases
nucleotides = ['A', 'C', 'G', 'T']
iupac_codes = {
    'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T',
    'AC': 'M', 'AG': 'R', 'AT': 'W', 'CG': 'S',
    'CT': 'Y', 'GT': 'K', 'ACG': 'V', 'ACT': 'H',
    'AGT': 'D', 'CGT': 'B', 'ACGT': 'N'
}

# Function to get consensus sequence with ambiguous nucleotides
def ppm_to_consensus(ppm, threshold=0.1):
    consensus = []
    for position in ppm:
        max_prob = np.max(position)
        bases = [nucleotides[i] for i, prob in enumerate(position) if max_prob - prob <= threshold]
        consensus.append(iupac_codes[''.join(sorted(bases))])
    return ''.join(consensus)

# Get consensus sequence with ambiguous nucleotides
#consensus_sequence_prte = ppm_to_consensus(ppm_PRTE, threshold=0.15)
#print("Consensus Sequence for PRTE:", consensus_sequence_prte)
#consensus_sequence_cert = ppm_to_consensus(ppm_CERT, threshold=0.15)
#print("Consensus Sequence for CERT:", consensus_sequence_cert)
consensus_sequence_polypurine = ppm_to_consensus(ppm_polypurine, threshold=0.15)
print("Consensus Sequence for CERT:", consensus_sequence_polypurine)

