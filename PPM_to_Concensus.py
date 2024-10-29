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
ppm_PRTE = np.array([
    [0.000000, 0.526882, 0.064516, 0.408602],
    [0.000000, 0.204301, 0.021505, 0.774194],
    [0.010753, 0.376344, 0.000000, 0.612903],
    [0.000000, 0.537634, 0.000000, 0.462366],
    [0.053763, 0.870968, 0.000000, 0.075269],
    [0.000000, 0.000000, 0.000000, 1.000000],
    [0.000000, 0.182796, 0.000000, 0.817204],
    [0.000000, 0.290323, 0.000000, 0.709677],
    [0.000000, 0.720430, 0.021505, 0.258065],
    [0.075269, 0.419355, 0.139785, 0.365591],
    [0.086022, 0.451613, 0.268817, 0.193548],
    [0.129032, 0.279570, 0.139785, 0.451613],
    [0.000000, 0.301075, 0.387097, 0.311828],
    [0.086022, 0.419355, 0.161290, 0.333333],
    [0.118280, 0.322581, 0.182796, 0.376344]
])

ppm_CERT = np.array([
    [0.150000, 0.200000, 0.475000, 0.175000],
    [0.000000, 0.775000, 0.200000, 0.025000],
    [0.175000, 0.650000, 0.175000, 0.000000],
    [0.000000, 0.075000, 0.450000, 0.475000],
    [0.000000, 1.000000, 0.000000, 0.000000],
    [0.000000, 0.350000, 0.425000, 0.225000],
    [0.000000, 1.000000, 0.000000, 0.000000],
    [0.150000, 0.425000, 0.200000, 0.225000],
    [0.000000, 0.225000, 0.775000, 0.000000],
    [0.000000, 0.875000, 0.000000, 0.125000],
    [0.000000, 0.675000, 0.200000, 0.125000],
    [0.075000, 0.275000, 0.425000, 0.225000],
    [0.000000, 0.725000, 0.275000, 0.000000],
    [0.025000, 0.900000, 0.000000, 0.075000],
    [0.150000, 0.275000, 0.375000, 0.200000]
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
consensus_sequence_prte = ppm_to_consensus(ppm_PRTE, threshold=0.15)
print("Consensus Sequence for PRTE:", consensus_sequence_prte)
consensus_sequence_cert = ppm_to_consensus(ppm_CERT, threshold=0.15)
print("Consensus Sequence for CERT:", consensus_sequence_cert)

