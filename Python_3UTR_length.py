import pandas as pd
import gffutils

# Path to your GFF file
gff_file = "/research_jude/rgs01_jude/groups/blancgrp/projects/rRNA_variation/common/EMT_analysis_Nitish/data/mouse/ensembl/Mus_musculus.GRCm39.104.rdna_rn18s.gtf"

# Create a database from the GFF file
db = gffutils.create_db(gff_file, dbfn=':memory:', force=True, keep_order=True, merge_strategy='merge', sort_attribute_values=True, disable_infer_transcripts=True)


# Extract 3' UTR features and calculate lengths for protein-coding genes
utr_lengths = []
for feature in db.features_of_type('three_prime_UTR'):
    parent_id = feature.attributes.get('Parent', [None])[0]
    parent_feature = db[parent_id]
    if parent_feature.featuretype == 'mRNA' and 'protein_coding' in parent_feature.attributes.get('gene_biotype', []):
        utr_length = len(feature)
        utr_lengths.append((parent_id, utr_length))

# Convert to DataFrame
utr_lengths_df = pd.DataFrame(utr_lengths, columns=['gene_id', 'utr_length'])

# Save to CSV
utr_lengths_df.to_csv("three_prime_utr_lengths_protein_coding.csv", index=False)

print("3' UTR lengths for protein-coding genes calculated and saved to 'three_prime_utr_lengths_protein_coding.csv'")

