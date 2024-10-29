import os

def read_homer_motif(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    
    motif_name = lines.strip().split()  # Extract the motif name from the first line
    print(motif_name)
    pwm = []
    for line in lines[1:]:
        if line.strip():  # Ensure the line is not empty
            pwm.append([float(x) for x in line.strip().split()])  # Convert each line to a list of floats
    
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
        f.write(f"MOTIF {motif_name}\n")
        f.write("letter-probability matrix: alength= 4 w= {} nsites= 20 E= 0\n".format(len(pwm)))
        for row in pwm:
            f.write(" ".join(f"{x:.6f}" for x in row) + "\n")
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

