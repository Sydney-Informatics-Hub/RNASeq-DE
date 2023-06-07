import os

# Directory path containing the quant.sf files
quant_dir = "../Batch_1_salmon"

# Search for quant.sf files within the directory and its subdirectories
quant_files = []
for root, dirs, files in os.walk(quant_dir):
    for file in files:
        if file == "quant.sf":
            quant_files.append(os.path.join(root, file))

# Output file for the combined matrix
output_file = "../Batch_1_salmon/combined_matrix.txt"

# Extract the common set of genes from the first quant.sf file
with open(quant_files[0], "r") as file:
    common_genes = [line.split("\t")[0] for line in file.readlines()[1:]]

# Create the header for the matrix using the quant.sf file names
header = "\t".join(quant_files)

# Write the header to the output file
with open(output_file, "w") as file:
    file.write(f"Gene\t{header}\n")

# Iterate through each common gene
for gene in common_genes:
    row = gene
    
    # Iterate through each quant.sf file
    for file in quant_files:
        with open(file, "r") as quant_file:
            expression_value = ""
            for line in quant_file:
                if line.split("\t")[0] == gene:
                    expression_value = line.split("\t")[-1].strip()
                    break

        row += f"\t{expression_value}"
    
    # Write the row to the output file
    with open(output_file, "a") as file:
        file.write(f"{row}\n")

print(f"Combined matrix is generated in {output_file}")

