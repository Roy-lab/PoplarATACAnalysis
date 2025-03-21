import re

# Define the input and output file paths
input_file = "ExampleData/gene.gtf"
output_file = "results/renamed_gene.gtf"

# Open the input and output files
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        # Skip comment lines (those starting with #)
        if line.startswith('#'):
            outfile.write(line)
            continue
        
        # Split the line into fields
        fields = line.strip().split('\t')
        
        # Rename chromosomes if necessary (matching patterns like Chr1, Chr2, ..., Chr19)
        if fields[0].startswith('Chr'):
            chrom = fields[0][3:]  # Remove 'Chr'
            # Convert the chromosome number to avoid leading zeros, then prepend 'chr'
            chrom = 'chr' + str(int(chrom))  # This will remove any leading zeros
            fields[0] = chrom  # Update the chromosome name
        
        # Join the fields back and write to the output file
        outfile.write('\t'.join(fields) + '\n')

print("Renaming complete, no missing fields!")
