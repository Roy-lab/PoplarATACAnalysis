#!/bin/bash

# Define the path to the faToTwoBit executable
FA_TO_2BIT="/mnt/dv/wid/projects5/Roy-LongRangeEvolution/ENSEMBL_Release_102_multiple_alignments/hal_chain_files/faToTwoBit"

# Define the input and output file names
FASTA_FILE="results/filtered_genome.fa"
TWOBIT_FILE="results/filtered_genome.2bit"

# Convert the filtered FASTA file to 2bit format
"$FA_TO_2BIT" "$FASTA_FILE" "$TWOBIT_FILE"

echo "Conversion complete! Output: $TWOBIT_FILE"
