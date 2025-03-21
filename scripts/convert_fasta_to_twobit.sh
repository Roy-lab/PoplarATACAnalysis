#!/bin/bash

# Define the path to the faToTwoBit executable
FA_TO_2BIT="/mnt/dv/wid/projects5/Roy-LongRangeEvolution/ENSEMBL_Release_102_multiple_alignments/hal_chain_files/faToTwoBit"

# Define the directory containing your FASTA files
SEQ_SRC_DIR="/mnt/dv//wid/projects7/Roy-plants/kirstlab/scripts/Poplar_ArchR"

# Change to the directory
cd "$SEQ_SRC_DIR"

# Loop through each .fa file and convert it to .2bit
#for fa_file in *.fa; do
#    twobit_file="${fa_file%.fa}.2bit"  # Change extension from .fa to .2bit
#    "$FA_TO_2BIT" "$fa_file" "$twobit_file"
#done

#echo "Conversion complete!"


# Define the input and output file names
FASTA_FILE="filtered_genome.fa"
TWOBIT_FILE="filtered_genome.2bit"

# Convert the filtered FASTA file to 2bit format
"$FA_TO_2BIT" "$FASTA_FILE" "$TWOBIT_FILE"

echo "Conversion complete! Output: $TWOBIT_FILE"
