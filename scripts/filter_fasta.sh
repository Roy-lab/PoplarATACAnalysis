#!/bin/bash

# Define the input and output file paths
input_file="genome.fa"
output_file="filtered_genome.fa"

# Process the file: remove scaffold sequences and rename chromosome headers
awk '
  /^>/ {
    if ($0 !~ /scaffold/) {
      chr_num = substr($0, 5) + 0  # Convert Chr01 -> 1, Chr10 -> 10 (removes leading zeros)
      print ">chr" chr_num
    }
    next
  }
  { print }
' "$input_file" > "$output_file"

echo "Conversion complete! Output saved in $output_file"
