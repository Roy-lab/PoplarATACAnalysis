#!/bin/bash

# Input files
FASTA_FILE="ExampleData/genome.fa"
TXT_FILE="ExampleData/fragments_corrected_dedup_count.tsv.gz"
OUTPUT_FILE="results/filtered_fragments_corrected_dedup_count.tsv.gz"
REMOVED_COUNT_FILE="removed_count.txt"

# Step 1: Create a chromosome length map from the FASTA file
declare -A chrom_lengths

awk '
  /^>/ {if (seq) print chrom, length(seq); chrom=substr($0,2); seq=""; next}
  {seq = seq $0}
  END {if (seq) print chrom, length(seq)}
' "$FASTA_FILE" > genome_lengths.txt

# Step 2: Filter the TSV file while reading it on the fly using zcat
removed_lines=0

zcat "$TXT_FILE" | awk -v OFS="\t" -v length_file="genome_lengths.txt" '
  BEGIN {
    while ((getline < length_file) > 0) {
      chrom_lengths[$1] = $2
    }
  }
  {
    if ($3 <= chrom_lengths[$1]) {
      print $0
    } else {
      removed++
    }
  }
  END { print removed > "'$REMOVED_COUNT_FILE'" }
' | gzip > "$OUTPUT_FILE"

# Print total removed lines
echo "Total removed lines: $(cat $REMOVED_COUNT_FILE)"
rm $REMOVED_COUNT_FILE  # Cleanup

echo "Filtered file saved as $OUTPUT_FILE"
