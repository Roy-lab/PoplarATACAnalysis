zcat filtered_fragments_corrected_dedup_count.tsv.gz | awk '
{
    chr_num = substr($1, 4) + 0  # Extract number and remove leading zeros
    print "chr" chr_num, $2, $3, $4, $5
}' OFS="\t" | gzip > renamed_filtered_fragments_corrected_dedup_count.tsv.gz
