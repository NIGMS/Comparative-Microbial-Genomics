#!/bin/bash

# Usage: ./download_sra.sh sra_list.txt /path/to/output_dir

sra_list="$1"
output_dir="$2"

mkdir -p "$output_dir"

while read -r accession; do
    echo "Downloading $accession..."
    fastq-dump --outdir "$output_dir" --gzip --split-files "$accession"
done < "$sra_list"
