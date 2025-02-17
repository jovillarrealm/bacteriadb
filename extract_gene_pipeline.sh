#!/bin/bash
: "${1:?Usage: $0 cnsg_script_dir <gene> <taxon> <filters> <feature_type> <limit>}"
script_dir=$(dirname "$0")/
cnsg_script_dir=$(realpath "$1")/
date_format="%d-%m-%Y"
today=$(date +"$date_format")
genes=${2:-"16S,23S"}
taxon=${3:-"eubacteria"}
feature_type=${4:-"rrna"}
max_len=${6:"8000"}

output_dir=$(echo "$genes" | awk 'BEGIN { FS="," ; OFS="-" } { gsub(",", OFS, $0); print }')

"$cnsg_script_dir"summary_download.sh -i "$taxon" -o "$taxon" -p GCF ${5:+ -l "$5"}
uv run --project "$script_dir" "$script_dir"filter.py "$taxon"/"$taxon"_"$today".tsv "$taxon"/"$taxon"_filtered_"$today".tsv
"$cnsg_script_dir"tsv_datasets_downloader.sh -i "$taxon"/"$taxon"_filtered_"$today".tsv -o "$taxon"/ --annotate=true
uv run --project "$script_dir" "$script_dir"extract_gene.py "$taxon"/GENOMIC1/ "$taxon"/GFF1/ "$output_dir"/ "$genes" "$feature_type" "$max_len"
