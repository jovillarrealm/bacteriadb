#!/bin/bash
: "${1:?cnsg_script_dir is unset or null, pls give path to cnsg-downloader-code}"
script_dir=$(dirname "$0")/
cnsg_script_dir=$(realpath "$1")/
date_format="%d-%m-%Y"
today=$(date +"$date_format")

gene=${2:-"16S"}
taxon=${3:-"eubacteria"}
filters=${4:-"$gene,ssu,ssr"}
feature_type=${5:-"rrna"}


"$cnsg_script_dir"summary_download.sh -i "$taxon" -o "$taxon" -p GCF ${6:+ -l "$6"}
uv run "$script_dir"filter.py "$taxon"/"$taxon"_"$today".tsv "$taxon"/"$taxon"_filtered_"$today".tsv
"$cnsg_script_dir"tsv_datasets_downloader.sh -i "$taxon"/"$taxon"_filtered_"$today".tsv -o "$taxon"/ --annotate=true
uv run "$script_dir"extract_gene.py "$taxon"/GENOMIC1/ "$taxon"/GFF1/ "$gene"/ "$gene" "$filters" "$feature_type"
