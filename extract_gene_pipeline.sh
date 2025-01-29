#!/bin/bash

cnsg_script_dir=$(realpath "$1")/
date_format="%d-%m-%Y"
today=$(date +"$date_format")

gene=${2:-"16S"}
taxon=${3:-"eubacteria"}
filters=${4:-"$gene,ssu,ssr"}
feature_type=${5:-"rrna"}


"$cnsg_script_dir"summary_download.sh -i "$taxon" -o "$taxon" -p GCF -l 10
uv run filter.py "$taxon"/"$taxon"_"$today".tsv "$taxon"/"$taxon"_filtered_"$today".tsv
"$cnsg_script_dir"tsv_datasets_downloader.sh -i "$taxon"/"$taxon"_filtered_"$today".tsv -o "$taxon"/ --annotate=true
uv run extract_gene.py "$taxon"/GENOMIC1/ "$taxon"/GFF1/ "$gene"/ "$gene" "$filters" "$feature_type"
