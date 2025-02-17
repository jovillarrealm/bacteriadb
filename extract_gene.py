from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os
from collections import defaultdict
from fssw import sliding_window



def take_input(min_len=None):
    import sys

    err_msg = """Usage:
    extract_gene.py <input_fna_dir> <input_gff_dir> <output_dir> <gene-list (comma separated)> <filter-list (comma separated)>  <feature type>"""
    if min_len and min_len > len(sys.argv) - 1:
        print(err_msg)
        sys.exit()
    _, input_fna_dir, input_gff_dir, output_dir, genes,  feature_type, MAX_LEN = (
        sys.argv
    )
    genes = genes.lower().split(",")
    MAX_LEN=int(MAX_LEN)
    return input_fna_dir, input_gff_dir, output_dir, genes,  feature_type, MAX_LEN


def build_record_id(fna_file: str) -> str:
    GCX, code, genusspecies, strain = fna_file.split("_")
    strain = strain.rstrip(".")
    return f"{GCX}_{code} {genusspecies} {strain}"


def extract_gene(
    gff_file: str,
    fna_file: str,
    output_file: str,
    filename: str,
    genes: list[str],
    search_feature_type: str,
    MAX_LEN:int,
):
    """
    Extract regions containing all specified genes exactly once.
    """
    
    features = []

    # Parse GFF file to collect relevant features
    with open(gff_file, "r") as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            cols = line.lower().strip().split("\t")
            if len(cols) < 9:
                continue
            seqid, source, feat_type, start_str, end_str, _, strand, _, attributes = cols
            if feat_type != search_feature_type:
                continue
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue

            attr_dict = parse_attributes(attributes)
            if "product" not in attr_dict:
                continue

            # Check if the feature matches any of the target genes
            for gene in genes:
                if gene in attr_dict["product"]:
                    features.append((seqid, start, end, strand, gene))

    # Group features by (seqid, strand) and sort by start position
    groups = defaultdict(list)
    for seqid, start, end, strand, gene in features:
        groups[(seqid, strand)].append((start, end, gene))

    # Load FNA sequences
    seq_dict = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))
    records = []
    required_genes = set(genes)

    # Process each group to find valid clusters
    for (seqid, strand), group_features in groups.items():
        # Sort features by start position
        sorted_features = sorted(group_features, key=lambda x: x[0])
        num_required = len(required_genes)
        
        
        # Extract valid clusters from sliding window
        clusters = []
        for chunk in sliding_window(sorted_features, num_required):
            chunk_genes = [gene for _, _, gene in chunk]
            
            # Check if chunk contains all required genes exactly once
            if chunk_genes == genes or chunk_genes == list(reversed(genes)):
                starts = [start for start, _, _ in chunk]
                ends = [end for _, end, _ in chunk]
                min_start = min(starts)
                max_end = max(ends)
                if abs(max_end-min_start) < MAX_LEN:
                    clusters.append((min_start, max_end))

        # Extract sequences for each valid cluster
        for min_start, max_end in clusters:
            seqid_upper = seqid.upper()
            if seqid_upper not in seq_dict:
                print(f"Sequence ID '{seqid_upper}' not found. Skipping.")
                continue
            
            seq_record = seq_dict[seqid_upper]
            sub_seq = seq_record.seq[min_start - 1:max_end]
            if strand == "-":
                sub_seq = sub_seq.reverse_complement()
            
            record_id = build_record_id(filename)
            description = f"{'-'.join(genes)} {strand}"
            record = SeqRecord(sub_seq, id=record_id, description=description)
            records.append(record)
            
            if len(sub_seq) > MAX_LEN:
                print(f"Sequence too long ({len(sub_seq)}).")

    # Write output
    if records:
        gene_part = "-".join(genes).upper()
        output_path = f"{output_file}_{gene_part}.fna"
        SeqIO.write(records, output_path, "fasta")
    else:
        print(f"No valid regions found in {gff_file}.")


def parse_attributes(attributes):
    attr_dict = {}
    for attr in attributes.split(";"):
        if "=" in attr:
            key, val = attr.split("=", 1)
            attr_dict[key.strip()] = val.strip()
    return attr_dict


def main():
    input_fna_dir, input_gff_dir, output_dir, genes,  feature_type, MAX_LEN = take_input(5)
    files = os.listdir(input_gff_dir)
    os.makedirs(output_dir, exist_ok=True)
    for file in files:
        if file.endswith(".gff"):
            filename = os.path.splitext(file)[0]
            fna_path = os.path.join(input_fna_dir, f"{filename}.fna")
            gff_path = os.path.join(input_gff_dir, file)
            output_path = os.path.join(output_dir, filename)
            extract_gene(
                gff_path,
                fna_path,
                output_path,
                f"{filename}.fna",
                genes,
                feature_type,
                MAX_LEN
            )


if __name__ == "__main__":
    main()