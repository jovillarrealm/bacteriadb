from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import os


# Read the data using Polars
def take_input(min_len=None):
    import sys

    err_msg = """Usage:
    extract_gene.py <input_fna_dir> <input_gff_dir> <output_dir> <gene> <filterlist (comma separated)>  <feature type>"""
    if min_len and min_len > len(sys.argv) - 1:
        print(err_msg)
        sys.exit()
    _, input_fna_dir, input_gff_dir, output_dir, gene, filterlist, feature_type = sys.argv
    filterlist = filterlist.split(",")
    return input_fna_dir, input_gff_dir, output_dir, gene, filterlist, feature_type


def build_record_id(fna_file: str) -> str:
    GCX, code, genusspecies, strain = fna_file.split("_")
    strain = strain.rstrip(".")
    out = f"{GCX}_{code} {genusspecies} {strain}"
    return out


def extract_gene(
    gff_file: str,
    fna_file: str,
    output_file: str,
    filename: str,
    gene: str,
    filters: list[str],
    search_feature_type: str,
):
    """
    Extract gene rRNA sequences from a GFF and FNA file.

    :param gff_file: Path to GFF file.
    :param fna_file: Path to FNA file.
    :param output_file: Path to output file.
    :param filters: List of strings to search for in the product attribute. Case insensitive.
    :param gene: Gene name to search for in the feature type. Case insensitive.
    :feature_type: Feature type to search for in the GFF file. For example, 'rRNA'. Case insensitive"""
    # Parse GFF file for gene features

    features = []
    with open(gff_file, "r") as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if len(cols) < 9:
                continue
            seqid, source, feat_type, start_str, end_str, _, strand, _, attributes = (
                cols
            )
            if feat_type.lower() != search_feature_type.lower():
                continue
            try:
                start = int(start_str)
                end = int(end_str)
            except ValueError:
                continue  # Skip invalid entries
            # Parse attributes to check for gene
            attr_dict = {}
            for attr in attributes.split(";"):
                if "=" in attr:
                    key, val = attr.split("=", 1)
                    attr_dict[key.strip().lower()] = val.strip().lower()
            # Check if product attribute indicates gene
            for attempt in filters:
                if "product" in attr_dict and (attempt.lower() in attr_dict["product"].lower()):
                    features.append((seqid, start, end, strand))

    # Load FNA sequences into a dictionary
    seq_dict = SeqIO.to_dict(SeqIO.parse(fna_file, "fasta"))


    # Extract sequences for each feature
    records = []
    for seqid, start, end, strand in features:
        if seqid not in seq_dict:
            print(f"Sequence ID '{seqid}' not found in FNA. Skipping.")
            continue
        seq_record = seq_dict[seqid]
        sub_seq = seq_record.seq[start - 1 : end]  # Convert to 0-based index
        if strand == "-":
            sub_seq = sub_seq.reverse_complement()
        # Create a new SeqRecord
        record_id = build_record_id(filename)
        record = SeqRecord(
            sub_seq,
            id=record_id,
            description=f"{gene} {strand}",
        )
        records.append(record)

    # Write output
    if records:
        SeqIO.write(records, f"{output_file}_{gene.upper()}.fna", "fasta")
        # print(f"Successfully wrote {len(records)} {gene} sequences to {output_file}")
    else:
        print(f"No {gene} {search_feature_type} features found in {fna_file}.")


def main():
    input_fna_dir, input_gff_dir, output_dir, gene, filters, feature_type = take_input(6)
    files = os.listdir(input_gff_dir)
    os.makedirs(output_dir, exist_ok=True)
    for file in files:
        if file.endswith(".gff"):
            filename, file_extension = os.path.splitext(file)
            fna_file, gff_file = (
                os.path.join(input_fna_dir, filename + ".fna"),
                os.path.join(input_gff_dir, file),
            )

            output_file = os.path.join(output_dir, filename)
            extract_gene(
                gff_file, fna_file, output_file, filename + ".fna", gene, filters, feature_type
            )


if __name__ == "__main__":
    main()
