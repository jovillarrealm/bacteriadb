import polars as pl
import re
import sys


def take_input(min_len=None) -> tuple[str, str]:
    err_msg = f"""Usage:
    {sys.argv[0]} <input_file> <output_file>"""
    if min_len and min_len > len(sys.argv) - 1:
        print(err_msg)
        sys.exit()
    _, input_file, output_file = sys.argv
    return input_file, output_file


def preliminary_bs(df: pl.DataFrame):
    print(df.dtypes)
    df = df.sort("Organism Name")
    print(df)
    sizes = df.group_by("Organism Name2").len()
    print(sizes.sort("len", descending=True))


def filter_and_select_best(df: pl.DataFrame, columns: list[str]) -> pl.DataFrame:
    # Use lazy evaluation for efficiency
    lf = df.lazy()

    lf = lf.with_columns(
        pl.max("Assembly Stats Total Sequence Length")
        .over("Organism Name2")
        .alias("max_length")
    )

    lf = lf.filter(
        pl.col("Assembly Stats Total Sequence Length") > 0.9 * pl.col("max_length")
    )

    lf = (
        lf.sort("Assembly Stats Contig N50", descending=True)
        .group_by("Organism Name2", maintain_order=True)
        .first()
    )

    result = lf.collect()  # Collect the lazyframe when needed
    preliminary_bs(result)
    result = result.drop("Organism Name2").select(columns)
    return result


def first_two_words(name: str) -> str:
    if name is None:
        return ""
    cleaned = re.sub(r"[^\w\s]", "", name)
    words = cleaned.split()
    return " ".join(words[:2])


if __name__ == "__main__":
    input_file, output_file = take_input(2)

    # Use lazy reading for initial CSV parsing
    columns = "Assembly Accession	Organism Name	Organism Infraspecific Names Strain	Assembly Stats Total Sequence Length	Assembly Stats Contig N50	Assembly Stats GC Count	Assembly Stats GC Percent".split(
        "\t"
    )  # Simplified split

    df = pl.scan_csv(input_file, separator="\t", with_column_names=columns)

    df = df.with_columns(
        df["Organism Name"]
        .map_elements(first_two_words, return_dtype=pl.String)
        .alias("Organism Name2")
    ).collect()  # collect here since map_elements is eager

    preliminary_bs(df)
    df2 = filter_and_select_best(df, columns)
    df2.write_csv(output_file, separator="\t")
