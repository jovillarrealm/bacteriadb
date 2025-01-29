import polars as pl
from polars import DataFrame
import re
# Read the data using Polars
def take_input(min_len=None):
    import sys

    err_msg = f"""Usage:
    {sys.argv[0]} <input_file> <output_file>"""
    if min_len and min_len > len(sys.argv) - 1:
        print(err_msg)
        sys.exit()
    _, input_file, output_file = sys.argv
    return input_file, output_file

def preliminary_bs(df:DataFrame):
    # Print the data types
    print(df.dtypes)

    # Sort by "Organism Name"
    df = df.sort("Organism Name")

    # Print the first 5 rows
    print(df)

    # Group by "Organism Name" and count occurrences
    sizes = df.group_by("Organism Name2").len()

    # Sort by count in descending order and print the top 12
    print(sizes.sort("len", descending=True))




def filter_and_select_best(df:DataFrame):
    columns="Assembly Accession	Organism Name	Organism Infraspecific Names Strain	Assembly Stats Total Sequence Length	Assembly Stats Contig N50	Assembly Stats GC Count	Assembly Stats GC Percent".split("\t")
    # Calculate the maximum total sequence length within each organism group

    df_with_max = df.with_columns(
        pl.max("Assembly Stats Total Sequence Length").over("Organism Name2").alias("max_length")
    )
    
    # Filter rows where the length is greater than 0.9 times the max in the group
    filtered = df_with_max.filter(
        pl.col("Assembly Stats Total Sequence Length") > 0.9 * pl.col("max_length")
    )
    
    # Sort by Contig N50 descending and select the top row per organism group
    result = filtered.sort("Assembly Stats Contig N50", descending=True).group_by(
        "Organism Name2", maintain_order=True
    ).first()
    
    # Drop the temporary max_length column
    #result = result.drop("max_length").select(columns)
    preliminary_bs(result)

    result = result.drop("Organism Name2").select(columns)
    
    
    return result

# Define a function to clean and extract the first two words
def first_two_words(name: str) -> str:
    if name is None:
        return ""
    # Remove non-alphanumeric characters except spaces
    cleaned = re.sub(r"[^\w\s]", "", name)
    # Extract first two words
    words = cleaned.split()
    return " ".join(words[:2])

if __name__ == "__main__":
    input_file, output_file=take_input(2)
    df = pl.read_csv(input_file, separator="\t")
    df = df.with_columns(
        df["Organism Name"].map_elements(first_two_words, return_dtype=pl.String).alias("Organism Name2")

    )
    preliminary_bs(df)
    df2 = filter_and_select_best(df)
    df2.write_csv(output_file, separator="\t")
    #df2.select(pl.col("Organism Name")).write_csv(reference_input, separator="\t")
            