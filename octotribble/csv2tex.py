#!/usr/bin/env python

"""Convert a CSV table into a latex tabular using pandas."""
import argparse
import pandas as pd


def convert(filename, transpose=False):
    """convert csv to tex table."""
    df = pd.read_csv(filename)

    if transpose:
        df = df.transpose()
        tex_name = filename.replace(".csv", "_transpose.tex")
        index = True
    else:
        tex_name = filename.replace(".csv", ".tex")
        index=False

    assert tex_name != filename, "This will overwrite the file, did you pass in a csv?"

    
    latex = df.to_latex(na_rep="-", index=index)

    with open(tex_name, "w") as f:
        f.write(r"\begin{table}")
        f.write("\n")
        f.write(r"\label{}")
        f.write("\n")
        f.write(r"\caption{}")
        f.write("\n")

        f.write(latex)
        f.write(r"\end{table}")
        f.write("\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert csv to latex tabular")
    parser.add_argument("filename", help="Name of csv file", type=str)
    parser.add_argument(
        "-t", "--transpose", help="Transpose table", action="store_true"
    )
    args = parser.parse_args()
    convert(args.filename, args.transpose)
