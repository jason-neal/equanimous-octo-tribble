#!/usr/bin/env python

"""Convert a CSV table into a latex tabular using pandas."""
import argparse
import pandas as pd


def convert(filename, transpose=False):
    """convert csv to tex table."""
    df = pd.from_csv(filename)

    if transpose:
        df = df.transpose()
        tex_name = filename.replace(".csv", "_transpose.tex")
    else:
        tex_name = filename.replace(".csv", ".tex")
    with open(tex_name, "w") as f:
        f.write(r"\begin{table}")
        f.write(r"\label{}")
        f.write(r"\caption{}")

        f.write(tex_name.df.to_latex(na_rep="-"))
        f.write(r"\end{table}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert csv to latex tabular")
    parser.add_argument("filename", help="Name of csv file", type=str)
    parser.add_argument(
        "-t", "--transpose", help="Transpose table", action="store_true"
    )
    args = parser.parse_args()
    convert(args.filename, args.transpose)
