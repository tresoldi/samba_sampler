#!/usr/bin/env python

"""
Interface with the samba_sampler library.
"""

# Import Python standard libraries
import argparse
import sys

# Import our library to leverage functions and classes
import samba_sampler as samba


def main():
    parser = argparse.ArgumentParser(
        description="Interface with the samba_sampler library"
    )

    parser.add_argument(
        "--algorithm",
        default="standard",
        help='Algorithm to use: "standard" or "progressive"',
    )
    parser.add_argument("k", type=int, help="Size of samples")
    parser.add_argument("n", type=int, help="Number of samples")
    parser.add_argument(
        "--freq_weight", type=float, default=1.0, help="A frequency weight factor"
    )
    parser.add_argument(
        "--matrices",
        default=None,
        help="List of one or more filenames separated by commas",
    )
    parser.add_argument(
        "--matrix_weights",
        default=None,
        help="List of matrix weights, given as floating points separated by commas",
    )
    parser.add_argument(
        "--tables",
        default=None,
        help="List of one or more filenames separated by commas",
    )
    parser.add_argument(
        "--table_weights",
        default=None,
        help="List of table weights, given as floating points separated by commas",
    )
    args = parser.parse_args()

    # Convert string arguments into corresponding Python types
    matrices = args.matrices.split(",") if args.matrices else None
    tables = args.tables.split(",") if args.tables else None
    matrix_weights = (
        list(map(float, args.matrix_weights.split(",")))
        if args.matrix_weights
        else None
    )
    table_weights = (
        list(map(float, args.table_weights.split(","))) if args.table_weights else None
    )

    # If matrices or tables were provided, make sure the corresponding weights have the same length
    if matrices and matrix_weights and len(matrices) != len(matrix_weights):
        parser.error("The number of matrices and matrix weights must be the same")
    if tables and table_weights and len(tables) != len(table_weights):
        parser.error("The number of tables and table weights must be the same")

    # Create an instance of the sampler
    sampler = samba.GenericSampler(
        matrix_files=matrices,
        table_files=tables,
        matrix_weights=matrix_weights,
        table_weights=table_weights,
    )

    # Print tuples yielded by method sample
    for taxa in sampler.sample(
        args.k, args.n, algorithm=args.algorithm, freq_weight=args.freq_weight
    ):
        sys.stdout.write(str(",".join(sorted(taxa))) + "\n")


if __name__ == "__main__":
    main()
