"""
Common functions for the library.
"""

# Import standard libraries
from collections import defaultdict
from pathlib import Path
from typing import *
import csv
import gzip
import math
import re

ETC_PATH = Path(__file__).parent / "etc"

# TODO: have a single matrix type


def read_splitstree_matrix(filename: Union[Path, str]) -> Dict[str, Dict[str, float]]:
    """
    Read a distance matrix in the SplitsTree format from a file.

    The distance matrix is returned as a dictionary of dictionaries,
    with each value as a dictionary to all other taxa. The function takes care

    Parameters
    ----------
    filename
        The file to read.

    Returns
    -------
    matrix
        A dictionary of dictionaries, where the first key is the taxon and the
        second key is the taxon to which the distance is computed.
    """

    # Read raw data
    header = True
    taxa = []
    matrix = {}
    with open(Path(filename), encoding="utf-8") as handler:
        for line in handler.readlines():
            if header:
                header = False
            else:
                line = re.sub(r"\s+", " ", line)
                tokens = line.split()
                taxon = tokens[0]
                taxa.append(taxon)
                dists = [float(dist) for dist in tokens[1:]]
                matrix[taxon] = dists

    # Make an actual dictionary matrix
    ret_matrix = {}
    for taxon_a, dists in matrix.items():
        ret_matrix[taxon_a] = {taxon_b: dist for dist, taxon_b in zip(dists, taxa)}

    return ret_matrix


def read_triangle_matrix(filename: Union[Path, str]) -> Dict[str, Dict[str, float]]:
    """
    Read a distance matrix in the triangle format from a file.

    The distance matrix is returned as a dictionary of dictionaries,
    with each value as a dictionary to all other taxa. The function takes care
    of opening gzipped files.

    Parameters
    ----------
    filename
        The file to read.

    Returns
    -------
    matrix
        A dictionary of dictionaries, where the first key is the taxon and the
        second key is the taxon to which the distance is computed.
    """

    # Make sure `filename` is a Path object, and open the file with the
    # gzip module if necessary
    filename = Path(filename)
    if filename.suffix == ".gz":
        handler = gzip.open(filename, "rt", encoding="utf-8")
    else:
        handler = open(filename, encoding="utf-8")

    # Read raw data, filling the matrix
    reader = csv.reader(handler, delimiter="\t")
    taxa = next(reader)[1:]
    matrix = defaultdict(dict)
    for line in reader:
        taxon_a = line[0]
        dists = [int(dist) if dist else None for dist in line[1:]]
        for taxon_b, dist in zip(taxa, dists):
            matrix[taxon_a][taxon_b] = dist
            matrix[taxon_b][taxon_a] = dist

    # Close the file when done
    handler.close()

    return matrix


def read_default_matrix() -> Dict[str, Dict[str, float]]:
    """
    Read the default global distance matrix.

    Returns
    -------
    matrix
        A dictionary of dictionaries, where the first key is the taxon and the
        second key is the taxon to which the distance is computed.
    """

    # TODO: read the latest version
    return read_splitstree_matrix(ETC_PATH / "gled_global.dst")
