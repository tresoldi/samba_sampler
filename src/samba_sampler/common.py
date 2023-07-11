"""
Common functions for the library.
"""

# Import standard libraries
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Union, Dict, Optional
import array
import bz2
import csv
import functools
import gzip
import itertools
import logging
import math
import pickle
import re

from .newick import Node

ETC_PATH = Path(__file__).parent / "etc"


class DistanceMatrix:
    """
    This class provides a memory-efficient implementation of a symmetric distance matrix.
    Rather than storing all elements, it only keeps half of the elements (excluding the diagonal),
    leveraging the symmetry of the matrix. The diagonal elements are assumed to be zeros.
    The matrix elements are flattened into a one-dimensional array.
    """

    data: array.array

    def __init__(
        self,
        keys: Optional[List[str]] = None,
        filename: Optional[Union[str, Path]] = None,
        datatype: str = "f",
    ):
        """
        Initializes the matrix. The matrix can be populated in one of two ways:

        1. By providing a list of keys. An empty matrix is initialized with these keys.
        2. By providing a filename to load the matrix from.

        Exactly one of 'keys' or 'filename' must be provided.

        @param keys: A list of unique keys that represent the elements of the matrix.
        @param filename: The name (or Path object) of the file from which to load the matrix.
        @param datatype: The datatype of the elements in the array. Must be a valid type code.
        """

        # Ensuring that either keys or filename is provided
        if (keys is None) == (filename is None):
            raise ValueError(
                "Either 'keys' or 'filename' must be provided, but not both."
            )

        # If filename is provided, read the matrix from file
        if filename:
            self.keys, self.data = self.read(filename)
        else:  # Else, create an empty matrix with given keys
            self.keys = sorted(keys)
            # Initializing the half-matrix with zeroes
            self.data = array.array(
                datatype, [0] * (len(self.keys) * (len(self.keys) - 1) // 2)
            )

        # Generate indices for each key for quick lookup
        self.indices = {key: i for i, key in enumerate(self.keys)}

    def read(self, filename: Union[str, Path]) -> Tuple[List[str], array.array]:
        """
        Reads a matrix from a compressed file.

        The method opens the compressed file, deserializes it, and returns the keys and data.

        @param filename: The filename or Path object to read the matrix from.
        @return: A tuple containing the keys and data of the matrix.
        """

        # Ensuring filename is a Path object
        filename = Path(filename) if isinstance(filename, str) else filename

        with bz2.open(filename, "rb") as file:
            keys, data = pickle.load(file)

        return keys, data

    def save(self, filename: Union[str, Path]):
        """
        Writes the matrix to a compressed file.

        This method serializes the keys and data and writes them to a compressed file.

        @param filename: The filename or Path object to write the matrix to.
        """

        # Ensuring filename is a Path object
        filename = Path(filename) if isinstance(filename, str) else filename

        with bz2.open(filename, "wb") as file:
            pickle.dump((self.keys, self.data), file)

    @functools.cache  # Caching results for repeated calls with the same arguments
    def _get_index(self, i: str, j: str) -> int:
        """
        Computes the index in the flattened array for the given pair of keys (i, j).

        This method calculates the index of each element in the one-dimensional array based on the keys.
        The computation is performed based on the lower-left triangle of the matrix,
        excluding the diagonal.

        @param i: The first key.
        @param j: The second key.
        @return: The index in the array corresponding to the pair of keys.
        """
        if self.indices[i] > self.indices[j]:
            i, j = j, i  # Ensuring i <= j for the calculation

        return self.indices[j] * (self.indices[j] - 1) // 2 + self.indices[i]

    def set(self, i: str, j: str, value: Union[int, float]):
        """
        Sets the value for a specific pair of keys (i, j) in the distance matrix.

        @param i: The first key.
        @param j: The second key.
        @param value: The value to be set for the pair of keys.
        """
        self.data[self._get_index(i, j)] = value

    def get(self, i: str, j: str) -> Union[int, float]:
        """
        Retrieves the value for a specific pair of keys (i, j) from the distance matrix.

        This method is symmetric, i.e., get(i, j) == get(j, i). If i == j, the returned value is 0,
        as the diagonal of the distance matrix is assumed to be 0.

        @param i: The first key.
        @param j: The second key.
        @return: The value corresponding to the pair of keys (i, j).
        """
        if i == j:
            return 0  # Diagonal values are assumed to be 0

        return self.data[self._get_index(i, j)]

    # TODO: join with get() above
    def __getitem__(self, item):
        """
        Retrieves the value for a specific pair of keys (i, j) from the distance matrix.

        This method is symmetric, i.e., get(i, j) == get(j, i). If i == j, the returned value is 0,
        as the diagonal of the distance matrix is assumed to be 0.

        @param item: A tuple containing the two keys.
        @return: The value corresponding to the pair of keys (i, j).
        """
        return self.get(*item)

    def rescale(self, scale_range=(0, 1), factor: float = 1.0):
        """
        Rescales the distance matrix to a given range and by the given factor.
        """

        # Calculating the minimum and maximum values
        min_value = min(self.data)
        max_value = max(self.data)

        # Rescaling the values
        for i, value in enumerate(self.data):
            self.data[i] = (value - min_value) / (
                max_value - min_value
            ) * factor + scale_range[0]


def tree2matrix(tree: Node) -> DistanceMatrix:
    """
    Converts a Newick tree into a symmetric distance matrix.

    @param tree: The input Newick tree to be converted.
    @return: The resulting distance matrix.
    """

    def most_recent_common_ancestor(
        anc_list1: List[Node], anc_list2: List[Node]
    ) -> Node:
        """
        Finds the most recent common ancestor of two lists of ancestors.

        @param anc_list1: The first list of ancestors.
        @param anc_list2: The second list of ancestors.
        @return: The label of the most recent common ancestor.
        @raise ValueError: Raised when no common ancestor is found.
        """
        for anc in anc_list1:
            if anc in anc_list2:
                return anc
        raise ValueError("No common ancestor found")

    @functools.cache
    def compute_distance(leaf1: str, leaf2: str) -> float:
        """
        Computes the distance between two leaves in the tree.

        The distance is computed as the sum of lengths from each leaf to their most recent common ancestor (MRCA).

        @param leaf1: The first leaf.
        @param leaf2: The second leaf.
        @return: The computed distance.
        """
        # Get the most recent common ancestor of the two leaves
        mrca = most_recent_common_ancestor(ancestors[leaf1], ancestors[leaf2])

        # Compute the lengths between leaves and the MRCA
        leaf1_length = sum(
            [n.length for n in ancestors[leaf1][: ancestors[leaf1].index(mrca)]]
        )
        leaf2_length = sum(
            [n.length for n in ancestors[leaf2][: ancestors[leaf2].index(mrca)]]
        )

        return leaf1_length + leaf2_length

    # Extract all leaves from the tree
    leaves = tree.get_leaves()

    # Initialize the distance matrix
    matrix = DistanceMatrix([leaf.name for leaf in leaves])

    # Build a dictionary mapping leaves to their ancestors; note that this currently
    # requires a complete traversal of the tree for each leaf, which is not efficient
    # (this could be improved by storing the ancestors in the tree nodes, but involves
    # changes to the `Node` class that are not urgent)
    ancestors = {leaf: leaf.ancestors for leaf in leaves}

    # Compute pairwise distances for each combination of leaves
    num_comb = math.comb(len(leaves), 2)
    for idx, (leaf1, leaf2) in enumerate(itertools.combinations(leaves, 2)):
        if idx % 1000 == 0:
            logging.info(
                f"Processed {idx} pairs of leaves (at `{leaf1.name},{leaf2.name}`) [{(idx/num_comb)*100:.2f}%]..."
            )
        matrix.set(leaf1.name, leaf2.name, compute_distance(leaf1, leaf2))

    return matrix


##########################


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
