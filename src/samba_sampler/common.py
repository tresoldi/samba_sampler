"""
Common functions for the library.
"""

# Import standard libraries
from collections import defaultdict
from pathlib import Path
from typing import *
import array
import bz2
import pickle
from functools import cache
import csv
import gzip
import itertools
import re
import math

ETC_PATH = Path(__file__).parent / "etc"


class DistanceMatrix:
    """
    This class represents a symmetric distance matrix in a memory-efficient way.
    Instead of storing all values, it stores only half of the elements excluding the diagonal,
    as the matrix is symmetric and the diagonal is assumed to be zeros. The matrix elements are
    mapped to a one-dimensional array.
    """

    def __init__(
        self,
        keys: Union[List[str], None] = None,
        filename: Union[str, None] = None,
        datatype: str = "f",
    ):
        """
        Initializes the matrix. The matrix can be initialized in two ways:

        1. By providing a list of keys. An empty matrix will be created with the given keys.
        2. By providing a filename from which to read the matrix.

        Exactly one of 'keys' or 'filename' must be provided.

        :param keys: A list of unique keys representing the matrix elements.
        :param filename: The filename from which to read the matrix.
        :param datatype: The datatype of the array elements. Must be a valid type code.
        """

        # Check that exactly one of keys or filename is provided
        if (keys is None) == (filename is None):
            raise ValueError(
                "Either 'keys' or 'filename' must be provided, but not both."
            )

        # Initialize the matrix
        if filename is not None:
            self.keys, self.data = self.read_matrix(filename)
        else:
            self.keys = sorted(keys)
            self.data = array.array(
                datatype, [0] * (len(self.keys) * (len(self.keys) - 1) // 2)
            )

        # Generate indices for quick lookup
        self.indices = {key: i for i, key in enumerate(self.keys)}

    def read_matrix(self, filename: Union[str, Path]) -> Tuple[List[str], array.array]:
        """
        Reads a matrix from a compressed file.

        This method reads the compressed file, deserializes it and returns the keys and data.

        :param filename: The filename or Path from which to read the matrix.
        :return: The keys and data as a tuple.
        """

        # Convert to pathlib.Path if necessary
        if isinstance(filename, str):
            filename = Path(filename)

        with bz2.open(filename, "rb") as file:
            keys, data = pickle.load(file)

        return keys, data

    def save_matrix(self, filename: Union[str, Path]):
        """
        Writes the matrix to a compressed file.

        This method serializes the keys and data and writes them to a compressed file.

        :param filename: The filename or Path to which to write the matrix.
        """

        # Convert to pathlib.Path if necessary
        if isinstance(filename, str):
            filename = Path(filename)

        with bz2.open(filename, "wb") as file:
            pickle.dump((self.keys, self.data), file)

    @cache
    def _get_index(self, i: str, j: str) -> int:
        """
        Computes the index in the array for the keys i and j.

        This method uses the fact that the distance matrix is stored in the form of a one-dimensional
        array. It calculates the index of each element in the array based on the keys. The calculation
        is based on the lower left triangle of the matrix, excluding the diagonal.
        """

        if self.indices[i] > self.indices[j]:
            i, j = j, i  # Ensure i <= j

        return self.indices[j] * (self.indices[j] - 1) // 2 + self.indices[i]

    def set(self, i: str, j: str, value: Union[int, float]):
        """
        Sets the value at keys i and j in the distance matrix.

        :param i: The first key.
        :param j: The second key.
        :param value: The value to be set.
        """
        self.data[self._get_index(i, j)] = value

    def get(self, i: str, j: str) -> Union[int, float]:
        """
        Retrieves the value at keys i and j from the distance matrix.

        This method is symmetric, i.e., get(i, j) == get(j, i). If i == j, the return value is 0,
        since the diagonal of the distance matrix is assumed to be 0.

        ::param i: The first key.
        :param j: The second key.
        :return: The value at keys i and j.
        """

        if i == j:
            return 0  # Diagonal values are assumed to be 0

        return self.data[self._get_index(i, j)]


def tree2matrix(tree):
    """
    Convert a newick tree to a distance matrix.

    :param tree:
    :return:
    """

    def most_recent_common_ancestor(anc_list1, anc_list2):
        for anc in anc_list1:
            if anc in anc_list2:
                return anc

        raise ValueError("No common ancestor found")

    @cache
    def compute_distance(leaf1, leaf2):
        # Get the most recent common ancestor of the two leaves
        mrca = most_recent_common_ancestor(ancestors[leaf1], ancestors[leaf2])

        # Sum the lengths between leaves and the mrca, grabbing the nodes by label
        # with the `get_node` method
        leaf1_length = sum(
            [n.length for n in ancestors[leaf1][: ancestors[leaf1].index(mrca)]]
        )
        leaf2_length = sum(
            [n.length for n in ancestors[leaf2][: ancestors[leaf2].index(mrca)]]
        )

        return leaf1_length + leaf2_length

    # Get all the leaves of the tree
    leaves = tree.get_leaves()

    # Build a dictionary with the leaves as keys and their ancestors as values
    ancestors = {leaf: leaf.ancestors for leaf in leaves}

    # For each pairwise combination of leaves, get their most recent common ancestor
    distances = {}
    num_comb = math.comb(len(leaves), 2)
    for idx, (leaf1, leaf2) in enumerate(itertools.combinations(leaves, 2)):
        if idx % 1000 == 0:
            print(
                f"Processed {idx} pairs of leaves (at `{leaf1.name},{leaf2.name}`) [{(idx/num_comb)*100:.2f}%]..."
            )

        #if idx == 10000:
        #    break

        distances[leaf1.name, leaf2.name] = compute_distance(leaf1, leaf2)

    return distances


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
