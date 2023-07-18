"""
Common functions for the library.
"""

# Import standard libraries
from collections import defaultdict
from pathlib import Path
from typing import List, Tuple, Union, Optional, Sequence
import array
import bz2
import csv
import functools
import itertools
import logging
import math
import pickle
import re

# Import local modules
from .newick import Node

# Define the path to the 'etc' directory, with supplementary material
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

    def __setitem__(self, key: Sequence[str], value: float):
        """
        Sets the value for a specific pair of keys (i, j) in the distance matrix.

        Note that the method does not check if the keys are valid (in particular, if
        only two keys are provided and if they are in the matrix).

        @param key: A sequence containing the two keys.
        @param value: The value to be set for the pair of keys.
        """
        if key[0] == key[1]:
            return 0  # Diagonal values are assumed to be 0

        self.data[self._get_index(key[0], key[1])] = value

    def __getitem__(self, item: Sequence[str]) -> float:
        """
        Returns the value for a specific pair of keys (i, j) in the distance matrix.

        Note that the diagonal values are assumed to be 0. Also note that the
        method does not check if the keys are valid (in particular, if
        only two keys are provided and if they are in the matrix).

        @param item: A sequence containing the two keys.
        @return: The value for the pair of keys.
        """
        if item[0] == item[1]:
            return 0  # Diagonal values are assumed to be 0

        return self.data[self._get_index(item[0], item[1])]

    def rescale(
        self, scale_range: Tuple[float, float] = (0.0, 1.0), factor: float = 1.0
    ):
        """
        Rescales the distance matrix to a given range and by the given factor.

        @param scale_range: The range to which the values should be rescaled.
        @param factor: The factor by which the values should be multiplied.
        """

        # Calculating the minimum and maximum values
        min_value = min(self.data)
        max_value = max(self.data)

        # Build a new temporary array of the same size of self.data, but
        # always of floating point type, to hold the new values
        temp_array = array.array("f", [0] * len(self.data))

        # Rescaling the values
        for i, value in enumerate(self.data):
            temp_array[i] = (value - min_value) / (
                max_value - min_value
            ) * factor + scale_range[0]

        # Use the new array as the data
        self.data = temp_array


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
        matrix[leaf1.name, leaf2.name] = compute_distance(leaf1, leaf2)

    return matrix


def dst2matrix(filename: Union[Path, str]) -> DistanceMatrix:
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
    mtx = DistanceMatrix(taxa)
    for taxon_a, dists in matrix.items():
        for dist, taxon_b in zip(dists, taxa):
            mtx[taxon_a, taxon_b] = dist

    return mtx


def build_table_from_file(
    filename: Union[str, Path],
    key: str,
    value: str,
    encoding: str = "utf-8",
    multiple: str = "average",
) -> dict:
    """
    Build a dictionary from a tabular file.

    @param filename: The name of the file to read.
    @param key: The name of the column to use as key.
    @param value: The name of the column to use as value.
    @param encoding: The encoding of the file.
    @param multiple: How to handle multiple values for the same key.
    @return: A dictionary with the values from the file.
    """

    # Make sure filename is a Path object
    filename = Path(filename)

    # Check that multiple is a valid option
    if multiple not in ["average", "max", "min"]:
        raise ValueError(f"Invalid value for 'multiple': {multiple}")

    # Open the file
    with filename.open(newline="", encoding=encoding) as csvfile:
        # Let the csv library sniff the dialect of the file
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        csvfile.seek(0)

        # Create a csv reader object
        reader = csv.DictReader(csvfile, dialect=dialect)

        # Create a defaultdict of lists to store all values associated with each key
        result_dict = defaultdict(list)

        for row in reader:
            # Append each value to the list of values for the appropriate key
            result_dict[row[key]].append(float(row[value]))

        # Calculate the final value to be associated with each key, based on the value of 'multiple'
        for key, values in result_dict.items():
            if multiple == "average":
                result_dict[key] = sum(values) / len(values)
            elif multiple == "max":
                result_dict[key] = max(values)
            elif multiple == "min":
                result_dict[key] = min(values)

    return dict(result_dict)


def filter_data(filter_string: str, data: List[dict]) -> List[dict]:
    """
    Apply filtering criteria to a list of dictionaries.

    @param filter_string: String that specifies the filtering criteria. The string format should be
    "feature1=-value1,value2;feature2=value3,-value4;feature3=>value5,<value6", where '-' in front of
    a value indicates that the dictionaries should not contain that value for the feature, '>' indicates
    that the dictionaries should contain a value that includes the given value (substring or list element),
    and '<' indicates that the dictionaries should not contain a value that includes the given value.

    @param data: List of dictionaries to apply the filtering criteria to.
    @return: List of dictionaries that meet all the filtering criteria.
    """

    # Parsing the filter string
    filters = filter_string.split(";")
    filter_dict = {}
    for filt in filters:
        key, values = filt.split("=")
        filter_dict[key] = {"neg": [], "pos": [], "part": [], "not_part": []}
        for value in values.split(","):
            if value.startswith("-"):
                filter_dict[key]["neg"].append(value[1:])
            elif value.startswith(">"):
                filter_dict[key]["part"].append(value[1:])
            elif value.startswith("<"):
                filter_dict[key]["not_part"].append(value[1:])
            else:
                filter_dict[key]["pos"].append(value)

    # Filtering data
    filtered_data = []
    for row in data:
        valid_row = True
        for key, value in filter_dict.items():
            if key in row:
                if row[key] in value["neg"] or (
                    value["pos"] and row[key] not in value["pos"]
                ):
                    valid_row = False
                    break
                for part_value in value["part"]:
                    if isinstance(row[key], str):
                        if part_value not in row[key]:
                            valid_row = False
                            break
                    elif isinstance(row[key], list):
                        if part_value not in row[key]:
                            valid_row = False
                            break
                for not_part_value in value["not_part"]:
                    if isinstance(row[key], str):
                        if not_part_value in row[key]:
                            valid_row = False
                            break
                    elif isinstance(row[key], list):
                        if not_part_value in row[key]:
                            valid_row = False
                            break
        if valid_row:
            filtered_data.append(row)

    return filtered_data


def filter_glottolog(filter_string: str) -> List[str]:
    """
    Reads the glottolog data and filters it according to the filter string.

    The function will load the glottolog data dump and filter it according to the
    filter string (using the `filter_data()` function). The function will then return
    a list of glottocodes that meet the filtering criteria.
    """

    # Obtain the path to the glottolog data dump and read it as a list of dictionaries
    glottolog_path = sorted(ETC_PATH.glob("glottolog.*.tsv"))[-1]
    glottolog_data = list(
        csv.DictReader(glottolog_path.open(encoding="utf-8"), delimiter="\t")
    )

    # Filter the glottolog data
    glottolog_data = filter_data(filter_string, glottolog_data)

    # Return a list of glottocodes
    return [row["glottocode"] for row in glottolog_data]
