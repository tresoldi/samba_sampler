#!/usr/bin/env python3

"""
Build supplementary data for the scaffolding.
"""

# Import Python standard libraries
from pathlib import Path
import csv
import glob
import itertools
import logging
import math
import os

# Import 3rd party libraries
from ete3 import Tree

# Import our library to leverage functions and classes
from samba_sampler import newick, tree2matrix, DistanceMatrix

ROOT_PATH = Path(__file__).parent.parent


def haversine(lat1: float, lon1: float, lat2: float, lon2: float) -> float:
    """
    Compute the Haversine distance between two coordinates.

    Parameters
    ----------
    lat1 : float
        The latitude of the first coordinate.
    lon1 : float
        The longitude of the first coordinate.
    lat2 : float
        The latitude of the second coordinate.
    lon2 : float
        The longitude of the second coordinate.

    Returns
    -------
    float
        The Haversine distance between the two coordinates.
    """

    radius = 6372.8  # km

    diff_lat = math.radians(lat2 - lat1)
    diff_lon = math.radians(lon2 - lon1)
    lat1 = math.radians(lat1)
    lat2 = math.radians(lat2)

    a = (
        math.sin(diff_lat / 2) ** 2
        + math.cos(lat1) * math.cos(lat2) * math.sin(diff_lon / 2) ** 2
    )
    c = 2 * math.asin(math.sqrt(a))

    return radius * c


def build_geodistance():
    """
    Build a CSV file with the distance between all pairs of languages.

    This function currently only computes the Haversine distance between
    the coordinates of the languages, but it could be extended to
    compute other distances, such as the distance along a road network,
    """

    # Temporary function specifying whether to accept or not a language
    # based on its Glottolog level and exceptions from GLED
    def _accept(entry):
        if entry["glottocode"] in [
            "alba1267",
            "geji1246",
            "nuuu1241",
            "jiar1240",
            "sout2986",
        ]:
            return True

        return entry["level"] in ["language", "dialect"]

    # Grab the path to the latest dump of the Glottolog database and
    # read the data into a dictionary
    glottolog_path = sorted(
        glob.glob(str(ROOT_PATH / "src" / "samba_sampler" / "etc" / "glottolog.*.tsv"))
    )[-1]
    with open(glottolog_path, encoding="utf-8") as f:
        # Read the data, filter entries whose `level` is not `language` or `dialect`,
        # filter to remove all entries without coordinates and then convert
        # the coordinates to floats
        # NOTE: filtering for languages/dialects is including the problematic items
        #       from ASJP that are propagated to GLED (families listed as languages)
        glottolog = {
            row["glottocode"]: row for row in csv.DictReader(f, delimiter="\t")
        }
        glottolog = {
            k: v
            for k, v in glottolog.items()
            if _accept(v) and v["latitude"] and v["longitude"]
        }
        glottolog = {
            k: {
                **v,
                "latitude": float(v["latitude"]),
                "longitude": float(v["longitude"]),
            }
            for k, v in glottolog.items()
        }

    # Extract the glottocodes
    glottocodes = sorted(glottolog.keys())

    # Compute the number of combinations of languages for reporting
    n = len(glottolog)
    n_combinations = n * (n - 1) / 2

    # Compute the Haversine distance between all pairs of languages and build the matrix
    logging.info("[build_geodistance] Building distance matrix")
    matrix = DistanceMatrix(glottocodes, datatype="h")  # C signed short type
    for idx, (l1, l2) in enumerate(itertools.combinations(glottocodes, 2)):
        if idx % 100000 == 0:
            logging.info(
                "[build_geodistance] Computing distance between {} and {} ({:.2f}%)".format(
                    l1, l2, 100 * idx / n_combinations
                )
            )
        d = haversine(
            glottolog[l1]["latitude"],
            glottolog[l1]["longitude"],
            glottolog[l2]["latitude"],
            glottolog[l2]["longitude"],
        )

        # Set the value in the distance matrix
        matrix.set(l1, l2, int(d))

    # Write the distance `matrix`` to a file
    logging.info("[build_geodistance] Writing distance matrix to MATRIX.BZ2 file")
    matrix.save(ROOT_PATH / "src" / "samba_sampler" / "etc" / "haversine.matrix.bz2")


def build_gled_tree():
    """
    Build a global tree from the GLED trees.
    """

    logging.info("[build_gled_tree] Building global tree from GLED trees")

    # Obtain a list of all tree files
    tree_files = glob.glob("gled.trees/*.tree")

    # iterate over all tree files to find the maximum branch length
    max_branch_lengths = []
    for tree_file in tree_files:
        max_branch_length = 0
        t = Tree(str(tree_file))
        for node in t.traverse("postorder"):
            if node.get_distance(t) > max_branch_length:
                max_branch_length = node.get_distance(t)
        max_branch_lengths.append(max_branch_length)

    # Obtain the maximum branch length (that will be used a length from all roots
    # to the root of the global tree) and the mean max branch length (that will
    # be added, on top of `max_branch_length`, to isolates)
    max_branch_length = max(max_branch_lengths)
    mean_max_branch_length = sum(max_branch_lengths) / len(max_branch_lengths)

    # Read the glottolog.tsv file
    glottolog_dump = sorted(ROOT_PATH.glob("src/samba_sampler/etc/glottolog.*.tsv"))[-1]
    isolates = []
    with open(glottolog_dump, "r", encoding="utf-8") as f:
        tsv_reader = csv.reader(f, delimiter="\t")
        headers = next(tsv_reader)  # get the headers
        for row in tsv_reader:
            row_dict = dict(zip(headers, row))
            if row_dict["isolate"] == "True":
                isolates.append(row_dict["glottocode"])

    # create the global tree
    global_tree = Tree()

    # add each tree from disk to the global tree; note that this will label all
    # internal nodes with unique names (while not necessary for the code, it
    # makes the tree easier to read and debug)
    for tree_file in tree_files:
        t = Tree(tree_file)

        # derive a base name for labeling internal nodes
        base_name = os.path.splitext(os.path.basename(tree_file))[0]

        # counter for unique node names
        node_counter = 1

        # Change the leaf labels to just the Glottocode and label internal nodes
        for node in t.traverse():
            if node.is_leaf():
                node.name = node.name.split("_")[1]
            else:  # node is internal
                node.name = f"__{base_name}_{node_counter}__"
                node_counter += 1

        t.dist = max_branch_length
        global_tree.add_child(t)

    # add isolates to the global tree
    for isolate in isolates:
        t = Tree(name=isolate)
        t.dist = max_branch_length + mean_max_branch_length
        global_tree.add_child(t)

    # count the number of leaves and log info
    num_taxa = len(global_tree.get_leaves())
    logging.info(f"[build_gled_tree] Global tree has {num_taxa} taxa")

    # write the global tree to a file
    output = ROOT_PATH / "src" / "samba_sampler" / "etc" / "global_tree.gled.newick"
    global_tree.write(outfile=output, format=1)


def build_gled_matrix():
    logging.info("[build_gled_matrix] Building global matrix from the global GLED tree")
    with open(
        ROOT_PATH / "src" / "samba_sampler" / "etc" / "global_tree.gled.newick",
        encoding="utf-8",
    ) as f:
        tree = newick.load(f)[0]

    matrix = tree2matrix(tree)

    # Write the matrix to a file
    logging.info("[build_gled_matrix] Writing global matrix to file")
    matrix.save(ROOT_PATH / "src" / "samba_sampler" / "etc" / "gled.matrix.bz2")


def main():
    """
    Main function.
    """

    # Build the supplementary geographic data for language sampling
    build_geodistance()

    # Build the GLED tree and matrix
    build_gled_tree()
    build_gled_matrix()


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
