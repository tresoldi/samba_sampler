#!/usr/bin/env python3

"""
Build supplementary data for the scaffolding.
"""

# Import Python standard libraries
from pathlib import Path
import csv
import glob
import itertools
import math
import gzip

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

    # Compute the Haversine distance between all pairs of languages
    dist = {}
    for idx, (l1, l2) in enumerate(itertools.combinations(glottocodes, 2)):
        if idx % 100000 == 0:
            print(
                "Computing distance between {} and {} ({:.2f}%)".format(
                    l1, l2, 100 * idx / n_combinations
                )
            )
        d = haversine(
            glottolog[l1]["latitude"],
            glottolog[l1]["longitude"],
            glottolog[l2]["latitude"],
            glottolog[l2]["longitude"],
        )
        dist[l1, l2] = d

    # Build the distance matrix
    print("Building distance matrix")
    matrix = {}
    for l1 in glottocodes:
        for l2 in glottocodes:
            if l1 == l2:
                matrix[l2, l1] = 0
            elif (l1, l2) in dist:
                matrix[l2, l1] = int(dist[l1, l2])
            else:
                matrix[l2, l1] = None

    # Write the distance `matrix`` to a file
    print("Writing distance matrix to DST.GZ file")
    with gzip.open(
        ROOT_PATH / "src" / "samba_sampler" / "etc" / "haversine.dst.gz",
        "wt",
        newline="",
    ) as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["glottocode"] + glottocodes)
        for l1 in glottocodes:
            writer.writerow([l1] + [matrix[l1, l2] for l2 in glottocodes])


def main():
    """
    Main function.
    """

    # Build the supplementary geographic data for language sampling
    build_geodistance()


if __name__ == "__main__":
    main()
