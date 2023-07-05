#!/usr/bin/env python

"""
Obtain the CLDF catalogs (CLTS, Concepticon, Glottolog) as local copies.
"""

# Import Python standard libraries
from pathlib import Path
from typing import *
import csv
import logging

# Import MPI-SHH libraries
from pyglottolog import Glottolog
from clldutils.path import git_describe

# Define the paths used by the script
BASE_PATH = Path(__file__).parent
ROOT_PATH = BASE_PATH.parent


def get_glottolog_path(glottolog_path: Optional[Union[Path, str]] = None) -> str:
    if not glottolog_path:
        glottolog_path = Path.home() / ".config" / "cldf" / "glottolog"
        glottolog_path = glottolog_path.as_posix()
    return glottolog_path


def get_languoids(glottolog_path: str):
    glottolog = Glottolog(glottolog_path)
    languoids = [languoid.glottocode for languoid in glottolog.languoids()]
    return glottolog, languoids


def process_languoid(languoid, corrections: dict):
    names = process_names(languoid)
    iso, closest_iso = process_iso_codes(languoid)
    latitude, longitude = process_lat_long(languoid)
    timespan = languoid.timespan if languoid.timespan else ""

    entry = {
        "ancestors": ";".join([anc.glottocode for anc in languoid.ancestors]),
        "category": languoid.category,
        "children": ";".join([child.glottocode for child in languoid.children]),
        "closest_iso": closest_iso,
        "countries": ";".join([country.name for country in languoid.countries]),
        "glottocode": languoid.glottocode,
        "iso": iso,
        "isolate": str(languoid.isolate),
        "latitude": latitude,
        "longitude": longitude,
        "level": languoid.level.id,
        "macroareas": ";".join([m.name for m in languoid.macroareas]),
        "name": languoid.name,
        "names": names,
        "timespan": timespan,
    }

    if languoid.glottocode in corrections:
        entry.update(
            corrections[languoid.glottocode]
        )  # overriding the values with corrections

    return entry


def process_names(languoid):
    if languoid.names:
        names = ";".join(languoid.names.get("multitree", []))
        names = names.replace('"', "")
    else:
        names = ""
    return names


def process_iso_codes(languoid):
    iso = languoid.iso if languoid.iso else ""
    closest_iso = languoid.closest_iso() if languoid.closest_iso() else ""
    return iso, closest_iso


def process_lat_long(languoid):
    latitude = (
        languoid.latitude if languoid.latitude and languoid.latitude != "None" else ""
    )
    longitude = (
        languoid.longitude
        if languoid.longitude and languoid.longitude != "None"
        else ""
    )
    return latitude, longitude


def process_family(data, glottocode2name):
    for row in data:
        if row["ancestors"]:
            ancestors = row["ancestors"].split(";")
            row["family"] = glottocode2name[ancestors[0]]
        else:
            row["family"] = ""


def write_results_to_disk(data, glottolog):
    version = git_describe(glottolog.repos)
    fields = [
        "glottocode",
        "name",
        "iso",
        "closest_iso",
        "family",
        "category",
        "level",
        "isolate",
        "macroareas",
        "countries",
        "latitude",
        "longitude",
        "ancestors",
        "children",
        "names",
        "timespan",
    ]
    with open(
        ROOT_PATH / "src" / "samba_sampler" / "etc" / f"glottolog.{version}.tsv",
        "w",
        newline="",
    ) as f:
        writer = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        writer.writeheader()
        writer.writerows(data)


def load_corrections():
    corrections = {}
    corrections_path = BASE_PATH / "etc" / "glottolog.corrections.csv"
    if corrections_path.exists():
        logging.info("Loading Glottolog corrections...")
        with open(corrections_path, "r") as file:
            reader = csv.DictReader(file)
            for row in reader:
                corrections[row["glottocode"]] = {
                    k: v for k, v in row.items() if v
                }  # exclude empty values

    return corrections


import concurrent.futures


def worker(glottocode, glottolog_path, corrections):
    glottolog = Glottolog(glottolog_path)
    languoid = glottolog.languoid(glottocode)
    return process_languoid(languoid, corrections)


def get_glottolog(glottolog_path: Optional[Union[Path, str]] = None):
    """
    Make a local tabular copy of the Glottolog data.
    Note that we intentionally exclude some fields, particularly those
    related to sources and classification comments: the goals is *not*
    to make a full copy or clone of the Glottolog data, but rather to
    make a tabular copy of the data that is relevant for the purposes
    of our projects.
    """
    glottolog_path = get_glottolog_path(glottolog_path)
    glottolog, languoids = get_languoids(glottolog_path)
    num_languoids = len(languoids)

    corrections = load_corrections()
    glottocode2name = {}
    data = []

    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = {
            executor.submit(worker, glottocode, glottolog_path, corrections): glottocode
            for glottocode in languoids
        }

        for i, future in enumerate(concurrent.futures.as_completed(futures)):
            data.append(future.result())
            print(f"Processed {i + 1} out of {num_languoids} entries")

    # Build glottocode2name
    for row in data:
        glottocode2name[row["glottocode"]] = row["name"]

    process_family(data, glottocode2name)
    data = sorted(data, key=lambda x: x["glottocode"])
    write_results_to_disk(data, glottolog)


def main():
    """
    Script entry point.
    """

    # Obtain the Glottolog data dump, if necessary
    # TODO: have a better approach that checks for the version
    glottolog_dumps = list(ROOT_PATH.glob("etc/glottolog.*.tsv"))
    if not glottolog_dumps:
        logging.info("No Glottolog dump found, generating it...")
        get_glottolog()
    else:
        logging.info("Glottolog dump already exists, skipping...")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
