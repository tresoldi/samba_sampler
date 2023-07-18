"""
test_common
==========

Tests for the `samba_sampler` package.
"""

# Import Python standard libraries
import csv
import pytest
from pathlib import Path

# Import the library being tested
import samba_sampler as samba


def test_trigger():
    assert 1 == 1


# def test_sampling():
#    sampler = samba.LanguageSampler()

# sampler = samba.GLED_Sampler()
# for idx, langset in enumerate(sampler.sample(4, 10)):
#    print(idx, langset)


def test_filter_data():
    def read_csv(filename):
        with open(filename, "r", encoding="utf-8") as f:
            return list(csv.DictReader(f, delimiter="\t"))

    source = (
        Path(__file__).parent.parent
        / "src"
        / "samba_sampler"
        / "etc"
        / "glottolog.v4.7-30-gafab84d5c4.tsv"
    )
    data = read_csv(source)

    filter_string = "family=Indo-European;ancestors=>germ1287,<nort3175"
    filtered_data = samba.filter_data(filter_string, data)

    # a. Test number of post-filtered entries
    assert len(filtered_data) == 345

    # b. Test presence of a single given entry
    assert any(d["glottocode"] == "swed1254" for d in filtered_data)  # Swedish

    # c. Test absence of a single given entry
    assert all(d["glottocode"] != "olde1238" for d in filtered_data)  # Old English
