"""
Tests for the `common` module: DistanceMatrix, table building, and tree2matrix.
"""

import math

import pytest

import samba_sampler as samba
from samba_sampler.common import DistanceMatrix


def _build_matrix():
    """Build a small 4x4 distance matrix with known values."""
    keys = ["a", "b", "c", "d"]
    matrix = DistanceMatrix(keys=keys)
    values = {
        ("a", "b"): 1.0,
        ("a", "c"): 2.0,
        ("a", "d"): 3.0,
        ("b", "c"): 4.0,
        ("b", "d"): 5.0,
        ("c", "d"): 6.0,
    }
    for (i, j), v in values.items():
        matrix[i, j] = v
    return matrix, values


def test_distance_matrix_symmetry_and_diagonal():
    matrix, values = _build_matrix()

    # Diagonal is always zero
    for key in matrix.keys:
        assert matrix[key, key] == 0

    # Values are symmetric
    for (i, j), v in values.items():
        assert matrix[i, j] == pytest.approx(v)
        assert matrix[j, i] == pytest.approx(v)


def test_distance_matrix_requires_exactly_one_source():
    with pytest.raises(ValueError):
        DistanceMatrix()
    with pytest.raises(ValueError):
        DistanceMatrix(keys=["a"], filename="foo.bz2")


def test_distance_matrix_save_load_roundtrip(tmp_path):
    matrix, values = _build_matrix()
    path = tmp_path / "matrix.bz2"
    matrix.save(path)

    loaded = DistanceMatrix(filename=path)
    assert loaded.keys == matrix.keys
    for (i, j), v in values.items():
        assert loaded[i, j] == pytest.approx(v)


def test_rescale_maps_to_range():
    matrix, _ = _build_matrix()
    matrix.rescale((0.0, 1.0))
    assert min(matrix.data) == pytest.approx(0.0)
    assert max(matrix.data) == pytest.approx(1.0)


def test_rescale_with_custom_range_and_factor():
    matrix, _ = _build_matrix()
    matrix.rescale((1.0, 2.0), factor=2.0)
    # Range [1, 2] scaled by factor 2 -> [2, 4]
    assert min(matrix.data) == pytest.approx(2.0)
    assert max(matrix.data) == pytest.approx(4.0)


def test_rescale_constant_matrix_does_not_divide_by_zero():
    matrix = DistanceMatrix(keys=["a", "b", "c"])
    for i, j in [("a", "b"), ("a", "c"), ("b", "c")]:
        matrix[i, j] = 5.0
    matrix.rescale((0.0, 1.0))
    assert all(v == pytest.approx(0.0) for v in matrix.data)


def test_build_table_from_file(tmp_path):
    path = tmp_path / "table.csv"
    path.write_text(
        "Glottocode,weight\nabcd1234,1.0\nabcd1234,3.0\nefgh5678,2.0\n",
        encoding="utf-8",
    )
    table = samba.build_table_from_file(path, "Glottocode", "weight")
    # Duplicate keys are averaged by default
    assert table["abcd1234"] == pytest.approx(2.0)
    assert table["efgh5678"] == pytest.approx(2.0)


def test_build_table_from_file_rejects_bad_multiple(tmp_path):
    path = tmp_path / "table.csv"
    path.write_text("Glottocode,weight\nabcd1234,1.0\n", encoding="utf-8")
    with pytest.raises(ValueError):
        samba.build_table_from_file(path, "Glottocode", "weight", multiple="median")


def test_tree2matrix_known_distances():
    tree = samba.loads(
        "(((culi1244:0.216844,deni1241:0.216844):0.0787837,"
        "jama1261:0.29563):1.07742,paum1247:1.37305);"
    )[0]
    matrix = samba.tree2matrix(tree)

    # Sister taxa: sum of the two terminal branch lengths to their shared parent
    assert matrix["culi1244", "deni1241"] == pytest.approx(0.433688, abs=1e-5)
    # Patristic distance across the tree: culi1244 -> parent -> jama1261
    assert matrix["culi1244", "jama1261"] == pytest.approx(0.591258, abs=1e-5)
    # Distances are symmetric
    assert matrix["culi1244", "jama1261"] == pytest.approx(
        matrix["jama1261", "culi1244"]
    )
    # A leaf is at distance zero from itself
    assert matrix["paum1247", "paum1247"] == 0
    assert not math.isnan(matrix["culi1244", "paum1247"])


def test_load_genera_bundled():
    genus_of = samba.load_genera()
    assert len(genus_of) > 5000
    # Swedish is an in-scope Germanic language present in the bundled map.
    assert "swed1254" in genus_of

    genus_of2, names = samba.load_genera(with_names=True)
    assert genus_of2 == genus_of
    assert names.get(genus_of["swed1254"])  # its genus carries a name
