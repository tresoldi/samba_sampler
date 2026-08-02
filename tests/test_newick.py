"""
Tests for the vendored `newick` module.
"""

import samba_sampler as samba

SIMPLE = "(A:1.0,B:2.0,(C:3.0,D:4.0)E:5.0)F;"


def test_loads_leaf_names():
    tree = samba.loads(SIMPLE)[0]
    assert sorted(tree.get_leaf_names()) == ["A", "B", "C", "D"]


def test_roundtrip_is_stable():
    tree = samba.loads(SIMPLE)[0]
    once = samba.dumps(tree)
    twice = samba.dumps(samba.loads(once)[0])
    assert once == twice


def test_leaf_structure_and_lengths():
    tree = samba.loads(SIMPLE)[0]
    lengths = {n.name: n.length for n in tree.get_leaves()}
    assert lengths["A"] == 1.0
    assert lengths["D"] == 4.0


def test_load_and_dump_file(tmp_path):
    path = tmp_path / "tree.newick"
    samba.dump(samba.loads(SIMPLE)[0], path.open("w", encoding="utf-8"))
    with path.open(encoding="utf-8") as fp:
        tree = samba.load(fp)[0]
    assert sorted(tree.get_leaf_names()) == ["A", "B", "C", "D"]
