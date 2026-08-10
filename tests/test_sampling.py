"""
Tests for the `sampling` module.
"""

import random

import pytest

from samba_sampler.common import DistanceMatrix
from samba_sampler.sampling import Sampler


def _matrix(keys):
    """Build a distance matrix whose distances grow with index difference."""
    matrix = DistanceMatrix(keys=keys)
    for i, ki in enumerate(matrix.keys):
        for j, kj in enumerate(matrix.keys):
            if j > i:
                matrix[ki, kj] = float(j - i)
    return matrix


KEYS = [f"t{i:02d}" for i in range(20)]


def test_requires_matrix_or_table():
    with pytest.raises(ValueError):
        Sampler()


def test_sampler_with_matrices_only():
    # Regression: constructing with only matrices used to raise TypeError
    sampler = Sampler(matrices=[_matrix(KEYS)])
    result = list(sampler.sample(5, seed=1))
    assert len(result) == 1
    assert len(result[0]) == 5
    assert len(set(result[0])) == 5  # no repeated taxa within a sample


def test_sampler_with_tables_only():
    table = {k: float(i) for i, k in enumerate(KEYS)}
    sampler = Sampler(tables=[table])
    result = list(sampler.sample(4, seed=1))
    assert len(result[0]) == 4


def test_sample_count_and_size():
    sampler = Sampler(matrices=[_matrix(KEYS)])
    result = list(sampler.sample(3, n=4, seed=7))
    assert len(result) == 4
    assert all(len(s) == 3 for s in result)


def test_progressive_algorithm():
    sampler = Sampler(matrices=[_matrix(KEYS)])
    result = list(sampler.sample(4, algorithm="progressive", seed=3))
    assert len(result[0]) == 4
    assert len(set(result[0])) == 4


def test_invalid_algorithm():
    sampler = Sampler(matrices=[_matrix(KEYS)])
    with pytest.raises(ValueError):
        list(sampler.sample(3, algorithm="does-not-exist"))


def test_invalid_summary():
    sampler = Sampler(matrices=[_matrix(KEYS)])
    with pytest.raises(ValueError):
        list(sampler.sample(3, summary="median"))


def test_reproducibility_with_seed():
    sampler = Sampler(matrices=[_matrix(KEYS)])
    first = list(sampler.sample(4, n=3, seed=123))
    second = list(sampler.sample(4, n=3, seed=123))
    assert first == second


def test_sampling_does_not_disturb_global_rng():
    sampler = Sampler(matrices=[_matrix(KEYS)])

    random.seed(1234)
    expected = [random.random() for _ in range(3)]

    random.seed(1234)
    first = random.random()
    # Drawing a sample must not consume or reset the global RNG state.
    list(sampler.sample(4, n=2, seed=999))
    rest = [random.random() for _ in range(2)]

    assert [first] + rest == expected


def test_taxa_subset_restricts_keys():
    subset = KEYS[:6]
    sampler = Sampler(matrices=[_matrix(KEYS)], taxa_subset=subset)
    result = list(sampler.sample(4, n=3, seed=5))
    sampled = {taxon for s in result for taxon in s}
    assert sampled.issubset(set(subset))
