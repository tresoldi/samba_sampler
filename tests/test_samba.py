"""
test_common
==========

Tests for the `samba_sampler` package.
"""


# Import the library being tested
import samba_sampler as samba


def test_trigger():
    assert 1 == 1


def test_sampling():
    sampler = samba.GLED_Sampler()