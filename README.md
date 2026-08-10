# SAMBA - Sampling Algorithms with Matrix-Based Weight Allocation

SAMBA is a Python package providing sampling methods via matrix-based distance measures to mitigate autocorrelation.

## Installation and Usage

The package can be installed via pip:

```bash
$ pip install samba_sampler
```

Detailed information on the usage of the package can be found in the [documentation](https://samba.readthedocs.io/en/latest/).
For a quick start, the following example shows how to use the package:

```python
import samba_sampler as samba
sampler = samba.LanguageSampler() # Default parameters
print(list(sampler.sample(5)))
```

To also spread a sample across genealogical **genera** (a mid-level unit
typologists use to control for genealogical bias), use `GenusSampler`, which
adds a bundled glottocode-to-genus map on top of the phylogenetic and geographic
matrices:

```python
import samba_sampler as samba
sampler = samba.GenusSampler(group_weight=1.0)
print(list(sampler.sample(5)))

# The mapping is also available directly:
genus_of = samba.load_genera()
```

## Changelog

Version 0.4.0
  - Added genus-aware sampling: a `GenusSampler` and an optional group-diversity
    term that spreads samples across genealogical genera to reduce genealogical
    autocorrelation
  - Bundled a glottocode-to-genus map (from the `genera` project: a calibrated
    Glottolog-depth cut, validated with WALS and supported by GLED and ASJP) and
    exposed it via `load_genera()`

Version 0.3.3
  - Fixed `LanguageSampler` so the quick-start example runs without external files
  - Fixed a bug in `tree2matrix` where terminal branch lengths were omitted from patristic distances
  - Fixed a crash when constructing a `Sampler` with only matrices or only tables
  - Isolated the random number generator so sampling no longer mutates global `random` state
  - Bounded internal caches, made table columns configurable, and corrected `DistanceMatrix.rescale`
  - Corrected the package's public API (`__all__`), raised the minimum Python to 3.9, and added a test suite and CI

Version 0.3.2 (2023-07-18)
  - Added taxa filtering, applying to the Glottolog data

Version 0.3.1 (2023-07-17)
  - Added first version of the walking distance matrix derived from Guzman Naranjo & Jäger (2023).

Version 0.3 (2023-07-13)
  - Initial release, following on the `arcaverborum` project.

## Community Guidelines

While the author can be contacted directly for support, it is recommended that third parties use GitHub standard features, such as issues and pull requests, to contribute, report problems, or seek support.

Contributing guidelines, including a code of conduct, can be found in the CONTRIBUTING.md file.

## Author, Citation, and Acknowledgements

The library is developed by Tiago Tresoldi (tiago@tresoldi.org).

The library is developed in the context of the Cultural Evolution of Texts project, with funding from the Riksbankens Jubileumsfond (grant agreement ID: MXM19-1087:1).

If you use `samba_sampler`, please cite it as:

> Tresoldi, Tiago (2023). SAMBA (Sampling Algorithms with Matrix-Based Weight Allocation): a Python package providing sampling methods via matrix-based distance measures to mitigate autocorrelation. Version 0.3. Uppsala: Uppsala University.

In BibTeX:

```bibtex
@misc{Tresoldi2023samba,
  author = {Tresoldi, Tiago},
  title = {SAMBA (Sampling Algorithms with Matrix-Based Weight Allocation): a Python package providing sampling methods via matrix-based distance measures to mitigate autocorrelation. Version 0.3.},
  howpublished = {\url{https://github.com/tresoldi/samba_sampler}},
  address = {Uppsala},
  published = {Upssala University},
  year = {2023}
}
```
