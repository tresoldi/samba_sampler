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
print(sampler.sample(5))
```

## Showcases

(...)

## Changelog

Version 0.1.0 (2021-03-01)
  - Initial release

## Community Guidelines

While the author can be contacted directly for support, it is recommended that third parties use GitHub standard features, such as issues and pull requests, to contribute, report problems, or seek support.

Contributing guidelines, including a code of conduct, can be found in the CONTRIBUTING.md file.

## Author, Citation, and Acknowledgements

The library is developed by Tiago Tresoldi (tiago@tresoldi.org).

The library is developed in the context of the Cultural Evolution of Texts project, with funding from the Riksbankens Jubileumsfond (grant agreement ID: MXM19-1087:1).

If you use `dafsa`, please cite it as:

> Tresoldi, Tiago (2020). DAFSA, a a library for computing Deterministic Acyclic Finite State Automata. Version 1.0. Jena. DOI: 10.5281/zenodo.3668870

In BibTeX:

```bibtex
@misc{Tresoldi2020dafsa,
  author = {Tresoldi, Tiago},
  title = {DAFSA, a a library for computing Deterministic Acyclic Finite State Automata. Version 1.0.},
  howpublished = {\url{https://github.com/tresoldi/dafsa}},
  address = {Jena},
  doi = {10.5281/zenodo.3668870}.
  year = {2020},
}
```

## TODO

- different distance for family in the gled tree (not over 1000)