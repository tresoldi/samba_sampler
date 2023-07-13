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

## Changelog

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
