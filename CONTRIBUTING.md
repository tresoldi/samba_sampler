# Contributing to SAMBA

Thank you for your interest in contributing to `samba_sampler`. Contributions of
all kinds are welcome: bug reports, documentation improvements, new features, and
test coverage.

## Reporting issues

Please use the [GitHub issue tracker](https://github.com/tresoldi/samba_sampler/issues)
to report bugs or request features. When reporting a bug, include:

- The version of `samba_sampler` and of Python you are using.
- A minimal, self-contained example that reproduces the problem.
- The full traceback, if an exception is raised.

## Development setup

Clone the repository and install the package in editable mode together with the
development and test extras:

```bash
git clone https://github.com/tresoldi/samba_sampler.git
cd samba_sampler
python -m pip install -e ".[dev,test]"
```

## Running the tests

The test suite is run with `pytest` and is configured to import the package
directly from `src/`, so no installation step is required:

```bash
pytest
```

Please make sure the full test suite passes and add tests covering any new
behaviour before opening a pull request.

## Code style

This project uses [black](https://black.readthedocs.io/) for formatting and
[flake8](https://flake8.pycqa.org/) for linting:

```bash
black src tests
flake8 src tests
```

## Pull requests

1. Fork the repository and create a topic branch for your change.
2. Keep changes focused; unrelated changes should go in separate pull requests.
3. Update documentation and the changelog in `README.md` when relevant.
4. Ensure the tests pass and the code is formatted before requesting review.

## Code of conduct

Everyone participating in this project is expected to be respectful and
considerate. Harassment or abusive behaviour of any kind will not be tolerated.
Please treat other contributors as you would like to be treated. Report any
unacceptable behaviour to the maintainer at tiago@tresoldi.org.
