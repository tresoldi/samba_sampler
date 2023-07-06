# Scaffolding

This directory hosts a collection of scripts responsible for downloading and processing external data, including CLDF datasets, and generating our own data, such as geodistances between languages.

It's important to note that these scripts are not necessary for end users to run, as the data is readily available in the repository. They are primarily beneficial for developers who wish to update the data or gain insights into its generation. Notably, these scripts eliminate the need to undergo the complicated process of installing the CLDF framework and its associated dependencies, a task currently incompatible with Windows systems.

Maintaining a local copy of the data, complete with metadata and especially version records, ensures our analysis pipeline's independence from external sources. This approach also simplifies automatic continuous deployment. Aiming towards this goal, we strive to minimize code and data dependencies required to compile our data, as the tracked versions will be stored in this location.

## Instructions

Create a virtual environment that supports CLDF catalogs and retrieve the necessary data. You can accomplish this with the commands provided below.
Please note, the `cldfbench catconfig -q` command creates local clones of the repositories on your hard drive, defaulting to `~/.config/cldf`.
This process may take several minutes and can require many gigabytes of storage space.

```bash
python -m venv env
source env/bin/activate
pip install -U pip setuptools wheel
pip install -r requirements.txt
cldfbench catconfig -q
```

Make local copies of the catalogues as stand-alone tabular files:

```bash
python process_catalogs.py
```

Build and distribute in the repository supplementary data (such as
geodistances):

```bash
python build_supplements.py
```

```bash
python build_trees.py
```