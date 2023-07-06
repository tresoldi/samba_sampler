# __init__.py for the `samba_sampler` package
# Path: src\samba_sampler\__init__.py

# Metadata
__version__ = "2.0"  # Remember to update this in setup.py
__author__ = "Tiago Tresoldi"
__email__ = "tiago@tresoldi.org"

# Local imports
from .common import read_splitstree_matrix, read_default_matrix, read_triangle_matrix, tree2matrix, DistanceMatrix
from .sampling import GLED_Sampler, MyLanguageSampler
from .newick import loads, dumps, load, dump

# Expose the functions
all = [
    "loads", "load", "dumps", "dump",
"newick2matrix", "DistanceMatrix",
    "read_matrix", "read_default_matrix", "read_triangle_matrix", "GLED_Sampler", "MyLanguageSampler", "Language_Sampler"]
