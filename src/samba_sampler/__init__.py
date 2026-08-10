# __init__.py for the `samba_sampler` package

# Metadata. This is the single source of truth for the version; `setup.cfg`
# reads it via `attr: samba_sampler.__version__`.
__version__ = "0.4.0"
__author__ = "Tiago Tresoldi"
__email__ = "tiago@tresoldi.org"

# Local imports
from .common import (
    tree2matrix,
    DistanceMatrix,
    build_table_from_file,
    filter_glottolog,
    filter_data,
    load_genera,
)
from .sampling import Sampler, LanguageSampler, GenericSampler, GenusSampler
from .newick import loads, dumps, load, dump

# Expose the public API
__all__ = [
    "tree2matrix",
    "DistanceMatrix",
    "build_table_from_file",
    "filter_glottolog",
    "filter_data",
    "load_genera",
    "Sampler",
    "LanguageSampler",
    "GenericSampler",
    "GenusSampler",
    "loads",
    "dumps",
    "load",
    "dump",
]
