"""
Functions for obtaining language samples.
"""

# Import from standard library
from typing import *
import collections
import itertools
import logging
import math
import random

# Import from other modules
from .common import read_splitstree_matrix, read_triangle_matrix, ETC_PATH


class Sampler:
    def __init__(self, matrices, weights=None, seed=None):
        """
        Initialize the sampler.
        :param matrices:
        :param weights:
        """

        # Store the matrices internally
        self._matrices = matrices

        # If `weights` is not provided, default to 1.0 for every item in `matrices`
        if not weights:
            self._weights = [1.0] * len(matrices)
        else:
            self._weights = weights

        # Seed the random number generator
        random.seed(seed)


class MyLanguageSampler(Sampler):
    def __init__(self):
        """
        Initialize the sampler.
        """

        phylomatrix = read_splitstree_matrix(ETC_PATH / "gled_global.dst")
        print(f"Read {len(phylomatrix)} entries from phylogenetic matrix")
        print(phylomatrix.keys())
        print(phylomatrix["Zaparo_zapa1253"].keys())

        print(len(phylomatrix.keys()))
        print(len(phylomatrix["Zaparo_zapa1253"].keys()))

        # geomatrix = read_triangle_matrix(ETC_PATH / "haversine.dst.gz")
        # print(f"Read {len(geomatrix)} entries from geographic matrix")

        # Initialize the parent class
        # super().__init__([phylomatrix, geomatrix])


class Language_Sampler:
    """
    A class for obtaining language samples.

    The usage of a class allows, particularly for complex models,
    to pre-compute some values that can be reused across different
    samples and to save loading time.
    """

    def __init__(self):
        """
        Initialize the sampler.
        """

        # Initialize the matrices for sampling
        self._pmatrix = None  # Phylogenetic matrix
        self._gmatrix = None  # Geographic matrix

        # List of the languages in the matrix
        self._lects = None

    # TODO: have yield instead of return
    def sample(
        self,
        k: int,
        t: int = 1,
        phylo_weight: float = 1.0,
        geo_weight: float = 1.0,
        freq_weight: float = 1.0,
        method: str = "mean",
        tries: int = 100,
        seed: Optional[Union[int, str]] = None,
    ):
        """
        Obtain a language sample from a distance matrix.

        Parameters
        ----------
        k
            The number of languages to sample in each set.
        t
            The number of sets to sample.
        phylo_weight
            The weight to give to the phylogenetic distance. Set to 0 to
            ignore the phylogenetic distance.
        geo_weight
            The weight to give to the geographic distance. Set to 0 to
            ignore the geographic distance.
        freq_weight
            The (penalty) weight to give to the frequency of the language in
            previous samples. Set to 0 to ignore the frequency.
        method
            The method to use to compute the weight of a sample. Can be either
            "mean" or "sum".
        tries
            The number of attempts to make to obtain a good sample.
        seed
            The seed to use for the random number generator. If None, the seed
            will be set to the current time.

        Returns
        -------
        sample
            A tuple of the sampled languages.
        """

        # Make sure we are using a valid method
        if method not in ["mean", "sum"]:
            raise ValueError("Invalid method in `sample()`: {}".format(method))

        # Set the seed for the random number generator
        random.seed(seed)

        # If the list of lects has not been extracted and cached (this is the
        # first call to the function), extract it
        if self._lects is None:
            self._lects = sorted(self._pmatrix.keys())

        # Initialize the set of sampled language sets
        sampled = set()
        lect_counter = collections.Counter()
        for i in range(t):
            # Run different attempts to obtain one good sample
            candidates = []
            for j in range(tries):
                # Obtain a random sample of k languages
                sampled_lects = tuple(random.sample(self._lects, k))

                # Compute the phylogenetic and geographic distances between each pair of languages
                pdists = []
                gdists = []
                for lang1, lang2 in itertools.combinations(sampled_lects, 2):
                    pdists.append(self._pmatrix[lang1][lang2])
                    gdists.append(self._gmatrix[lang1][lang2])

                # Compute the frequency of each language in the sample
                # TODO: use a Witten-Bell smoothing or something similar from lpngrams,
                #       also accounting for pairs, etc.
                fdists = [lect_counter[lang] for lang in sampled_lects]

                # Store the current candidate along with its combined weights
                if method == "mean":
                    pdist = sum(pdists) / k
                    gdist = sum(gdists) / k
                    if not lect_counter.most_common(1):
                        fdist = 1.0
                    else:
                        fdist = 1 - (
                            (sum(fdists) / k) / (lect_counter.most_common(1)[0][1])
                        )
                else:  # sum
                    pdist = sum(pdists)
                    gdist = sum(gdists)
                    if not lect_counter.most_common(1):
                        fdist = 1.0
                    else:
                        fdist = 1 - (
                            sum(fdists)
                            / sum([e[1] for e in lect_counter.most_common(k)])
                        )

                candidates.append(
                    (
                        sampled_lects,
                        pdist * phylo_weight + gdist * geo_weight + fdist * freq_weight,
                    )
                )

            # Perform a weighted sampling of the candidates, using the second
            # element in each tuple
            random_set = random.choices(*zip(*candidates), k=1)[0]
            sampled.add(random_set)

            # Update the tracking of the sampled languages across different sets
            lect_counter.update(random_set)

        # Return the sampled language sets
        return sampled

    # TODO: move to common? use `range`
    def _rescale_matrix(
        self, matrix: Dict[str, Dict[str, float]], scale_range=(0, 1), factor=1.0
    ):
        """
        Rescale a matrix to a given range.

        Parameters
        ----------
        matrix
            The matrix to rescale.
        scale_range
            The range to rescale to.
        factor
            An exponential factor to apply to the values in the matrix.

        Returns
        -------
        rescaled
            The rescaled matrix.
        """

        # Obtain the minimum and maximum values in the matrix
        minval, maxval = None, -1
        for inner_dict in matrix.values():
            for value in inner_dict.values():
                if minval is None or value < minval:
                    minval = value
                if value > maxval:
                    maxval = value

        # Rescale the matrix
        rescaled = {
            lang1: {
                lang2: math.pow(
                    (matrix[lang1][lang2] - minval) / (maxval - minval), factor
                )
                for lang2 in matrix[lang1].keys()
            }
            for lang1 in matrix.keys()
        }

        # Return the rescaled matrix
        return rescaled


class GLED_Sampler(Language_Sampler):
    """
    The default sampler using data from the GLED project.
    """

    def __init__(self):
        """
        Initialize the sampler.
        """

        # Initialize the parent class
        super().__init__()

        # Load the phylogenetic data
        logging.warning("Loading the phylogenetic matrix from GLED...")
        self._pmatrix = read_splitstree_matrix(ETC_PATH / "gled_global.dst")

        # Load the geographic data
        # NOTE: for now, we are post-processing the haversine data to
        #       select the languages in the GLED matrix; this will be
        #       improved in future editions also by incorporating Glottolog
        #       database dump
        logging.warning("Loading the geographic matrix from GLED...")
        haversine = read_triangle_matrix(ETC_PATH / "haversine.dst.gz")

        # Obtain a dictionary mapping the glottocodes in `self._pmatrix`
        # to their keys (which are in the format "{language}_{glottocode}"),
        # filter the `haversine` matrix accordingly, and rename its keys
        # using such a mapping to match the keys in the phylogenetic matrix
        mapping = {lect.split("_")[1]: lect for lect in self._pmatrix.keys()}
        self._gmatrix = {
            mapping[glottocode_a]: {
                mapping[glottocode_b]: dist
                for glottocode_b, dist in haversine[glottocode_a].items()
                if glottocode_b in mapping
            }
            for glottocode_a in mapping.keys()
        }

        # Rescale all values in the phylogenetic and geographic matrices in
        # the range [0, 1]
        logging.warning("Rescaling the phylogenetic matrix...")
        self._pmatrix = self._rescale_matrix(self._pmatrix)
        logging.warning("Rescaling the geographic matrix...")
        self._gmatrix = self._rescale_matrix(self._gmatrix, factor=0.5)
