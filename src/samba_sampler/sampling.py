"""
Functions for obtaining language samples.
"""

# Import from standard library
from typing import Optional, Union, List
import collections
import itertools
import random

# Import from other modules
from .common import ETC_PATH, DistanceMatrix


class Sampler:
    def __init__(self, matrices: List[DistanceMatrix], weights=None):
        """
        Initialize the sampler.
        :param matrices:
        :param weights:
        """

        # Store the matrices internally
        self._matrices = matrices

        # Store the set of overlapping keys among the matrices as a sorted list
        self._keys = sorted(set.intersection(*[set(m.keys) for m in matrices]))

        # If `weights` is not provided, default to 1.0 for every item in `matrices`
        if not weights:
            self._weights = [1.0] * len(matrices)
        else:
            self._weights = weights

    # TODO: have yield instead of return
    def sample(
        self,
        k: int,
        t: int = 1,
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

        # Initialize the set of sampled sets
        sampled = set()
        taxa_counter = collections.Counter()
        for i in range(t):
            # Run different attempts to obtain one good sample
            candidates = []
            for j in range(tries):
                # Obtain a random sample of k taxa
                sampled_taxa = tuple(random.sample(self._keys, k))

                # Compute the distances between each pair of taxa
                # pdists = []
                # gdists = []
                dists = [[] for i in range(len(self._matrices))]
                for lang1, lang2 in itertools.combinations(sampled_taxa, 2):
                    # pdists.append(self._matrices[0][lang1,lang2])
                    # gdists.append(self._matrices[1][lang1,lang2])
                    for m in range(len(self._matrices)):
                        dists[m].append(self._matrices[m][lang1, lang2])

                # Compute the frequency of each taxon in the sample
                # TODO: use a Witten-Bell smoothing or something similar from lpngrams,
                #       also accounting for pairs, etc.
                fdists = [taxa_counter[taxon] for taxon in sampled_taxa]

                # Initialize a list to store the calculated distances
                calculated_dists = []

                for dists_list in dists:
                    # Calculate the distance depending on the method
                    if method == "mean":
                        calculated_dist = sum(dists_list) / k
                    else:  # sum
                        calculated_dist = sum(dists_list)

                    calculated_dists.append(calculated_dist)

                # Compute the fdist
                if not taxa_counter.most_common(1):
                    fdist = 1.0
                else:
                    if method == "mean":
                        fdist = 1 - (
                            (sum(fdists) / k) / (taxa_counter.most_common(1)[0][1])
                        )
                    else:  # sum
                        fdist = 1 - (
                            sum(fdists)
                            / sum([e[1] for e in taxa_counter.most_common(k)])
                        )

                # Now, calculated_dists is a list containing all the calculated distances.
                # We just need to use these in the calculation that is appended to candidates.
                score = (
                    sum(
                        dist * weight
                        for dist, weight in zip(calculated_dists, self._weights)
                    )
                    + fdist * freq_weight
                )
                candidates.append((sampled_taxa, score))

            # Perform a weighted sampling of the candidates, using the second
            # element in each tuple
            random_set = random.choices(*zip(*candidates), k=1)[0]
            sampled.add(random_set)

            # Update the tracking of the sampled languages across different sets
            taxa_counter.update(random_set)

        # Return the sampled language sets
        return sampled


class LanguageSampler(Sampler):
    def __init__(self):
        """
        Initialize the sampler.
        """

        phylomatrix = DistanceMatrix(filename=ETC_PATH / "gled.matrix.bz2")
        print(f"Read {len(phylomatrix.data)} entries from phylogenetic matrix")
        print(phylomatrix.keys[:20])

        geomatrix = DistanceMatrix(filename=ETC_PATH / "haversine.matrix.bz2")
        print(f"Read {len(geomatrix.data)} entries from geographic matrix")
        print(geomatrix.keys[:20])

        # TODO: decide on factors
        phylomatrix.rescale()
        # geomatrix.rescale()

        # Initialize the parent class
        super().__init__([phylomatrix, geomatrix])
