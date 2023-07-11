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

    def sample(
        self,
        k: int,
        t: int = 1,
        freq_weight: float = 1.0,
        method: str = "mean",
        tries: int = 100,
        seed: Optional[Union[int, str]] = None,
    ):
        # Validate method
        if method not in ["mean", "sum"]:
            raise ValueError("Invalid method in `sample()`: {}".format(method))

        # Set the seed
        random.seed(seed)

        # Initialize the set of sampled sets
        taxa_counter = collections.Counter()

        for _ in range(t):
            candidates = []
            for _ in range(tries):
                # Obtain a random sample of k taxa
                sampled_taxa = tuple(random.sample(self._keys, k))

                # Compute the distances between each pair of taxa
                matrix_length = len(self._matrices)
                dists = [
                    [
                        self._matrices[m][lang1, lang2]
                        for lang1, lang2 in itertools.combinations(sampled_taxa, 2)
                    ]
                    for m in range(matrix_length)
                ]

                # Compute the frequency of each taxon in the sample
                fdists = [taxa_counter[taxon] for taxon in sampled_taxa]

                # Calculate the distances
                calculated_dists = [
                    sum(dist) / k if method == "mean" else sum(dist) for dist in dists
                ]

                # Compute fdist
                most_common = taxa_counter.most_common(1)
                if most_common:
                    most_common_count = most_common[0][1]
                    fdist = 1 - (
                        sum(fdists)
                        / (
                            k * most_common_count
                            if method == "mean"
                            else sum([e[1] for e in taxa_counter.most_common(k)])
                        )
                    )
                else:
                    fdist = 1.0

                # Calculate the final score
                score = (
                    sum(
                        dist * weight
                        for dist, weight in zip(calculated_dists, self._weights)
                    )
                    + fdist * freq_weight
                )
                candidates.append((sampled_taxa, score))

            # Perform a weighted sampling of the candidates
            candidates, weights = zip(*candidates)
            selected = random.choices(candidates, weights=weights, k=1)[0]

            # Update the tracking of the sampled languages across different sets
            taxa_counter.update(selected)

            # Yield the results
            yield selected


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
