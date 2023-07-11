"""
Functions for obtaining language samples.
"""

# Import from standard library
from typing import Optional, Union, List, Tuple, Sequence
import collections
import itertools
import random
import functools

# Import from other modules
from .common import ETC_PATH, DistanceMatrix, build_dict_from_file


class Sampler:
    def __init__(
        self,
        matrices: List[DistanceMatrix],
        tables: List[dict],
        mweights=None,
        tweights=None,
    ):
        """
        Initialize the sampler.
        :param matrices:
        :param mweights:
        """

        # Store the matrices internally
        self._matrices = matrices
        self._num_matrices = len(matrices)

        # Store the tables internally
        self._tables = tables
        self._num_tables = len(tables)

        # Raise an error if no matrix and no table was provided
        if self._num_matrices == 0 and self._num_tables == 0:
            raise ValueError("No matrix and no table provided.")

        # Obtain the set of keys for all matrices and for all tables, if they
        # were provided, and then compute the intersection of the sets,
        # which will be the set of keys to be used for sampling. The intersection
        # is stored as a sorted list in `self._keys`.
        if self._num_matrices > 0 and self._num_tables > 0:
            # Obtain the intersection of the sets of keys for all matrices
            mkeys = set(matrices[0].keys)
            for matrix in matrices[1:]:
                mkeys = mkeys.intersection(set(matrix.keys))

            # Obtain the intersection of the sets of keys for all tables
            tkeys = set(tables[0].keys())
            for table in tables[1:]:
                tkeys = tkeys.intersection(set(table.keys()))

            # Compute the intersection of the sets of keys
            self._keys = sorted(mkeys.intersection(tkeys))
        elif self._num_matrices > 0:
            # Obtain the intersection of the sets of keys for all matrices
            mkeys = set(matrices[0].keys)
            for matrix in matrices[1:]:
                mkeys = mkeys.intersection(set(matrix.keys))
            self._keys = sorted(mkeys)
        else:
            # Obtain the intersection of the sets of keys for all tables
            tkeys = set(tables[0].keys)
            for table in tables[1:]:
                tkeys = tkeys.intersection(set(table.keys()))
            self._keys = sorted(tkeys)

        # If `mweights` is not provided, default to 1.0 for every item in `matrices`
        if not mweights:
            self._mweights = [1.0] * len(matrices)
        else:
            self._mweights = mweights

        # If `tweights` is not provided, default to 1.0 for every item in `tables`
        if not tweights:
            self._tweights = [1.0] * len(tables)
        else:
            self._tweights = tweights

    @functools.cache  # Caching results for repeated calls with the same arguments
    def compute_score(self, sampled_taxa: Tuple, method: str) -> float:
        """
        Compute the score for a given set of sampled taxa.

        The score is computed as the weighted average of the distances between each pair of taxa,
        where the weights are given by the `weights` parameter.

        :param sampled_taxa: The set of sampled taxa.
        :param method: The method to use for computing the distance.
        :return: The score.
        """

        # Initialize the score
        score = 0.0

        # If there are matrices, compute their contribution to the score
        if self._matrices:
            # Compute the distances between each pair of taxa
            dists = [
                [
                    self._matrices[m][lang1, lang2]
                    for lang1, lang2 in itertools.combinations(sampled_taxa, 2)
                ]
                for m in range(self._num_matrices)
            ]

            # Calculate the distances
            calculated_dists = [
                sum(dist) / len(sampled_taxa) if method == "mean" else sum(dist)
                for dist in dists
            ]

            # Add to the score
            score += sum(
                dist * weight for dist, weight in zip(calculated_dists, self._mweights)
            )

        # If there are tables, compute their contribution to the score
        if self._tables:
            # Compute the values for each sampled taxon
            table_vals = [
                [self._tables[t][taxon] for taxon in sampled_taxa]
                for t in range(self._num_tables)
            ]

            # Calculate the values
            calculated_vals = [
                sum(val) / len(sampled_taxa) if method == "mean" else sum(val)
                for val in table_vals
            ]

            # Add to the score
            score += sum(
                val * weight for val, weight in zip(calculated_vals, self._tweights)
            )

        return score

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

                # Compute the score without the frequency part
                score = self.compute_score(sampled_taxa, method)

                # Compute the frequency of each taxon in the sample
                fdists = [taxa_counter[taxon] for taxon in sampled_taxa]

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

                # Add the fdist to the final score
                score += fdist * freq_weight
                candidates.append((sampled_taxa, score))

            # Perform a weighted sampling of the candidates
            candidates, weights = zip(*candidates)
            selected = random.choices(candidates, weights=weights, k=1)[0]

            # Update the tracking of the sampled languages across different sets
            taxa_counter.update(selected)

            # Yield the results
            yield selected

    def progressive_sample(
        self,
        k: int,
        t: int = 1,
        candidates_ratio: float = 0.1,
        freq_weight: float = 1.0,
        method: str = "mean",
        seed: Optional[Union[int, str]] = None,
    ):
        # Validate method
        if method not in ["mean", "sum"]:
            raise ValueError(
                "Invalid method in `progressive_sample()`: {}".format(method)
            )

        # Set the seed
        random.seed(seed)

        # Determine number of initial random candidates
        num_candidates = max(1, int(len(self._keys) * candidates_ratio))

        # Initialize the set of sampled sets
        taxa_counter = collections.Counter()

        for _ in range(t):
            # Start with a single random taxon
            sampled_taxa = [random.choice(self._keys)]

            while len(sampled_taxa) < k:
                # Randomly pick remaining taxa to be considered
                remaining_taxa = random.sample(
                    list(set(self._keys) - set(sampled_taxa)), num_candidates
                )

                scores = []
                for taxon in remaining_taxa:
                    new_sample = tuple(sampled_taxa + [taxon])
                    score = self.compute_score(new_sample, method)

                    # Penalize taxa that have been sampled already
                    fdists = [taxa_counter[taxon] for taxon in new_sample]
                    if not taxa_counter.most_common(1):
                        fdist = 1.0
                    else:
                        most_common_count = taxa_counter.most_common(1)[0][1]
                        fdist = 1 - (
                            sum(fdists)
                            / (
                                k * most_common_count
                                if method == "mean"
                                else sum([e[1] for e in taxa_counter.most_common(k)])
                            )
                        )
                    score += fdist * freq_weight

                    scores.append((new_sample, score))

                # Find the set with the highest score
                best_sample = max(scores, key=lambda x: x[1])[0]
                sampled_taxa = list(best_sample)

            # Update the tracking of the sampled languages across different sets
            taxa_counter.update(sampled_taxa)

            # Yield the results
            yield tuple(sampled_taxa)


class LanguageSampler(Sampler):
    def __init__(self):
        """
        Initialize the sampler.
        """

        # TODO: add rescaler
        from pathlib import Path

        p = Path(__file__).parent.parent.parent / "maja.csv"
        tok_weights = build_dict_from_file(p, "Glottocode", "motion hits")
        print(f"Read {len(tok_weights)} entries from token weights")

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
        super().__init__([phylomatrix, geomatrix], [tok_weights])
