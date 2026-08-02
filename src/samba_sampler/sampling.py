"""
Functions for obtaining taxa samples.
"""

# Import from standard library
from typing import Optional, Union, List, Tuple, Generator
import collections
import functools
import itertools
import logging
import random

# Import from other modules
from .common import ETC_PATH, DistanceMatrix, build_table_from_file


class Sampler:
    def __init__(
        self,
        matrices: Optional[List[DistanceMatrix]] = None,
        tables: Optional[List[dict]] = None,
        mweights: Optional[List[float]] = None,
        tweights: Optional[List[float]] = None,
        taxa_subset: Optional[List[str]] = None,
    ):
        """
        Initialize the sampler.

        @param matrices: A list of matrices to be used for sampling.
        @param tables: A list of tables to be used for sampling.
        @param mweights: A list of weights for the matrices.
        @param tweights: A list of weights for the tables.
        @param taxa_subset: A list of taxa to be used for sampling.
        """

        # Raise an error if no matrix and no table was provided
        if not matrices and not tables:
            raise ValueError("No matrix and no table provided.")

        # Store the matrices and tables internally, normalizing a missing
        # argument to an empty list so that either one may be omitted.
        self._matrices = matrices = list(matrices) if matrices else []
        self._num_matrices = len(self._matrices)

        self._tables = tables = list(tables) if tables else []
        self._num_tables = len(self._tables)

        # Obtain the set of common keys for all matrices and for all tables, if they
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
            tkeys = set(tables[0].keys())
            for table in tables[1:]:
                tkeys = tkeys.intersection(set(table.keys()))
            self._keys = sorted(tkeys)

        # If a subset of taxa was provided, filter the keys
        if taxa_subset:
            self._keys = [key for key in self._keys if key in taxa_subset]

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

    # Bounded cache for repeated calls with the same arguments. The score
    # depends only on the (immutable) matrices/tables and the taxa, so caching
    # is safe; the size cap keeps memory bounded during long sampling runs.
    @functools.lru_cache(maxsize=2**16)
    def compute_score(self, sampled_taxa: Tuple, method: str) -> float:
        """
        Compute the score for a given set of sampled taxa.

        The score is computed as the weighted average of the distances between each pair of taxa,
        where the weights are given by the `weights` parameter.

        @param sampled_taxa: The set of sampled taxa.
        @param method: The method to use for computing the distance.
        @return: The score.
        """

        # Initialize the score
        score = 0.0

        # If there are matrices, compute their contribution to the score
        if self._matrices:
            # Compute the distances between each pair of taxa
            dists = [
                [
                    self._matrices[m][taxon1, taxon2]
                    for taxon1, taxon2 in itertools.combinations(sampled_taxa, 2)
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
        n: int = 1,
        freq_weight: float = 1.0,
        candidates_ratio: float = 0.1,  # only used for progressive sampling
        algorithm: str = "standard",
        summary: str = "mean",
        subpop: int = 100,  # only used for standard sampling
        seed: Optional[Union[int, str]] = None,
    ) -> Generator[Tuple[str], None, None]:
        """
        Sample a set of taxa.
        """

        if algorithm == "standard":
            yield from self.standard_sampling(
                k, n, freq_weight, summary, subpop, seed=seed
            )
        elif algorithm == "progressive":
            yield from self.progressive_sampling(
                k, n, freq_weight, summary, candidates_ratio, seed=seed
            )  # TODO: check parameters
        else:
            raise ValueError(f"Invalid sampling method: {algorithm}")

    def standard_sampling(
        self,
        k: int,
        n: int = 1,
        freq_weight: float = 1.0,
        summary: str = "mean",
        subpop: int = 100,
        seed: Optional[Union[int, str]] = None,
    ) -> Generator[Tuple[str], None, None]:
        """
        Sample a set of taxa.

        This is the default sampling method, using a non-progressive approach in which full sets
        of taxa are sampled at random and picked based on their score. For a progressive approach,
        see `sample_progressive()`.

        @param k: The number of taxa to sample.
        @param n: The number of samples to generate.
        @param freq_weight: The weight to be given to the frequency of each taxon in the sample.
        @param summary: The method to use for computing the distance.
        @param subpop: The number of samples to generate for each subpopulation.
        @param seed: The seed to be used for the random number generator.
        @return: A generator of samples.
        """

        # Use a local random generator so sampling is reproducible without
        # mutating the global `random` state shared with the rest of the process.
        rng = random.Random(seed)

        # Validate method
        if summary not in ["mean", "sum"]:
            raise ValueError("Invalid method in `sample()`: {}".format(summary))

        # Initialize the set of sampled sets
        taxa_counter = collections.Counter()

        for _ in range(n):
            candidates = []
            for _ in range(subpop):
                # Obtain a random sample of k taxa
                sampled_taxa = tuple(rng.sample(self._keys, k))

                # Compute the score without the frequency part
                score = self.compute_score(sampled_taxa, summary)

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
                            if summary == "mean"
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
            selected = rng.choices(candidates, weights=weights, k=1)[0]

            # Update the tracking of the sampled languages across different sets
            taxa_counter.update(selected)

            # Yield the results
            yield selected

    def progressive_sampling(
        self,
        k: int,
        n: int = 1,
        freq_weight: float = 1.0,
        summary: str = "mean",
        candidates_ratio: float = 0.1,
        seed: Optional[Union[int, str]] = None,
    ) -> Generator[Tuple[str], None, None]:
        """
        Sample a set of taxa using a progressive approach.

        This method starts with a single random taxon and progressively adds new taxa to the
        sample, picking the one that maximizes the score. For a non-progressive approach, see
        `sample()`.

        @param k: The number of taxa to sample.
        @param n: The number of samples to generate.
        @param candidates_ratio: The ratio of candidates to consider for each new taxon.
        @param freq_weight: The weight to be given to the frequency of each taxon in the sample.
        @param summary: The method to use for computing the distance.
        @param seed: The seed to be used for the random number generator.
        @return: A generator of samples.
        """

        # Use a local random generator so sampling is reproducible without
        # mutating the global `random` state shared with the rest of the process.
        rng = random.Random(seed)

        # Validate method
        if summary not in ["mean", "sum"]:
            raise ValueError(
                "Invalid method in `progressive_sample()`: {}".format(summary)
            )

        # Determine number of initial random candidates
        num_candidates = max(1, int(len(self._keys) * candidates_ratio))

        # Initialize the set of sampled sets
        taxa_counter = collections.Counter()

        for _ in range(n):
            # Start with a single random taxon
            sampled_taxa = [rng.choice(self._keys)]

            while len(sampled_taxa) < k:
                # Randomly pick remaining taxa to be considered, capped at the
                # number of taxa still available to avoid a sample larger than
                # the population.
                available = sorted(set(self._keys) - set(sampled_taxa))
                remaining_taxa = rng.sample(
                    available, min(num_candidates, len(available))
                )

                scores = []
                for taxon in remaining_taxa:
                    new_sample = tuple(sampled_taxa + [taxon])
                    score = self.compute_score(new_sample, summary)

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
                                if summary == "mean"
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
    """
    A ready-to-use sampler over the bundled language data.

    It combines the bundled phylogenetic (GLED) and geographic (haversine)
    distance matrices, both rescaled to a common [0, 1] range. An optional
    table of per-language weights may be supplied to bias the sampling.
    """

    def __init__(
        self,
        matrix_weights: Optional[List[float]] = None,
        weights_table: Optional[dict] = None,
        table_weights: Optional[List[float]] = None,
    ):
        """
        Initialize the sampler.

        @param matrix_weights: Optional weights for the phylogenetic and
            geographic matrices, in that order. Defaults to equal weights.
        @param weights_table: Optional mapping of glottocode to a scalar weight
            to be added to the sampling score.
        @param table_weights: Optional weights for the supplied table.
        """

        phylomatrix = DistanceMatrix(filename=ETC_PATH / "gled.matrix.bz2")
        logging.info("Read %d entries from phylogenetic matrix", len(phylomatrix.data))

        geomatrix = DistanceMatrix(filename=ETC_PATH / "haversine.matrix.bz2")
        logging.info("Read %d entries from geographic matrix", len(geomatrix.data))

        # Rescale both matrices to a common range so that neither dominates the
        # score purely because of its native units.
        phylomatrix.rescale()
        geomatrix.rescale()

        tables = [weights_table] if weights_table else None

        # Initialize the parent class
        super().__init__(
            [phylomatrix, geomatrix],
            tables,
            mweights=matrix_weights,
            tweights=table_weights,
        )


class GenericSampler(Sampler):
    def __init__(
        self,
        matrix_files,
        table_files,
        matrix_weights,
        table_weights,
        taxa_subset=None,
        table_key: str = "Glottocode",
        table_value: str = "weight",
    ):
        """
        Initialize the sampler.

        @param matrix_files: Filenames of distance matrices to load from the
            bundled data directory.
        @param table_files: Filenames of tabular files to load as weight tables.
        @param matrix_weights: Weights for the matrices.
        @param table_weights: Weights for the tables.
        @param taxa_subset: Optional subset of taxa to restrict sampling to.
        @param table_key: Name of the column to use as the key when reading tables.
        @param table_value: Name of the column to use as the value when reading tables.
        """

        # Load all matrices and tables
        if matrix_files:
            matrices = [DistanceMatrix(filename=ETC_PATH / m) for m in matrix_files]
        else:
            matrices = []

        if table_files:
            tables = [
                build_table_from_file(t, table_key, table_value) for t in table_files
            ]
        else:
            tables = []

        # Initialize the parent class
        super().__init__(
            matrices,
            tables,
            mweights=matrix_weights,
            tweights=table_weights,
            taxa_subset=taxa_subset,
        )
