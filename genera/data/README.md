# genera data

Generated artifacts. Rebuild with the `pipeline/` scripts (see `../README.md`).

## `genera.tsv` — one row per in-scope Glottolog language

| column | description |
|--------|-------------|
| `glottocode` | Glottolog glottocode of the language |
| `name` | Glottolog language name |
| `family` | Glottolog family name |
| `family_gc` | glottocode of the family root (isolates: the language itself) |
| `genus_id` | glottocode of the Glottolog node defining the genus (stable id) |
| `genus_name` | Glottolog name of that node |
| `method` | `glottolog_depth` or `isolate` |
| `wals_genus` | majority WALS/Dryer genus among the genus's members, if any |
| `confidence` | `wals-anchored` \| `gled-supported` \| `asjp-supported` \| `glottolog-only` \| `isolate` |
| `asjp_n` | phase 2: ASJP-covered members compared for this genus (blank if untested) |
| `asjp_score` | phase 2: ASJP separation score (intra vs nearest-relative distance) |

## `genera_summary.tsv` — one row per genus

`genus_id`, `genus_name`, `family`, `family_gc`, `n_languages`, `method`,
`wals_genus`, `gled_monophyletic` (blank if untestable), `gled_crown_height`
(GLED relative divergence), `example_languages`.

## `asjp_report.md`, `asjp_subdivisions.tsv` (phase 2)

`asjp_report.md` summarizes the ASJP upgrades of `glottolog-only` genera and the
over-large-genus subdivision proposals. `asjp_subdivisions.tsv` has one row per
over-large genus: `genus_id`, `genus_name`, `family`, `n_languages`, `n_asjp`,
`proposed_k`, `silhouette`, `cluster_sizes` (advisory; genus ids are unchanged).

## `params.json`, `calibration_report.md`, `validation_report.md`

Calibrated parameters (with a `phase2` block) and the calibration/validation
reports.

## Sources & attribution

Derived from **Glottolog** (classification; bundled with samba_sampler) and
**WALS** (Dryer genera; CLDF `languages.csv`, v2020.1) — both CC-BY 4.0 — plus
the bundled **GLED** phylogenies (validation) and, in phase 2, **ASJP**
wordlists (lexibank/asjp CLDF, pinned commit; CC-BY 4.0) for independent lexical
support. The fetched WALS and ASJP caches live under `data/sources/` and are not
vendored here; regenerate them with `pipeline/01_fetch_sources.py`.
