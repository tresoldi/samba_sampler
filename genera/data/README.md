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
| `confidence` | `wals-anchored` \| `gled-supported` \| `glottolog-only` \| `isolate` |

## `genera_summary.tsv` — one row per genus

`genus_id`, `genus_name`, `family`, `family_gc`, `n_languages`, `method`,
`wals_genus`, `gled_monophyletic` (blank if untestable), `gled_crown_height`
(GLED relative divergence), `example_languages`.

## `params.json`, `calibration_report.md`, `validation_report.md`

Calibrated parameters and the calibration/validation reports.

## Sources & attribution

Derived from **Glottolog** (classification; bundled with samba_sampler) and
**WALS** (Dryer genera; CLDF `languages.csv`, v2020.1) — both CC-BY 4.0 — plus
the bundled **GLED** phylogenies (validation only). The fetched WALS cache lives
under `data/sources/` and is not vendored here; regenerate it with
`pipeline/01_fetch_sources.py`.
