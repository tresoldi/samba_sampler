# Language Genera — project plan

A side project of `samba_sampler` that extends Dryer's (1989) notion of the
language **genus** — the mid-level genealogical unit used in WALS to control for
genealogical bias — with (1) a more operational, reproducible definition and
(2) a mapping of **every** Glottolog language to a genus, coining names where
none exist.

## Motivation

Dryer's genus is a genealogical grouping of roughly comparable time depth — on
the order of the subfamilies of Indo-European (Germanic, Romance, Slavic,
Indo-Aryan, Iranian, …), i.e. a divergence of no more than ~3,500–4,000 years.
It exists so that typological samples can be balanced by counting units of
comparable age rather than families, whose depths vary by an order of magnitude.

The WALS/Dryer genera have two limitations we address here:

1. **Impressionistic definition.** They were assigned by expert judgement, not
   from an explicit, reproducible criterion.
2. **Partial coverage.** They cover only the ~2,600 WALS languages (~500
   genera). Most of Glottolog's ~8,600 languages have no genus.

## What "scientific definition" means here

> A **genus** is a maximal clade whose crown divergence does not exceed a fixed
> threshold **T**, where T is calibrated so the induced partition best
> reproduces (a) Indo-European's traditional subfamilies and (b) WALS genera on
> the families for which we have both a phylogeny and WALS labels.

Because our phylogenies (GLED) carry **relative** divergence (accumulated
lexical change), T is a calibrated *distance*, not an absolute age. We report it
as such and record the Indo-European yardstick explicitly.

## Coverage strategy (hybrid)

Only ~2,770 of 8,591 Glottolog language-level nodes sit in a GLED tree, so a
single criterion cannot cover everything. We therefore combine, in decreasing
confidence:

| method       | applies to                                   | confidence |
|--------------|----------------------------------------------|------------|
| `tree`       | families with a GLED tree, cut at T          | high       |
| `asjp`       | GLED-less families with enough ASJP wordlists| medium     |
| `glottolog`  | everything else, via a calibrated depth rule | lower      |
| `isolate`    | Glottolog isolates                           | n/a (1 lang)|

The Glottolog-depth rule is itself *calibrated* against the tree-based genera:
we measure where tree-based genus boundaries fall in Glottolog's nesting and
learn the depth that best matches, so the fallback stays consistent with the
scientific cut instead of being an independent guess.

## Units in scope

`level == language`, excluding non-genetic categories: Bookkeeping, Artificial
Language, Sign Language, Unattested, Unclassifiable, Speech Register, Mixed
Language, Pidgin. Genuinely extinct-but-attested languages (Latin, Sanskrit, …)
are kept. Isolates each form their own genus (Dryer's rule).

## Naming coined genera

1. Reuse the WALS genus name when the induced genus matches a WALS genus.
2. Otherwise the Glottolog subgroup name spanning the genus, + disambiguator.
3. Otherwise a representative (most-documented) member language.
4. Last resort: `<family>-genus-NNN`.

Every genus also carries a stable, name-independent `genus_id`.

## Sources

- **Glottolog** dump (bundled `src/samba_sampler/etc/glottolog.*.tsv`) — the
  classification backbone and language metadata.
- **GLED** family trees (bundled `scaffolding/gled.trees/*.tree` and
  `global_tree.gled.newick`) — relative-divergence phylogenies.
- **WALS** CLDF `languages.csv` (fetched, cached under `data/sources/`) —
  Dryer genera as calibration + validation ground truth.
- **ASJP** (lexibank, fetched) — wordlists for GLED-less families (phase 2).

## Build order

- **Phase 1 (bundled + WALS):** scaffolding, WALS fetch, calibrate T + the
  Glottolog-depth rule, assign genera to *all* in-scope languages via
  tree+glottolog+isolate, validate against WALS. Produces a complete artifact.
- **Phase 2 (ASJP):** upgrade GLED-less families from `glottolog` to `asjp`
  where wordlists allow; re-validate.

## Outputs

- `data/genera.tsv` — one row per language: `glottocode, name, family,
  genus_id, genus_name, method, confidence`.
- `data/genera_summary.tsv` — one row per genus.
- `data/calibration_report.md` — chosen T, agreement metrics, IE case study.
- `methods.md` — the definition, calibration, and caveats, prose form.

## Explicitly out of scope (for now)

No changes to the shipped `samba_sampler` package or its data. Wiring genera
into the sampler (a genus-aware sampling mode) is a documented future step.
