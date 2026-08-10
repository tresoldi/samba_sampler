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

> A **genus** is the Glottolog clade obtained by cutting a family's
> classification at a calibrated taxonomic **depth d** below the family root,
> where d is chosen to best reproduce WALS/Dryer genera (with Indo-European as
> the yardstick).

**Note — method pivot (evidence-driven).** The project began intending to cut
GLED phylogenies at a divergence threshold T (see `methods.md` §2). Calibration
showed GLED branch lengths are not on a common scale across families
(per-family optimal cut heights span ~3 orders of magnitude), so no single T
reproduces genera — it tops out at ARI ≈ 0.48 vs WALS. A **scale-free Glottolog
taxonomic-depth cut reaches ARI ≈ 0.87** and recovers Indo-European's genera
almost exactly, so it is the primary criterion. GLED is retained for
independent validation and phase-2 refinement. See `methods.md` for the full
argument and numbers.

## Coverage strategy (hybrid)

The Glottolog taxonomic-depth cut covers **all** in-scope languages uniformly
(it needs only the classification, which Glottolog has for everything). Each
assignment carries a confidence tier recording how much independent evidence
backs it:

| confidence       | meaning                                             | phase-1 count |
|------------------|-----------------------------------------------------|---------------|
| `wals-anchored`  | genus contains WALS-labeled languages               | 6,502 |
| `gled-supported` | no WALS, but the genus is monophyletic in GLED      | 277 |
| `asjp-supported` | no WALS/GLED, but ASJP recovers it as a cluster (ph2)| 187 |
| `glottolog-only` | rests on Glottolog classification alone             | 500 |
| `isolate`        | Glottolog isolate (one language)                    | 184 |

Phase 2 (done) upgrades `glottolog-only` genera using ASJP wordlists and proposes
subdivisions for the few over-large genera (Bantoid, Oceanic, …).

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
  where wordlists allow; re-validate. **Done** — `pipeline/05_refine_asjp.py`
  upgrades 63 genera (187 languages) to a new `asjp-supported` tier and proposes
  ASJP-driven subdivisions for the over-large genera.

## Outputs

- `data/genera.tsv` — one row per language: `glottocode, name, family,
  genus_id, genus_name, method, confidence`.
- `data/genera_summary.tsv` — one row per genus.
- `data/calibration_report.md` — chosen T, agreement metrics, IE case study.
- `methods.md` — the definition, calibration, and caveats, prose form.

## Wiring into samba_sampler

The genus map is bundled into the package (`src/samba_sampler/etc/genera.tsv`)
and exposed as a genus-aware sampling mode: `samba.GenusSampler` and an optional
group-diversity term on the base `Sampler` reward samples that span distinct
genera, directly targeting genealogical autocorrelation. `samba.load_genera()`
returns the mapping for custom use.
