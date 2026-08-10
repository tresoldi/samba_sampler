# Methods — defining and assigning language genera

This document states the operational genus definition and the calibration behind
it. Numbers are produced by `pipeline/02_calibrate.py` and
`pipeline/04_validate.py` and mirrored in `data/calibration_report.md` and
`data/validation_report.md`.

## 1. The unit

A **genus** approximates Dryer's (1989) genus: a genealogical group of roughly
the time depth of the subfamilies of Indo-European. Operationally:

> A **genus** is the Glottolog clade obtained by cutting a family's
> classification at a fixed taxonomic depth **d** below the family root, where d
> is calibrated to best reproduce WALS/Dryer genera. Shallow families are cut at
> their deepest available group (never split into individual languages);
> isolates each form one genus.

The genus is thus a *named Glottolog node*: its stable id is that node's
glottocode and its name is the node's Glottolog name.

## 2. Why not a GLED branch-length threshold?

The project began intending to cut GLED phylogenies at a divergence threshold T.
The data refuted this as a *primary* criterion: GLED branch lengths encode
relative lexical divergence that is **not on a common scale across families**.
Each family's own agreement-optimal cut height varies across roughly three
orders of magnitude, so no single T reproduces genus structure — a global T that
splits Germanic correctly merges Romance, Slavic, Celtic, Greek and Armenian into
one blob. A scale-free taxonomic-depth cut avoids this entirely.

Concretely, thresholding GLED distances tops out at ARI ≈ 0.48 against WALS,
while the depth cut reaches ARI ≈ 0.87 (see §4).

## 3. Scope

Glottolog nodes with `level == language`, excluding non-genetic categories
(Bookkeeping, Artificial Language, Sign Language, Unattested, Unclassifiable,
Speech Register, Mixed Language, Pidgin): **7,650 languages**. Genuinely
extinct-but-attested languages (Latin, Sanskrit, …) are kept.

## 4. Calibration

d is chosen by maximizing the adjusted Rand index between the induced partition
and WALS genera over the 2,261 in-scope languages carrying a WALS genus.

| depth d | ARI vs WALS | V-measure |
|---------|-------------|-----------|
| 1 | 0.447 | 0.899 |
| 2 | 0.773 | 0.931 |
| **3** | **0.867** | 0.925 |
| 4 | 0.650 | 0.896 |
| 5 | 0.507 | 0.879 |

**d = 3** is selected. Indo-European then resolves into 10 genera (Indo-Aryan,
Iranian, Nuristani, Slavic, Eastern Baltic, Nuclear Celtic, Northwest Germanic,
Latino-Faliscan, Greek, Armenian) against WALS's 11 — the intended yardstick.

## 5. Independent validation with GLED

GLED phylogenies are not used to define genera, but to check them. Of the 512
genera testable in GLED trees (≥2 leaves), **95.5% are monophyletic** — strong
independent support that the depth cut carves real clades. GLED crown depths also
quantify how *unequal* the genera remain (median 0.36, long right tail): the
residual that phase-2 refinement targets.

## 6. Assignment and confidence

Every in-scope language is assigned to its depth-d genus node. Each carries a
confidence tier:

- `wals-anchored` — the genus contains WALS-labeled languages (majority WALS
  genus recorded).
- `gled-supported` — no WALS labels, but the genus is monophyletic in GLED.
- `glottolog-only` — neither; rests on Glottolog classification alone.
- `isolate` — a Glottolog isolate (one language).

## 7. Naming

A genus is a Glottolog node, so its name is that node's Glottolog name (e.g.
*Indo-Aryan*, *Bantoid*, *Zapotecan*) and its id is the node's glottocode.
Where a genus aligns with a WALS genus, the WALS name is recorded alongside for
cross-reference. No opaque `family-genus-NNN` placeholders are needed.

## 8. Known limitations

- The cut inherits Glottolog's uneven classification depth. A few genera are very
  large (Bantoid 730, Oceanic 521, Indo-Aryan 220) — faithful to Dryer, whose own
  genera are equally uneven, but not equal-depth. Subdividing these into more
  uniform units with GLED/ASJP is the phase-2 improvement over Dryer.
- d is a single global depth; a minority of families would be better cut one
  level shallower or deeper. Per-family depth calibration is left to future work.
- WALS genus agreement (ARI 0.87) is high but not perfect; ~40% of WALS genera
  are split across our (usually finer) genera.

## 9. Phase 2 — ASJP support for `glottolog-only` genera

Phase 1 leaves 687 languages in `glottolog-only` genera: groups with no WALS
labels and no testable GLED monophyly. Phase 2 brings a second, independent
source — **ASJP** 40-item wordlists — to bear on them.

**Distance.** Language-pair distance is the mean, over shared concepts, of the
normalized Levenshtein distance (LDN) between ASJP-transcribed forms (minimum
over multiple forms per concept), requiring ≥10 shared concepts.

**Criterion.** A `glottolog-only` genus is upgraded to `asjp-supported` when it
is *lexically cohesive and separated*: its mean intra-genus ASJP distance is
smaller than the mean distance from its members to their **nearest same-family
relatives** (a positive silhouette-style separation). This is a *relative*
criterion, so — like the depth cut, and unlike a global distance threshold — it
does not assume ASJP divergences are on a common scale across families
(the same objection that ruled out a global GLED threshold in §2).

**Result.** Of the 394 `glottolog-only` genera, 70 have enough ASJP coverage to
test (≥2 comparable wordlists inside the genus and ≥1 outside); **63 separate and
are upgraded** (187 languages), leaving 500 `glottolog-only`. Examples span
families: Chug-Lish (Sino-Tibetan), Kolami-Naiki (Dravidian), Subanen and
Sambalic (Austronesian), Albanian-Tosk (Indo-European).

**Over-large genera.** The uneven-depth genera flagged in §8 (Bantoid, Oceanic,
Indo-Aryan, …) are clustered with average-linkage over their ASJP distances; the
best-silhouette split (k and cluster sizes) is recorded in
`data/asjp_subdivisions.tsv` as an **advisory** proposal — genus ids are left
unchanged, since coining subgenera is a larger step left to future work.
