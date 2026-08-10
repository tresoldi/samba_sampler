# Methods — defining and assigning language genera

> Status: skeleton. The calibrated threshold **T**, the Glottolog-depth rule,
> and the agreement metrics are filled in by `pipeline/02_calibrate.py` and
> written to `data/calibration_report.md`; this document narrates the final,
> calibrated method once phase 1 runs.

## 1. The unit

A **genus** is a maximal clade whose crown divergence does not exceed a fixed
threshold **T**. It approximates Dryer's (1989) genus: a genealogical group of
roughly the time depth of the subfamilies of Indo-European.

## 2. Scope

Glottolog nodes with `level == language`, excluding the non-genetic categories
Bookkeeping, Artificial Language, Sign Language, Unattested, Unclassifiable,
Speech Register, Mixed Language, and Pidgin. Isolates each form a singleton
genus.

## 3. Calibrating T

Threshold T is chosen on the families that have **both** a GLED phylogeny and
WALS genus labels, as the tree-cut height that maximises agreement (adjusted
Rand index / V-measure) between the induced partition and the WALS partition.
Indo-European is the explicit yardstick: at the chosen T, Germanic, Romance,
Slavic, Indic, Iranian, Baltic, Celtic, Greek, Armenian, Albanian and Nuristani
must each emerge as separate genera. GLED branch lengths encode *relative*
lexical divergence, so T is a distance, not an age.

## 4. Assignment (hybrid, by decreasing confidence)

1. **`tree`** — families with a GLED tree: cut at T.
2. **`asjp`** — GLED-less families with enough ASJP wordlists: cut an ASJP tree
   at a threshold cross-calibrated to T (phase 2).
3. **`glottolog`** — everything else: cut Glottolog's nesting at the depth that,
   on tree-covered families, best reproduces the tree-based genus boundaries.
4. **`isolate`** — Glottolog isolates: one genus each.

## 5. Naming

WALS name on match → spanning Glottolog subgroup name + disambiguator →
representative member language → `<family>-genus-NNN`. A stable `genus_id`
accompanies every genus regardless of name.

## 6. Caveats

- T is a relative-divergence distance, not a calibrated age.
- `glottolog`-method genera inherit Glottolog's uneven subclassification depth
  and are the lowest-confidence tier.
- Coverage and per-method counts are reported in `data/calibration_report.md`.
