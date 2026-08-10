# genera

Extending Dryer's (1989) language **genus** with an operational, reproducible
definition and full Glottolog coverage. See [`PLAN.md`](PLAN.md) for the design
and [`methods.md`](methods.md) for the definition and calibration.

This is a **side project** of `samba_sampler`. It does not modify the shipped
package or its data; it reads the bundled Glottolog dump and GLED trees and
writes its own artifacts under `data/`.

## Reproduce

```bash
pip install -r requirements.txt
python pipeline/01_fetch_sources.py     # cache WALS (+ ASJP in phase 2)
python pipeline/02_calibrate.py         # calibrate cut depth d against WALS
python pipeline/03_assign_genera.py     # build data/genera.tsv
python pipeline/04_validate.py          # agreement vs WALS, size dist, naming
```

## Method in one line

A genus is a Glottolog clade cut at a WALS-calibrated taxonomic **depth d = 3**;
GLED phylogenies validate the result and drive phase-2 refinement. See
[`methods.md`](methods.md) for why a GLED distance threshold was rejected.

## Result (phase 1)

- **7,650 languages → 1,551 genera**, full coverage.
- **ARI 0.867 / V 0.925** against WALS/Dryer genera; Indo-European recovered
  as 10 genera vs WALS's 11.
- **95.5%** of GLED-testable genera are monophyletic.

## Outputs

| file | contents |
|------|----------|
| `data/genera.tsv` | one row per in-scope language → its genus (see `data/README.md`) |
| `data/genera_summary.tsv` | one row per genus |
| `data/calibration_report.md` | depth sweep, IE yardstick, GLED validation |
| `data/validation_report.md` | WALS agreement, size distribution, naming |
| `data/params.json` | calibrated parameters |

## Status

- Phase 1 (GLED + Glottolog + WALS) — **complete**.
- Phase 2 (ASJP extension; subdivide over-large genera) — planned.
