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
python pipeline/02_calibrate.py         # calibrate threshold T + depth rule
python pipeline/03_assign_genera.py     # build data/genera.tsv
python pipeline/04_validate.py          # agreement vs WALS, IE case study
```

## Outputs

| file | contents |
|------|----------|
| `data/genera.tsv` | one row per in-scope language → its genus |
| `data/genera_summary.tsv` | one row per genus |
| `data/calibration_report.md` | chosen T, metrics, Indo-European check |

## Status

Phase 1 (GLED + Glottolog + WALS) — in progress.
Phase 2 (ASJP extension) — planned.
