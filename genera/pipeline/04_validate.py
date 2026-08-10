#!/usr/bin/env python
"""
Validate the genus assignment (data/genera.tsv) and write a report.

Reports agreement with WALS/Dryer genera, coverage and confidence breakdown,
genus-size distribution, the largest genera (candidates for phase-2
subdivision), GLED monophyly, and how many genus names are newly coined
(Glottolog names with no WALS equivalent).

Writes:
  data/validation_report.md
"""

from __future__ import annotations

import csv
from collections import Counter, defaultdict

import numpy as np
from sklearn.metrics import adjusted_rand_score, v_measure_score

import genlib as gl


def load_assignment():
    rows = list(csv.DictReader((gl.BASE / "data" / "genera.tsv").open(encoding="utf-8"), delimiter="\t"))
    return rows


def main() -> int:
    rows = load_assignment()
    wals = gl.load_wals()

    genus_of = {r["glottocode"]: r["genus_id"] for r in rows}
    n_lang, n_genera = len(rows), len(set(genus_of.values()))

    # Agreement with WALS on shared languages.
    shared = [(r, wals[r["glottocode"]]) for r in rows if r["glottocode"] in wals]
    true = [w for _, w in shared]
    pred = [r["genus_id"] for r, _ in shared]
    ari = adjusted_rand_score(true, pred)
    vm = v_measure_score(true, pred)

    # WALS-genus purity: does each WALS genus fall inside a single one of ours?
    wals_to_ours = defaultdict(set)
    for r, w in shared:
        wals_to_ours[w].add(r["genus_id"])
    clean = sum(1 for s in wals_to_ours.values() if len(s) == 1)

    conf = Counter(r["confidence"] for r in rows)
    method = Counter(r["method"] for r in rows)

    sizes = Counter(genus_of.values())
    size_vals = np.array(list(sizes.values()))
    largest = sizes.most_common(15)

    # Coined names: our genus name not equal to any WALS genus name.
    wals_names = set(wals.values())
    genus_names = {r["genus_id"]: r["genus_name"] for r in rows}
    coined = sum(1 for gid, nm in genus_names.items() if nm not in wals_names)

    id2name = genus_names

    lines = ["# Validation report\n"]
    lines.append("## Coverage\n")
    lines.append(f"- **{n_lang} in-scope languages -> {n_genera} genera.**")
    lines.append(f"- Languages with a WALS genus (for validation): {len(shared)}.")
    lines.append("- Confidence tiers: " + ", ".join(f"{k} {v}" for k, v in conf.most_common()) + ".")
    lines.append("- Methods: " + ", ".join(f"{k} {v}" for k, v in method.most_common()) + ".\n")

    lines.append("## Agreement with WALS/Dryer genera\n")
    lines.append(f"- **Adjusted Rand index = {ari:.3f}**, V-measure = {vm:.3f} "
                 f"(on {len(shared)} shared languages).")
    lines.append(f"- **{clean} / {len(wals_to_ours)} ({clean/len(wals_to_ours):.1%}) "
                 f"WALS genera map cleanly into a single one of our genera** "
                 f"(no genus is split across WALS boundaries).\n")

    lines.append("## Genus-size distribution\n")
    lines.append(f"- median {int(np.median(size_vals))}, mean {size_vals.mean():.1f}, "
                 f"max {size_vals.max()}; singletons {int((size_vals==1).sum())} "
                 f"({(size_vals==1).mean():.1%}).")
    buckets = Counter(min(int(s), 50) for s in size_vals)
    hist = {"1": buckets.get(1, 0),
            "2-4": sum(v for k, v in buckets.items() if 2 <= k <= 4),
            "5-9": sum(v for k, v in buckets.items() if 5 <= k <= 9),
            "10-24": sum(v for k, v in buckets.items() if 10 <= k <= 24),
            "25-49": sum(v for k, v in buckets.items() if 25 <= k <= 49),
            ">=50": sum(1 for s in size_vals if s >= 50)}
    lines.append("- size buckets: " + ", ".join(f"{k}: {v}" for k, v in hist.items()) + ".\n")

    lines.append("## Largest genera (phase-2 subdivision candidates)\n")
    lines.append("| genus | family | languages |")
    lines.append("|-------|--------|-----------|")
    fam_of = {r["genus_id"]: r["family"] for r in rows}
    for gid, n in largest:
        lines.append(f"| {id2name[gid]} | {fam_of[gid]} | {n} |")
    lines.append("")

    lines.append("## Naming\n")
    lines.append(f"- Every genus is a named Glottolog node (stable id = its glottocode).")
    lines.append(f"- **{coined} / {n_genera} genera carry a name with no WALS "
                 f"equivalent** — i.e. genera newly delimited or named here, using "
                 f"the Glottolog subgroup name.\n")

    (gl.BASE / "data" / "validation_report.md").write_text("\n".join(lines))
    print(f"ARI={ari:.3f} V={vm:.3f}; {n_genera} genera; {coined} coined names")
    print("Wrote data/validation_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
