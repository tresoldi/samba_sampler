#!/usr/bin/env python
"""
Phase 2 — ASJP refinement of the genus assignment.

Phase 1 leaves a tier of `glottolog-only` genera: groups that rest on the
Glottolog classification alone, with no WALS labels and no testable GLED
monophyly. This script brings a second independent source — ASJP wordlists — to
bear on them:

  1. Upgrade `glottolog-only` genera to `asjp-supported` where ASJP wordlists
     recover them as a lexically cohesive, separated cluster within their family
     (see genlib.asjp_genus_support).
  2. Flag over-large genera and, where ASJP coverage allows, propose a data-driven
     subdivision (advisory; genus ids are left unchanged in this release).

Re-run `04_validate.py` afterwards to refresh the validation report.

Reads:
  data/genera.tsv              (phase-1 assignment)
  data/sources/asjp_*.csv      (fetched by 01_fetch_sources.py)
Writes:
  data/genera.tsv              (adds asjp_n / asjp_score, upgrades confidence)
  data/asjp_subdivisions.tsv   (over-large genus subdivision proposals)
  data/asjp_report.md
  data/params.json             (adds a "phase2" block)
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict

import numpy as np
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score

import genlib as gl

# ASJP-responsible confidence tiers: recomputed from scratch on every run so the
# script is idempotent (a genus that stops separating reverts to glottolog-only).
ASJP_TIERS = {"glottolog-only", "asjp-supported"}

LARGE_GENUS = 100  # size at/above which a genus is a subdivision candidate
MIN_SUB_ASJP = 8  # min ASJP-covered members needed to attempt a subdivision
MAX_SUB_ASJP = 250  # cap on members entering the clustering, to bound cost


def read_genera():
    path = gl.BASE / "data" / "genera.tsv"
    return list(csv.DictReader(path.open(encoding="utf-8"), delimiter="\t"))


def propose_subdivision(members_gcs, asjp):
    """
    Propose a subdivision of an over-large genus from ASJP distances.

    Clusters the genus's ASJP-covered members with average-linkage agglomerative
    clustering over the precomputed distance matrix and picks the number of
    clusters (2..7) with the best silhouette. Missing pairwise distances are
    imputed as the maximum normalized distance (1.0). Returns ``None`` when
    coverage is too thin or no split beats a single cluster.
    """
    covered = sorted(asjp.covered(members_gcs))[:MAX_SUB_ASJP]
    n = len(covered)
    if n < MIN_SUB_ASJP:
        return None

    dmat = np.ones((n, n), dtype=float)
    np.fill_diagonal(dmat, 0.0)
    finite = 0
    for i in range(n):
        for j in range(i + 1, n):
            d = asjp.distance(covered[i], covered[j])
            if d is not None:
                dmat[i, j] = dmat[j, i] = d
                finite += 1
    if finite < n:  # too sparse to cluster meaningfully
        return None

    best = None
    for k in range(2, min(7, n - 1) + 1):
        labels = AgglomerativeClustering(
            n_clusters=k, metric="precomputed", linkage="average"
        ).fit_predict(dmat)
        if len(set(labels)) < 2:
            continue
        score = silhouette_score(dmat, labels, metric="precomputed")
        if best is None or score > best[0]:
            best = (score, k, sorted(Counter(labels).values(), reverse=True))
    if best and best[0] > 0:
        return {"k": best[1], "silhouette": round(best[0], 4), "sizes": best[2], "n_asjp": n}
    return None


def main() -> int:
    rows = read_genera()
    glottolog = gl.Glottolog.load()

    genus_of = {r["glottocode"]: r["genus_id"] for r in rows}
    members = defaultdict(list)
    base_conf = {}
    for r in rows:
        members[r["genus_id"]].append(r["glottocode"])
        base_conf.setdefault(r["genus_id"], r["confidence"])

    target_genera = [g for g, c in base_conf.items() if c in ASJP_TIERS]
    print(f"ASJP-responsible (glottolog-only) genera to test: {len(target_genera)}")

    print("Loading ASJP wordlists (in-scope languages only) ...")
    asjp = gl.Asjp.load(keep=set(genus_of))
    print(f"ASJP wordlists cover {len(asjp.words)} in-scope languages")

    support = gl.asjp_genus_support(genus_of, glottolog, asjp, target_genera)
    upgraded = {g for g, s in support.items() if s.get("supported")}
    tested = {g for g, s in support.items() if "score" in s}
    print(f"testable: {len(tested)}  upgraded -> asjp-supported: {len(upgraded)}")

    # ----- rewrite genera.tsv: add asjp columns, refresh confidence -----
    fields = list(rows[0].keys())
    for extra in ("asjp_n", "asjp_score"):
        if extra not in fields:
            fields.append(extra)
    lang_path = gl.BASE / "data" / "genera.tsv"
    with lang_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        for r in rows:
            g = r["genus_id"]
            s = support.get(g, {})
            r["asjp_n"] = s.get("n_asjp", "")
            r["asjp_score"] = s.get("score", "")
            if g in target_genera:
                r["confidence"] = "asjp-supported" if g in upgraded else "glottolog-only"
            w.writerow(r)
    print(f"Wrote {lang_path}")

    # ----- over-large genus subdivision proposals -----
    large = sorted(
        ((g, len(m)) for g, m in members.items() if len(m) >= LARGE_GENUS),
        key=lambda x: -x[1],
    )
    proposals = []
    sub_path = gl.BASE / "data" / "asjp_subdivisions.tsv"
    with sub_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(
            ["genus_id", "genus_name", "family", "n_languages", "n_asjp",
             "proposed_k", "silhouette", "cluster_sizes"]
        )
        for g, n in large:
            fam = glottolog.family_name(glottolog.family_gc(members[g][0]))
            prop = propose_subdivision(members[g], asjp)
            if prop:
                proposals.append((g, n, fam, prop))
                w.writerow([g, glottolog.name(g), fam, n, prop["n_asjp"],
                            prop["k"], prop["silhouette"],
                            ";".join(map(str, prop["sizes"]))])
            else:
                w.writerow([g, glottolog.name(g), fam, n,
                            len(asjp.covered(members[g])), "", "", ""])
    print(f"Wrote {sub_path} ({len(large)} over-large genera, "
          f"{len(proposals)} with a proposed split)")

    # ----- report -----
    n_tested = len(tested)
    untestable = len(target_genera) - n_tested
    lines = ["# ASJP refinement (phase 2)\n"]
    lines.append("## Upgrading `glottolog-only` genera\n")
    lines.append(f"- Target `glottolog-only` genera: **{len(target_genera)}**.")
    lines.append(f"- Testable with ASJP (>=2 comparable wordlists inside the genus "
                 f"and >=1 outside it): **{n_tested}**.")
    lines.append(f"- **Upgraded to `asjp-supported`: {len(upgraded)}** — mean "
                 f"intra-genus ASJP distance below the distance to nearest "
                 f"same-family relatives.")
    lines.append(f"- Not upgraded: {n_tested - len(upgraded)} tested but not "
                 f"separated; {untestable} untestable (insufficient ASJP "
                 f"coverage).\n")

    best_up = sorted(
        ((g, support[g]) for g in upgraded), key=lambda kv: -kv[1]["score"]
    )[:12]
    if best_up:
        lines.append("### Best-separated upgraded genera\n")
        lines.append("| genus | family | n(ASJP) | intra | inter | score |")
        lines.append("|-------|--------|---------|-------|-------|-------|")
        for g, s in best_up:
            fam = glottolog.family_name(glottolog.family_gc(members[g][0]))
            lines.append(f"| {glottolog.name(g)} | {fam} | {s['n_asjp']} | "
                         f"{s['intra']} | {s['inter']} | {s['score']} |")
        lines.append("")

    lines.append("## Over-large genus subdivision (advisory)\n")
    lines.append(f"- Genera with >= {LARGE_GENUS} languages: **{len(large)}**; "
                 f"ASJP proposes a subdivision (positive silhouette) for "
                 f"**{len(proposals)}** of them.")
    lines.append("- Proposals are advisory: `genus_id`s are unchanged in this "
                 "release (subdividing would coin new genus nodes, left to future "
                 "work). Full detail in `data/asjp_subdivisions.tsv`.\n")
    if proposals:
        lines.append("| genus | family | languages | proposed k | silhouette | cluster sizes |")
        lines.append("|-------|--------|-----------|-----------|-----------|---------------|")
        for g, n, fam, prop in proposals[:15]:
            lines.append(f"| {glottolog.name(g)} | {fam} | {n} | {prop['k']} | "
                         f"{prop['silhouette']} | {';'.join(map(str, prop['sizes']))} |")
        lines.append("")

    (gl.BASE / "data" / "asjp_report.md").write_text("\n".join(lines))
    print("Wrote data/asjp_report.md")

    # ----- params.json phase2 block -----
    params_path = gl.BASE / "data" / "params.json"
    params = json.loads(params_path.read_text())
    params["phase2"] = {
        "asjp_ref": "012795349540ba0dabfdcf2be16f2e77622f62d6",
        "asjp_languages_in_scope": len(asjp.words),
        "glottolog_only_targets": len(target_genera),
        "asjp_testable": n_tested,
        "upgraded_to_asjp_supported": len(upgraded),
        "large_genera": len(large),
        "subdivision_proposals": len(proposals),
    }
    params_path.write_text(json.dumps(params, indent=2))
    print("Updated data/params.json (phase2 block)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
