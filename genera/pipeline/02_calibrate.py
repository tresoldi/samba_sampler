#!/usr/bin/env python
"""
Calibrate the genus definition.

Primary criterion: a genus is a Glottolog clade obtained by cutting each
family's classification at a fixed taxonomic depth d below the family root. d is
chosen to maximize agreement (adjusted Rand index) with WALS/Dryer genera.

Why not a GLED branch-length threshold? GLED divergences are not on a common
scale across families (per-family optimal thresholds span ~1000x), so no single
distance cut reproduces genus structure; a scale-free taxonomic-depth cut does.
This script quantifies both facts.

GLED trees are used here to VALIDATE the resulting genera (monophyly and
relative crown depth), not to define them.

Writes:
  data/calibration_report.md
  data/params.json   {"method": "glottolog_depth", "depth": d, ...}
"""

from __future__ import annotations

import json
from collections import defaultdict

import numpy as np
from sklearn.metrics import adjusted_rand_score, v_measure_score

import genlib as gl

DEPTHS = list(range(1, 7))


def wals_eval_set(glottolog, wals):
    inscope = {r["glottocode"] for r in glottolog.in_scope_languages()}
    return {gc: g for gc, g in wals.items() if gc in inscope and gc in glottolog.by_gc}


def depth_agreement(glottolog, wl, depth):
    true = [wl[gc] for gc in wl]
    pred = [glottolog.genus_node(gc, depth) for gc in wl]
    return adjusted_rand_score(true, pred), v_measure_score(true, pred)


def branch_scale_heterogeneity(glottolog, wals, trees):
    """Spread of each family's own ARI-optimal GLED cut height (the anti-motivation)."""
    grid = None
    opt = []
    cal = [(t, {gc: wals[gc] for gc in t.leaf_gcs if gc in wals}) for t in trees]
    cal = [(t, l) for t, l in cal if len(set(l.values())) >= 2 and len(l) >= 4]
    all_h = np.array([h for t, _ in cal for h in t.internal_heights()])
    grid = np.quantile(all_h, np.linspace(0.02, 0.98, 60))
    for t, l in cal:
        best = max(
            grid,
            key=lambda T: adjusted_rand_score(
                [l[g] for g in l], [t.cut(T).get(g, -1) for g in l]
            ),
        )
        opt.append(float(best))
    opt = np.array(opt)
    return {
        "families": len(opt),
        "median": float(np.median(opt)),
        "q25": float(np.percentile(opt, 25)),
        "q75": float(np.percentile(opt, 75)),
        "min": float(opt.min()),
        "max": float(opt.max()),
    }


def indo_european_case(glottolog, wals, depth):
    ie_family = next(
        (
            r["glottocode"]
            for r in glottolog.rows
            if r["name"] == "Indo-European" and r["level"] == "family"
        ),
        None,
    )
    clusters = defaultdict(list)
    for gc in wals:
        if (
            gc in glottolog.by_gc
            and glottolog.family_name(glottolog.family_gc(gc)) == "Indo-European"
        ):
            clusters[glottolog.genus_node(gc, depth)].append(gc)
    rows = []
    for node, members in sorted(clusters.items(), key=lambda kv: -len(kv[1])):
        wals_genera = sorted({wals[g] for g in members})
        rows.append((glottolog.name(node), len(members), wals_genera))
    return rows


def main() -> int:
    print("Loading Glottolog, WALS, GLED trees ...")
    glottolog = gl.Glottolog.load()
    wals = gl.load_wals()
    trees = gl.load_gled_trees(glottolog)
    wl = wals_eval_set(glottolog, wals)
    print(f"WALS-labeled in-scope languages: {len(wl)}")

    sweep = [(d, *depth_agreement(glottolog, wl, d)) for d in DEPTHS]
    best_d, best_ari, best_vm = max(sweep, key=lambda r: r[1])
    print(f"Best depth d = {best_d}  (ARI={best_ari:.3f}, V={best_vm:.3f})")

    # Assign at best depth and validate against GLED.
    genus_of = gl.assign_genera(glottolog, best_d)
    n_genera = len(set(genus_of.values()))
    mono, testable, crown = gl.gled_monophyly(genus_of, trees)
    ie_rows = indo_european_case(glottolog, wals, best_d)
    scale = branch_scale_heterogeneity(glottolog, wals, trees)

    params = {
        "method": "glottolog_depth",
        "depth": best_d,
        "wals_agreement": {"ari": round(best_ari, 4), "v_measure": round(best_vm, 4)},
        "n_genera": n_genera,
        "n_languages": len(genus_of),
        "gled_validation": {
            "monophyletic": mono,
            "testable": testable,
            "rate": round(mono / testable, 4) if testable else None,
        },
        "branch_scale_heterogeneity": scale,
    }
    (gl.BASE / "data" / "params.json").write_text(json.dumps(params, indent=2))

    crown_vals = np.array(list(crown.values())) if crown else np.array([])

    lines = ["# Calibration report\n"]
    lines.append("## Primary criterion: Glottolog taxonomic-depth cut\n")
    lines.append(f"- **Selected depth d = {best_d}** levels below the family root.")
    lines.append(
        f"- Agreement with WALS/Dryer genera: **ARI = {best_ari:.3f}**, "
        f"V-measure = {best_vm:.3f} (on {len(wl)} WALS-labeled languages)."
    )
    lines.append(
        f"- Yields **{n_genera} genera** over {len(genus_of)} in-scope " f"languages.\n"
    )

    lines.append("### Depth sweep (agreement with WALS)\n")
    lines.append("| depth d | ARI | V-measure |")
    lines.append("|---------|-----|-----------|")
    for d, ari, vm in sweep:
        mark = "  **<- selected**" if d == best_d else ""
        lines.append(f"| {d} | {ari:.3f} | {vm:.3f} |{mark}")
    lines.append("")

    lines.append("## Indo-European yardstick\n")
    lines.append(
        f"At d = {best_d}, Indo-European resolves into {len(ie_rows)} genera "
        f"(WALS has 11):\n"
    )
    lines.append("| genus (Glottolog name) | size | WALS genera covered |")
    lines.append("|------------------------|------|---------------------|")
    for name, size, wals_genera in ie_rows:
        lines.append(f"| {name} | {size} | {', '.join(wals_genera)} |")
    lines.append("")

    lines.append("## GLED validation\n")
    if testable:
        lines.append(
            f"- Of {testable} genera testable in GLED trees (>=2 leaves), "
            f"**{mono} ({mono/testable:.1%}) are monophyletic**."
        )
    if crown_vals.size:
        lines.append(
            f"- Genus crown depth in GLED (relative divergence): median "
            f"{np.median(crown_vals):.3f}, IQR "
            f"[{np.percentile(crown_vals,25):.3f}, "
            f"{np.percentile(crown_vals,75):.3f}], max "
            f"{crown_vals.max():.3f}."
        )
        lines.append(
            "  The spread quantifies how *unequal* in depth the genera "
            "still are — the residual that phase-2 GLED/ASJP refinement "
            "targets (e.g. subdividing very deep genera)."
        )
    lines.append("")

    lines.append("## Why not a GLED branch-length threshold?\n")
    lines.append(
        "Each family's own ARI-optimal GLED cut height varies enormously, "
        "so no single distance threshold reproduces genera across "
        "families:\n"
    )
    lines.append(
        f"- Per-family optimal cut height (n={scale['families']} families): "
        f"median {scale['median']:.3f}, IQR [{scale['q25']:.3f}, "
        f"{scale['q75']:.3f}], range [{scale['min']:.3f}, "
        f"{scale['max']:.3f}]."
    )
    lines.append(
        f"- That ~{scale['max']/max(scale['min'],1e-6):.0f}x spread is why "
        "GLED divergences are used for validation and refinement, not as "
        "the primary cut.\n"
    )

    (gl.BASE / "data" / "calibration_report.md").write_text("\n".join(lines))
    print("Wrote data/params.json and data/calibration_report.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
