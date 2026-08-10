#!/usr/bin/env python
"""
Assign every in-scope Glottolog language to a genus, using the calibrated
Glottolog taxonomic-depth cut (data/params.json).

A genus is identified by the Glottolog node that defines it: ``genus_id`` is that
node's glottocode (stable) and ``genus_name`` its Glottolog name. WALS genus
labels are attached for reference where available.

Writes:
  data/genera.tsv           one row per language
  data/genera_summary.tsv   one row per genus
"""

from __future__ import annotations

import csv
import json
from collections import Counter, defaultdict

import genlib as gl


def main() -> int:
    params = json.loads((gl.BASE / "data" / "params.json").read_text())
    depth = params["depth"]
    print(f"Assigning genera at depth d = {depth} ...")

    glottolog = gl.Glottolog.load()
    wals = gl.load_wals()
    trees = gl.load_gled_trees(glottolog)

    genus_of = gl.assign_genera(glottolog, depth)
    support = gl.gled_genus_support(genus_of, trees)

    # Per-genus members and majority WALS label.
    members = defaultdict(list)
    for gc, genus in genus_of.items():
        members[genus].append(gc)
    wals_label = {}
    for genus, gcs in members.items():
        labels = [wals[g] for g in gcs if g in wals]
        wals_label[genus] = Counter(labels).most_common(1)[0][0] if labels else ""

    def method_of(gc: str, genus: str) -> str:
        return "isolate" if glottolog.is_isolate(gc) and members[genus] == [gc] else "glottolog_depth"

    def confidence_of(gc: str, genus: str) -> str:
        if method_of(gc, genus) == "isolate":
            return "isolate"
        if wals_label[genus]:
            return "wals-anchored"
        if support.get(genus, {}).get("monophyletic") and support[genus]["n_leaves"] >= 2:
            return "gled-supported"
        return "glottolog-only"

    # ----- per-language table -----
    lang_path = gl.BASE / "data" / "genera.tsv"
    with lang_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "glottocode", "name", "family", "family_gc",
            "genus_id", "genus_name", "method", "wals_genus", "confidence",
        ])
        for r in sorted(glottolog.in_scope_languages(), key=lambda r: r["glottocode"]):
            gc = r["glottocode"]
            genus = genus_of[gc]
            fam_gc = glottolog.family_gc(gc)
            w.writerow([
                gc, r["name"], glottolog.family_name(fam_gc), fam_gc,
                genus, glottolog.name(genus), method_of(gc, genus),
                wals_label[genus], confidence_of(gc, genus),
            ])
    print(f"Wrote {lang_path} ({len(genus_of)} languages)")

    # ----- per-genus summary -----
    sum_path = gl.BASE / "data" / "genera_summary.tsv"
    with sum_path.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow([
            "genus_id", "genus_name", "family", "family_gc", "n_languages",
            "method", "wals_genus", "gled_monophyletic", "gled_crown_height",
            "example_languages",
        ])
        for genus, gcs in sorted(members.items(), key=lambda kv: -len(kv[1])):
            rep = gcs[0]
            fam_gc = glottolog.family_gc(rep)
            sup = support.get(genus, {})
            mono = "" if sup.get("n_leaves", 0) < 2 else str(sup.get("monophyletic"))
            crown = round(sup["crown_height"], 4) if sup.get("crown_height") else ""
            examples = ", ".join(glottolog.name(g) for g in gcs[:4])
            is_iso = glottolog.is_isolate(rep) and gcs == [rep]
            w.writerow([
                genus, glottolog.name(genus), glottolog.family_name(fam_gc), fam_gc,
                len(gcs), "isolate" if is_iso else "glottolog_depth",
                wals_label[genus], mono, crown, examples,
            ])
    print(f"Wrote {sum_path} ({len(members)} genera)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
