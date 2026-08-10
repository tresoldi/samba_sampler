"""
Shared library for the genera pipeline: data loading, tree cutting, and
Glottolog-structure extension. Imported by the numbered pipeline scripts.

Design notes
------------
* GLED branch lengths encode *relative* lexical divergence, so a genus is a
  maximal clade whose crown height (max branch-length distance from a node down
  to its leaves) does not exceed a calibrated threshold ``T``.
* Not every language of a GLED-covered family is a tree leaf, and some families
  have no tree at all. Those languages are placed by *Glottolog-structure
  extension*: a language joins the genus of its closest classified relative
  (deepest shared Glottolog ancestry).
"""

from __future__ import annotations

import csv
import re
import sys
from collections import defaultdict
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

# Make the bundled samba_sampler importable without installing it.
_REPO = Path(__file__).resolve().parent.parent.parent
sys.path.insert(0, str(_REPO / "src"))
from samba_sampler import newick  # noqa: E402

BASE = Path(__file__).resolve().parent.parent
ETC = _REPO / "src" / "samba_sampler" / "etc"
GLED_TREES = _REPO / "scaffolding" / "gled.trees"
SOURCES = BASE / "data" / "sources"

# Glottolog categories that are not genealogical language units.
EXCLUDED_CATEGORIES = {
    "Bookkeeping",
    "Artificial Language",
    "Sign Language",
    "Unattested",
    "Unclassifiable",
    "Speech Register",
    "Mixed Language",
    "Pidgin",
}

_GLOTTOCODE_RE = re.compile(r"([a-z0-9]{4}[0-9]{4})$")


def leaf_glottocode(name: Optional[str]) -> Optional[str]:
    """Extract a trailing glottocode from a GLED leaf label (``Name_gccc1234``)."""
    if not name:
        return None
    m = _GLOTTOCODE_RE.search(name)
    return m.group(1) if m else None


# --------------------------------------------------------------------------- #
# Glottolog
# --------------------------------------------------------------------------- #
class Glottolog:
    def __init__(self, rows: List[dict]):
        self.rows = rows
        self.by_gc: Dict[str, dict] = {r["glottocode"]: r for r in rows}

    @classmethod
    def load(cls) -> "Glottolog":
        path = sorted(ETC.glob("glottolog.*.tsv"))[-1]
        with path.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        return cls(rows)

    def ancestors(self, gc: str) -> List[str]:
        """Glottocodes from family root down to the immediate parent (top-first)."""
        r = self.by_gc.get(gc)
        if not r or not r["ancestors"]:
            return []
        return r["ancestors"].split(";")

    def family_gc(self, gc: str) -> str:
        """Top-level family glottocode; isolates/top-level nodes return themselves."""
        anc = self.ancestors(gc)
        return anc[0] if anc else gc

    def family_name(self, gc: str) -> str:
        r = self.by_gc.get(gc, {})
        return r.get("family") or r.get("name") or ""

    def name(self, gc: str) -> str:
        return self.by_gc.get(gc, {}).get("name", gc)

    def is_isolate(self, gc: str) -> bool:
        r = self.by_gc.get(gc, {})
        return r.get("isolate") == "True" or not r.get("ancestors")

    def in_scope_languages(self) -> List[dict]:
        return [
            r
            for r in self.rows
            if r["level"] == "language" and r["category"] not in EXCLUDED_CATEGORIES
        ]

    def lineage(self, gc: str) -> List[str]:
        """Full top-down path including the node itself (family..node)."""
        return self.ancestors(gc) + [gc]

    def genus_node(self, gc: str, depth: int) -> str:
        """
        The Glottolog node defining ``gc``'s genus under a taxonomic-depth cut.

        Cuts the classification ``depth`` levels below the family root. Shallow
        families are capped at their deepest *group* (never the language itself),
        so they are not shattered into singletons; isolates and top-level
        families are their own genus.
        """
        lin = self.lineage(gc)
        if len(lin) == 1:
            return gc  # isolate / top-level family
        idx = min(depth, len(lin) - 1)
        if idx == len(lin) - 1:
            idx -= 1  # cap at the parent group, never the language node
        return lin[idx]

    def shared_depth(self, gc1: str, gc2: str) -> int:
        """Length of the shared top-down ancestor prefix of two languoids."""
        a, b = self.lineage(gc1), self.lineage(gc2)
        n = 0
        for x, y in zip(a, b):
            if x != y:
                break
            n += 1
        return n


# --------------------------------------------------------------------------- #
# WALS
# --------------------------------------------------------------------------- #
def load_wals() -> Dict[str, str]:
    """Return {glottocode: WALS genus} for WALS rows carrying both."""
    path = SOURCES / "wals_languages.csv"
    out: Dict[str, str] = {}
    with path.open(encoding="utf-8") as fh:
        for w in csv.DictReader(fh):
            gc, genus = w.get("Glottocode"), w.get("Genus")
            if gc and genus:
                out[gc] = genus
    return out


# --------------------------------------------------------------------------- #
# GLED trees
# --------------------------------------------------------------------------- #
class GledTree:
    def __init__(self, path: Path, glottolog: Glottolog):
        self.path = path
        self.root = newick.read(str(path))[0]
        self.leaves = [
            (lf, gc)
            for lf in self.root.get_leaves()
            if (gc := leaf_glottocode(lf.name))
        ]
        self.leaf_gcs = [gc for _, gc in self.leaves]
        # Associate the tree with a Glottolog family via majority vote of leaves.
        fams = [
            glottolog.family_gc(gc) for gc in self.leaf_gcs if gc in glottolog.by_gc
        ]
        self.family_gc = max(set(fams), key=fams.count) if fams else None

    def node_heights(self) -> Dict[int, float]:
        h: Dict[int, float] = {}
        for n in self.root.walk(mode="postorder"):
            if n.is_leaf:
                h[id(n)] = 0.0
            else:
                h[id(n)] = max(h[id(c)] + (c.length or 0.0) for c in n.descendants)
        return h

    def internal_heights(self) -> List[float]:
        h = self.node_heights()
        return [h[id(n)] for n in self.root.walk() if not n.is_leaf]

    def cut(self, threshold: float) -> Dict[str, int]:
        """
        Partition leaves into genera: maximal clades of crown height <= threshold.
        Returns {glottocode: local genus index}.
        """
        h = self.node_heights()
        genus_of: Dict[str, int] = {}
        gid = 0
        stack = [self.root]
        while stack:
            node = stack.pop()
            if h[id(node)] <= threshold:
                for lf in node.get_leaves():
                    gc = leaf_glottocode(lf.name)
                    if gc:
                        genus_of[gc] = gid
                gid += 1
            else:
                stack.extend(node.descendants)
        return genus_of


def load_gled_trees(glottolog: Glottolog) -> List[GledTree]:
    return [GledTree(p, glottolog) for p in sorted(GLED_TREES.glob("*.tree"))]


# --------------------------------------------------------------------------- #
# Glottolog-structure extension
# --------------------------------------------------------------------------- #
def extend_by_glottolog(
    targets: List[str],
    anchors: Dict[str, str],
    glottolog: Glottolog,
) -> Dict[str, str]:
    """
    Assign each target glottocode to a genus by nearest classified relative.

    A target joins the genus of the anchor with which it shares the deepest
    Glottolog ancestry. Ties break toward the genus most common among the
    deepest-tied anchors (a proxy for the larger, better-attested genus).

    @param targets: glottocodes needing a genus (may include anchors; those keep
        their own label).
    @param anchors: {glottocode: genus} already-assigned reference languages.
    @param glottolog: a loaded Glottolog.
    @return: {glottocode: genus} covering every target.
    """
    anchor_items = list(anchors.items())
    out: Dict[str, str] = {}
    for gc in targets:
        if gc in anchors:
            out[gc] = anchors[gc]
            continue
        best_depth = -1
        tied: List[str] = []
        for a_gc, genus in anchor_items:
            d = glottolog.shared_depth(gc, a_gc)
            if d > best_depth:
                best_depth, tied = d, [genus]
            elif d == best_depth:
                tied.append(genus)
        if tied and best_depth > 0:
            out[gc] = max(set(tied), key=tied.count)
        # best_depth == 0 (no shared family) should not happen within a family;
        # such a target is left unassigned for the caller to handle.
    return out


def assign_genera(glottolog: "Glottolog", depth: int) -> Dict[str, str]:
    """Assign every in-scope language to a genus node at the given cut depth."""
    return {
        r["glottocode"]: glottolog.genus_node(r["glottocode"], depth)
        for r in glottolog.in_scope_languages()
    }


def gled_monophyly(
    genus_of: Dict[str, str], trees: List["GledTree"]
) -> Tuple[int, int, Dict[str, float]]:
    """
    Check, per GLED tree, whether assigned genera form monophyletic groups.

    A genus is monophyletic in a tree if the smallest clade containing all its
    leaves present in that tree contains no leaves of another genus.

    @return: (monophyletic_count, testable_count, {genus: max_crown_height}) where
        crown height is the genus subtree height in GLED (relative divergence).
    """
    mono, testable = 0, 0
    crown: Dict[str, float] = {}
    for tree in trees:
        heights = tree.node_heights()
        # Leaves present in this tree, grouped by assigned genus.
        present = [(lf, gc) for lf, gc in tree.leaves if gc in genus_of]
        genera_here = {}
        for lf, gc in present:
            genera_here.setdefault(genus_of[gc], []).append(lf)
        for genus, leaves in genera_here.items():
            if len(leaves) < 2:
                continue
            testable += 1
            mrca = _mrca(leaves)
            clade_leaves = set(id(x) for x in mrca.get_leaves())
            assigned = set(id(lf) for lf, _ in present)
            # Monophyletic iff every leaf under the MRCA that is itself assigned
            # belongs to this same genus.
            others = [
                gc
                for lf, gc in present
                if id(lf) in clade_leaves and genus_of[gc] != genus
            ]
            if not others:
                mono += 1
            crown[genus] = max(crown.get(genus, 0.0), heights[id(mrca)])
    return mono, testable, crown


def gled_genus_support(
    genus_of: Dict[str, str], trees: List["GledTree"]
) -> Dict[str, dict]:
    """Per-genus GLED evidence: {genus: {monophyletic, n_leaves, crown_height}}."""
    status: Dict[str, dict] = {}
    for tree in trees:
        heights = tree.node_heights()
        present = [(lf, gc) for lf, gc in tree.leaves if gc in genus_of]
        by_genus: Dict[str, list] = {}
        for lf, gc in present:
            by_genus.setdefault(genus_of[gc], []).append(lf)
        present_ids = {id(lf) for lf, _ in present}
        for genus, leaves in by_genus.items():
            prev = status.get(
                genus, {"monophyletic": True, "n_leaves": 0, "crown_height": 0.0}
            )
            prev["n_leaves"] += len(leaves)
            if len(leaves) >= 2:
                mrca = _mrca(leaves)
                clade = {id(x) for x in mrca.get_leaves()}
                intruders = any(
                    id(lf) in clade and genus_of[gc] != genus for lf, gc in present
                )
                if intruders:
                    prev["monophyletic"] = False
                prev["crown_height"] = max(prev["crown_height"], heights[id(mrca)])
            status[genus] = prev
    return status


def _mrca(leaves: list):
    """Most recent common ancestor of a list of newick leaf Nodes."""
    anc_sets = [set(id(a) for a in ([lf] + lf.ancestors)) for lf in leaves]
    common = set.intersection(*anc_sets)
    # Choose the deepest common ancestor: the one that is not an ancestor of
    # another common ancestor. Walk up from the first leaf until we hit a common id.
    for node in [leaves[0]] + leaves[0].ancestors:
        if id(node) in common:
            return node
    return leaves[0].ancestors[-1]


# --------------------------------------------------------------------------- #
# ASJP wordlists (phase 2)
# --------------------------------------------------------------------------- #
# Minimum number of shared concepts required to compute a language-pair distance.
ASJP_MIN_SHARED = 10

# ASJP modifier / annotation symbols stripped before comparing word forms. The
# ASJP base alphabet is letters plus the digits 3 5 7 8 (which encode phonemes),
# so digits are kept; only combiners and loan/uncertainty markers are removed.
_ASJP_STRIP = str.maketrans("", "", '~$*"%. ')


def _clean_asjp(form: str) -> str:
    """Normalize an ASJP word form: drop modifiers, keep base phoneme symbols."""
    return form.strip().translate(_ASJP_STRIP)


def _ldn(a: str, b: str) -> float:
    """Levenshtein distance between two strings, normalized by the longer length."""
    if a == b:
        return 0.0
    la, lb = len(a), len(b)
    if la == 0 or lb == 0:
        return 1.0
    prev = list(range(lb + 1))
    for i, ca in enumerate(a, 1):
        cur = [i]
        for j, cb in enumerate(b, 1):
            cost = 0 if ca == cb else 1
            cur.append(min(prev[j] + 1, cur[j - 1] + 1, prev[j - 1] + cost))
        prev = cur
    return prev[lb] / max(la, lb)


class Asjp:
    """
    ASJP wordlists keyed by Glottocode.

    Each glottocode maps to ``{concept_id: [asjp_form, ...]}``, aggregated over
    every ASJP doculect that Glottolog assigns to that glottocode.
    """

    def __init__(self, words: Dict[str, Dict[str, List[str]]]):
        self.words = words
        self._cache: Dict[Tuple[str, str], Optional[float]] = {}

    @classmethod
    def load(cls, keep: Optional[set] = None) -> "Asjp":
        """
        Load ASJP forms from the cached CLDF dump.

        @param keep: optional set of glottocodes to restrict loading to (keeps the
            in-memory footprint small when only a subset of languages is needed).
        """
        gc_of: Dict[str, str] = {}
        with (SOURCES / "asjp_languages.csv").open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                gc = r.get("Glottocode")
                if gc:
                    gc_of[r["ID"]] = gc

        words: Dict[str, Dict[str, List[str]]] = {}
        with (SOURCES / "asjp_forms.csv").open(encoding="utf-8") as fh:
            for r in csv.DictReader(fh):
                gc = gc_of.get(r["Language_ID"])
                if not gc or (keep is not None and gc not in keep):
                    continue
                form = _clean_asjp(r.get("Form") or "")
                if form:
                    words.setdefault(gc, {}).setdefault(r["Parameter_ID"], []).append(
                        form
                    )
        return cls(words)

    def has(self, gc: str) -> bool:
        return gc in self.words

    def covered(self, gcs: Iterable[str]) -> List[str]:
        return [gc for gc in gcs if gc in self.words]

    def distance(self, gc1: str, gc2: str) -> Optional[float]:
        """
        ASJP lexical distance: mean normalized Levenshtein distance over shared
        concepts. For a concept with several forms on either side the minimum
        pairwise distance is taken. Returns ``None`` when fewer than
        ``ASJP_MIN_SHARED`` concepts are shared.
        """
        key = (gc1, gc2) if gc1 <= gc2 else (gc2, gc1)
        if key in self._cache:
            return self._cache[key]
        w1, w2 = self.words.get(gc1), self.words.get(gc2)
        result: Optional[float] = None
        if w1 and w2:
            shared = set(w1) & set(w2)
            if len(shared) >= ASJP_MIN_SHARED:
                total = 0.0
                for c in shared:
                    total += min(_ldn(f1, f2) for f1 in w1[c] for f2 in w2[c])
                result = total / len(shared)
        self._cache[key] = result
        return result


# Caps that bound the cost of the separation test on very large families while
# keeping it meaningful (a genus is compared against its nearest relatives).
_ASJP_MAX_GENUS = 60  # genus members entering the intra-genus comparison
_ASJP_MAX_OTHERS = 150  # nearest same-family non-members entering the inter comparison


def asjp_genus_support(
    genus_of: Dict[str, str],
    glottolog: "Glottolog",
    asjp: "Asjp",
    target_genera: Iterable[str],
) -> Dict[str, dict]:
    """
    Test whether target genera are lexically cohesive and separated in ASJP.

    A genus is ``supported`` when its mean intra-genus ASJP distance is smaller
    than the mean distance from its members to their nearest same-family
    relatives (a positive silhouette-style separation): ASJP wordlists
    independently recover it as a tighter cluster than its surroundings. This is
    a *relative* criterion, so — unlike a global distance threshold — it does not
    assume ASJP divergences are on a common scale across families.

    Comparison against "the rest of the family" is capped at the nearest
    relatives (by shared Glottolog depth) to bound cost and to make the test
    stricter: a genus must separate from the languages most likely to be confused
    with it, not merely from distant family members.

    @return: {genus_id: {n_asjp, intra, inter, score, supported, reason}} for
        every target genus (supported=False with a reason where untestable).
    """
    target_set = set(target_genera)

    fam_members: Dict[str, List[Tuple[str, str]]] = defaultdict(list)
    for gc, genus in genus_of.items():
        fam_members[glottolog.family_gc(gc)].append((gc, genus))

    results: Dict[str, dict] = {}
    for members in fam_members.values():
        present_targets = {g for _, g in members if g in target_set}
        if not present_targets:
            continue
        covered = [(gc, g) for gc, g in members if asjp.has(gc)]

        for genus in present_targets:
            g_gcs = [gc for gc, g in covered if g == genus][:_ASJP_MAX_GENUS]
            others = [gc for gc, g in covered if g != genus]
            if len(others) > _ASJP_MAX_OTHERS:
                others.sort(
                    key=lambda o: max(
                        (glottolog.shared_depth(o, m) for m in g_gcs), default=0
                    ),
                    reverse=True,
                )
                others = others[:_ASJP_MAX_OTHERS]

            intra = [
                d
                for a, b in combinations(g_gcs, 2)
                if (d := asjp.distance(a, b)) is not None
            ]
            if not intra:
                results[genus] = {
                    "n_asjp": len(g_gcs),
                    "supported": False,
                    "reason": "fewer than 2 comparable ASJP wordlists in genus",
                }
                continue
            inter = [
                d
                for a in g_gcs
                for b in others
                if (d := asjp.distance(a, b)) is not None
            ]
            if not inter:
                results[genus] = {
                    "n_asjp": len(g_gcs),
                    "supported": False,
                    "reason": "no comparable ASJP wordlists outside genus in family",
                }
                continue
            a = sum(intra) / len(intra)
            b = sum(inter) / len(inter)
            results[genus] = {
                "n_asjp": len(g_gcs),
                "intra": round(a, 4),
                "inter": round(b, 4),
                "score": round((b - a) / max(a, b, 1e-9), 4),
                "supported": a < b,
                "reason": "",
            }
    return results


def genus_glottolog_depths(
    members_by_genus: Dict[str, List[str]], glottolog: Glottolog
) -> List[int]:
    """
    For each genus, the Glottolog depth (levels below family root) of the MRCA of
    its members. Used to calibrate the Glottolog-depth fallback rule.
    """
    depths: List[int] = []
    for members in members_by_genus.values():
        if len(members) < 2:
            continue
        lineages = [glottolog.lineage(m) for m in members]
        depth = 0
        for level in zip(*lineages):
            if len(set(level)) == 1:
                depth += 1
            else:
                break
        depths.append(max(0, depth - 1))  # 0-based levels below the family root
    return depths
