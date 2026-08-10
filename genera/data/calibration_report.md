# Calibration report

## Primary criterion: Glottolog taxonomic-depth cut

- **Selected depth d = 3** levels below the family root.
- Agreement with WALS/Dryer genera: **ARI = 0.867**, V-measure = 0.925 (on 2261 WALS-labeled languages).
- Yields **1551 genera** over 7650 in-scope languages.

### Depth sweep (agreement with WALS)

| depth d | ARI | V-measure |
|---------|-----|-----------|
| 1 | 0.447 | 0.899 |
| 2 | 0.773 | 0.931 |
| 3 | 0.867 | 0.925 |  **<- selected**
| 4 | 0.650 | 0.896 |
| 5 | 0.507 | 0.879 |
| 6 | 0.273 | 0.864 |

## Indo-European yardstick

At d = 3, Indo-European resolves into 10 genera (WALS has 11):

| genus (Glottolog name) | size | WALS genera covered |
|------------------------|------|---------------------|
| Indo-Aryan | 45 | Creoles and Pidgins, Indic |
| Northwest Germanic | 32 | Creoles and Pidgins, Germanic |
| Latino-Faliscan | 29 | Creoles and Pidgins, Romance |
| Iranian | 24 | Iranian |
| Slavic | 16 | Slavic |
| Nuclear Celtic | 6 | Celtic |
| Nuristani | 3 | Nuristani |
| Eastern Baltic | 2 | Baltic |
| Eastern-Western Armenian | 1 | Armenian |
| Greek | 1 | Greek |

## GLED validation

- Of 512 genera testable in GLED trees (>=2 leaves), **489 (95.5%) are monophyletic**.
- Genus crown depth in GLED (relative divergence): median 0.359, IQR [0.149, 1.073], max 70.227.
  The spread quantifies how *unequal* in depth the genera still are — the residual that phase-2 GLED/ASJP refinement targets (e.g. subdividing very deep genera).

## Why not a GLED branch-length threshold?

Each family's own ARI-optimal GLED cut height varies enormously, so no single distance threshold reproduces genera across families:

- Per-family optimal cut height (n=36 families): median 1.340, IQR [0.619, 2.976], range [0.227, 8.799].
- That ~39x spread is why GLED divergences are used for validation and refinement, not as the primary cut.
