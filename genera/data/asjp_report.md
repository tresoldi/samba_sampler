# ASJP refinement (phase 2)

## Upgrading `glottolog-only` genera

- Target `glottolog-only` genera: **394**.
- Testable with ASJP (>=2 comparable wordlists inside the genus and >=1 outside it): **70**.
- **Upgraded to `asjp-supported`: 63** — mean intra-genus ASJP distance below the distance to nearest same-family relatives.
- Not upgraded: 7 tested but not separated; 324 untestable (insufficient ASJP coverage).

### Best-separated upgraded genera

| genus | family | n(ASJP) | intra | inter | score |
|-------|--------|---------|-------|-------|-------|
| Chug-Lish | Sino-Tibetan | 2 | 0.1255 | 0.7849 | 0.8401 |
| East Seram | Austronesian | 2 | 0.2 | 0.7064 | 0.7169 |
| Kolami-Naiki | Dravidian | 2 | 0.2069 | 0.6506 | 0.6819 |
| One | Nuclear Torricelli | 2 | 0.2609 | 0.7622 | 0.6578 |
| Albanian-Tosk | Indo-European | 3 | 0.3092 | 0.8362 | 0.6302 |
| Bullom | Atlantic-Congo | 3 | 0.3153 | 0.8374 | 0.6234 |
| Benyadu-Bekati | Austronesian | 3 | 0.2865 | 0.7329 | 0.6091 |
| Subanen | Austronesian | 3 | 0.263 | 0.6321 | 0.5839 |
| Sabakor | Nuclear Trans New Guinea | 2 | 0.3555 | 0.8485 | 0.5811 |
| Sambalic | Austronesian | 8 | 0.3126 | 0.7142 | 0.5622 |
| Moken-Moklen | Austronesian | 2 | 0.3596 | 0.7825 | 0.5404 |
| Limba | Atlantic-Congo | 2 | 0.4071 | 0.8599 | 0.5266 |

## Over-large genus subdivision (advisory)

- Genera with >= 100 languages: **6**; ASJP proposes a subdivision (positive silhouette) for **6** of them.
- Proposals are advisory: `genus_id`s are unchanged in this release (subdividing would coin new genus nodes, left to future work). Full detail in `data/asjp_subdivisions.tsv`.

| genus | family | languages | proposed k | silhouette | cluster sizes |
|-------|--------|-----------|-----------|-----------|---------------|
| Bantoid | Atlantic-Congo | 730 | 2 | 0.1394 | 249;1 |
| Oceanic | Austronesian | 521 | 2 | 0.1242 | 249;1 |
| Indo-Aryan | Indo-European | 220 | 2 | 0.4157 | 142;8 |
| Cameroun-Ubangian | Atlantic-Congo | 118 | 2 | 0.1975 | 107;1 |
| Northwest Germanic | Indo-European | 104 | 2 | 0.3018 | 42;25 |
| Loloish | Sino-Tibetan | 101 | 2 | 0.1506 | 27;1 |
