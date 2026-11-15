# Cavity redo and orthologue clean-up (2026-11-13)

Scripts: `scripts/21_cavity_robustness.py`, `scripts/22_orthologues.py`; shared code in `src/abhd/structops.py` (scripts 19 and 20 now
import it instead of loading each other by path). Tables: `data/derived/cavity_*.csv`, `orthologue_sets.csv`, `orth_*.csv`.
They answer the first two items of reports/11 and supersede the pocket-size and orthologue statements in reports/10.

## Cavity redo
Design. For every gene: model {CATH-trimmed domain, full-length AlphaFold model with atoms of pLDDT < 50 removed} x anchor {native:
nucleophile atom, else elbow CA; ca_trunc: elbow CA with the elbow residue's side chain removed for every gene, so Gly and Ser
genes are compared like for like} x 3 grid orientations (identity + 2 random rotations with sub-voxel shifts) x probe clearance
{3.0, 3.2, 3.6 A}; grid half-width 18 A (was 14). 3,960 variant runs. A contrast is called **consistent** only if the sign agrees in all
18 (model x orientation x clearance) combinations.

Results.
- Genome-wide the cavity ranking is stable within a model (Spearman 0.80-0.96 across clearance and anchor) but **only 0.61-0.70 between the
  trimmed-domain and full-length model**: the choice of model matters.
- 23% of like-for-like variants still touch the box (ABHD5 61%, ABHD6 44%, ABHD13 39%, ABHD4 33%), so their volumes are lower bounds.
- Median like-for-like volumes (3.2 A / 3.6 A; range across variants in `cavity_summary.csv`): ABHD6 436/141, ABHD5 375/184, ABHD3 309/51,
  ABHD2 283/34, MGLL 264/100, ABHD10 249/86, ABHD1 199/38, ABHD4 188/14, ABHD12 161/14, ABHD16A 155/64, ABHD13 151/70, ABHD15 119/9,
  ABHD17A 188/55, LYPLA1 39/11, ABHD14B 15/1.

| Contrast | Like-for-like (ca_trunc) | Native anchor |
|---|---|---|
| ABHD5 > ABHD4 | 18/18, median +189 | 18/18 |
| ABHD15 < ABHD3 | 18/18, median -116 | 17/18 |
| ABHD15 < ABHD2 | 18/18, median -51 | 14/18 |
| ABHD15 vs ABHD1 | 12/18 smaller (not consistent) | 12/18 |
| ABHD2 < ABHD6 | 18/18, median -119 | 14/18 |
| ABHD12 < MGLL | 18/18, median -81 | 18/18 |
| ABHD6 > MGLL | 16/18 | 17/18 |
| ABHD2 vs MGLL | 8 larger / 10 smaller (no consistent order) | 9/9 |
| ABHD10 vs MGLL, ABHD3 vs ABHD2 | not consistent | not consistent |
| ABHD16A < MGLL | 17/18 | 17/18 |

Interpretation.
- **Corrects reports/10.** ABHD15's cavity is smaller than ABHD3's and ABHD2's in every like-for-like variant (the earlier "normal-sized
  pocket" came from non-like-for-like anchors and one model). The ABHD5 cavity is larger than ABHD4's in all variants (and is
  truncated by the box, so the gap is, if anything, underestimated).
- Robust ordering in the panel: ABHD6 and ABHD5 have the largest cavities; ABHD2 is smaller than ABHD6; ABHD12 is smaller than MGLL;
  ABHD2 and MGLL cannot be ordered. No claim is made for the other pairs.
- These are volumes of apo models with a crude grid method; a size difference is a descriptor, not evidence about substrate or
  inhibitor selectivity.

## Orthologue clean-up
Design. Paginated UniProt retrieval (`gene:<NAME> AND taxonomy_id:7742`, all pages); each candidate is scored (BLOSUM62 local) against
all 110 human sequences and kept only if its best human hit is the target gene with at least a 10% margin over the runner-up, and
>=40% identity and >=70% coverage to the human domain; one sequence per organism (best score); human excluded; taxonomic class from the UniProt lineage.

| Target | Candidates | Kept | Not kept (other paralog, low score, duplicate organism) | Kept by class |
|---|---|---|---|---|
| ABHD15 | 399 | 363 | 31 | birds 235, mammals 94, ray-finned fish 16, amphibians 5, lepidosaurs 5, turtles 5, crocodilians 3 |
| ABHD5 | 720 | 434 | 277 | birds 282, mammals 95, fish 31, lepidosaurs 13, turtles 5, amphibians 5, crocodilians 3 |
| ABHD4 | 380 | 165 | 212 | mammals 99, fish 48, lepidosaurs 11, amphibians 4, crocodilians 2, turtles 1 |
| ABHD3, ABHD2, ABHD1 | 618, 606, 237 | 429, 240, 121 | | in `orthologue_sets.csv` |

(Half of the ABHD4 and ABHD5 candidates fail the best-hit rule: the gene-name query returns many entries of the other paralog.)

Results.
- **ABHD5 Asn153 is Asn in 434/434 retained sequences** (earlier "one Ser" was a mislabelled ABHD4-like entry). All classes.
- **ABHD15 elbow Gly is not universal**: Gly in birds (233/235), mammals (92/94), amphibians (5/5), lepidosaurs (5/5), turtles (4/5),
  teleost fish (12/16); **Ser (followed by Thr) in all 3 crocodilians and 3 early-diverging ray-finned fish (alligator gar, arapaima, bichir)**, and
  Ser in 2 mammals (hedgehog, desert hamster; the same residue follows, so these may be sequence errors). Birds are the crocodilians' sister group and are
  Gly, so the crocodilian Ser is either a reversion or an ancestral state retained; the data cannot choose without a tree. The statement
  "the nucleophile loss is vertebrate-wide" is replaced by: "Gly at the elbow in the large majority of sampled amniotes and amphibians; Ser in a
  few basal fish and crocodilians; polarity unresolved". Whether those Ser-bearing proteins have a complete triad was not tested.
- ABHD15 His391 (100%) and Asp360 (99.7%) remain at mid-rank percentiles 0.96 and 0.92 of domain positions (ties no longer penalised);
  the elbow Gly itself is at 0.68. Active paralogs: nucleophile Ser 100% in ABHD2/3/4, 98% in ABHD1.
- **ABHD4 180-198 segment** (the one segment displaced between ABHD4 and ABHD5): in ABHD4 orthologues it is conserved (median 0.76; P187, P188, W190 and K192 at 0.93-0.99),
  and 7 of its 19 positions have **no aligned residue in ABHD5**. The ABHD4-ABHD5 difference is therefore a conserved ABHD4-specific insertion (PPxW...K), not a
  remodelled ABHD5 surface; in the ABHD5 orthologues the aligned positions are conserved at median 0.61. This replaces the "ABHD5 pocket
  remodelled" and "lipid-droplet surface" speculation in reports/10 and also gives a candidate reason for the smaller ABHD4 cavity (a
  hypothesis: the insertion sits against the pocket, 6 shell residues).

## Remaining caveats
- Orthologue sets are 55-65% birds (taxon sampling), 99% unreviewed; paralog rejection uses only the 110 human sequences.
- No phylogeny was built, so ancestral states and number of losses are not inferred.
- Cavity volumes remain model- and method-dependent; ABHD6, ABHD5, ABHD13 and ABHD4 are box-truncated.
