# Stage 6: global comparison, FunFam concordance, and the ABHD name

Scripts: `scripts/12_global_comparison.py`, `scripts/13_abhd_label_analysis.py`. Outputs: `data/derived/dist_*.csv`,
`global_label_tests.csv`, `global_clusters.csv`, `global_concordance.csv`, `discordant_pairs.csv`,
`abhd_label_summary.csv`, `abhd_label_by_gene.csv`, `funfam_purity.csv`; figures `stage6_label_purity.html`,
`stage6b_abhd_label.html`.

## Spaces
TM-align structure distance (tmtools, on trimmed domains), BLOSUM62 Smith-Waterman local-alignment sequence distance (the first run used a global score clipped at zero, which tied 94% of pairs at 1.0; replaced after review), TM-align normalised by the shorter chain (`tmshort`), ProtBert embedding,
FunFam HMM-HMM E-values (103 genes), SOAP, triad RMSD. Foldseek was not used: a download of its binary was blocked, so
TM-align (via the `tmtools` Python package) provides the structural comparison.

## Do ABHD genes cluster? (central claim)
| Space | ABHD genes' nearest neighbour is ABHD (95% CI) | ... skipping same-stem paralogs | ... skipping the closest 3.4% of pairs (TM < 0.2 equivalent) | R-precision | Chance |
|---|---|---|---|---|---|
| TM-align | 80% (60-95) | 75% | 45% | 0.40 | 0.17 |
| TM-align, shorter-chain norm. | 85% (70-100) | 85% | 60% | 0.35 | 0.17 |
| Sequence (local) | 80% (60-95) | 70% | 60% | 0.48 | 0.17 |
| FunFam (17 genes) | 82% (65-100) | 71% | 71% | 0.42 | 0.16 |
| ProtBert | 50% | 25% | | 0.26 | 0.17 |
| SOAP | 60% | 55% | | 0.24 | 0.17 |
(R-precision = fraction ABHD among the R = 19 nearest neighbours. Permutation test of "ABHD genes closer to each other than to
others": p = 0.0002 (TM-align), 0.0006 (sequence), 0.009 (FunFam), 0.01 (ProtBert), 0.19 (SOAP); permutations treat genes as
independent. A length-stratified permutation, shuffling labels within domain-length quartiles, raises the chance R-precision
from 0.17 to 0.27 for TM-align but the observed 0.40 stays outside the null 95% range (0.19-0.36, p=0.004).
`scripts/17_review_controls.py`; `review_label_robustness.csv`.)

- The strong form of the claim - that the name carries no structural or sequence information - is **rejected** (with the caveat that much of the structure signal is carried by close paralog sub-groups: ABHD1/2/3/15, ABHD4/5, ABHD17A/B/C, ABHD14A/B, ABHD16A/B): ABHD genes
  usually sit next to other ABHD genes. This holds with same-stem paralogs (ABHD17A/B/C, ABHD14A/B, ABHD16A/B) skipped and
  with neighbours restricted to domain length within 20% (ABHD median 244 residues vs 310 for others), but dropping the closest 3.4% of pairs (the share TM < 0.2 removes) lowers the nearest-neighbour rate in every space (TM 0.45, sequence 0.60, FunFam 0.71): the signal is paralog-driven in sequence as well as structure. The stem-based paralog filter misses ABHD1/2/3/15 and ABHD4/5 pairs.
- The weaker form (the set is not compact) is **supported by R-precision, but not uniformly by clustering**. 60% of the 19
  nearest neighbours of an ABHD gene are not ABHD in structure space (52% in local-sequence space). In a 20-cluster average-
  linkage TM-align cut ABHD genes occupy 4 clusters: 12 ABHD genes share one 30-member cluster with 18 non-ABHD proteins (MGLL,
  LYPLA1/2, LYPLAL1, APEH, ESD, ACOT1/2/4, BAAT, CMBL, AFMID, the DPP family and FAP); five (ABHD4/5/6/8/11) share a
  25-member cluster with lipases, epoxide hydrolases and NDRGs. **Occupying 4 of 20 clusters is evidence for clustering**:
  random 20-gene sets occupy 8.8 (p=0.001; the number is linkage-dependent: 2 single, 4 average, 6 complete). Only 3 ABHD
  genes sit in majority-ABHD TM-align clusters. With the corrected local-sequence distance a 20-cluster cut puts 19 of 20
  ABHD genes in majority-ABHD clusters (3 clusters), so by 20-cluster sequence clustering the name looks more cohesive than by structure (but see the quantile-matched exclusion above); this
  is the opposite of the first run, which used a degenerate seqid. ABHD6, ABHD8, ABHD10 and ABHD11 have a non-ABHD nearest
  neighbour in structure space (BPHL, MEST, MGLL, BPHL), and EPHX3 is among the 3 nearest neighbours of ABHD4/5/6/8/11.
- Label definition: the strict flag is "gene symbol starts with ABHD" (20 genes). UniProt calls EPHX3 and EPHX4 "Abhydrolase
  domain-containing protein 9/7" and LIPJ/K/M/N "lipase-like abhydrolase domain-containing 1-4". Under this broad definition
  (`abhd_family`, 26 genes) the TM nearest-neighbour rate is 85% and R-precision 0.43, so the conclusion holds. ABHD12B
  (Q7Z5M8, reviewed) exists in human but is in neither source, so it is not in the universe.
- FunFams: all 15 FunFams that contain an ABHD gene contain only ABHD genes, but 14 of them have a single human member, so
  this is nearly vacuous; ABHD17A/B/C is the only real test (they share one FunFam); 86 of 103 FunFam members are not ABHD-named.

## Concordance between spaces (Spearman of pairwise distances)
TM-align vs local sequence 0.52, vs FunFam 0.38; sequence vs FunFam 0.56; triad RMSD vs TM-align 0.48; ProtBert is nearly
uncorrelated with the others (0.02-0.10) and SOAP is weak. Agreement with FunFam at 20 clusters (adjusted Rand): sequence
0.23, TM-align 0.12, ProtBert 0.06, SOAP 0.02: FunFams are far finer than a 20-cluster cut.

## Discordant pairs
With the corrected sequence distance only 12 pairs are discordant at the 5%/50% rank thresholds. One is structure-close but
sequence-far (MEST-NDRG1, TM distance 0.26); 11 are sequence-close but structure-far, dominated by SERHL (poorly modelled)
and ACOT6 (a pseudo-enzyme-like outlier), which look like model or domain-trimming artefacts. The earlier list of 49 pairs
(NDRGs vs SPG21/BPHL/MEST) was largely a product of the tied sequence distances.

## Caveats
- Human-only, 110 genes: "closest neighbour" is relative to this set, not the full superfamily.
- TM-align and sequence distances are computed on CATH domain spans of different lengths (PREP/PREPL include a beta-propeller;
  NLGN/CES/ACHE domains are 520-610 residues). Length matching, a shorter-chain TM normalisation and a length-stratified
  permutation were added as controls. Re-trimming PREP/PREPL to the hydrolase fold has not been done.
- ProtBert embeddings are mean-pooled and compress distances; ranks are meaningful, values are not.
