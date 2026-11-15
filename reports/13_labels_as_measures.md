# Stage 13: CATH, FunFam, EC and UniProt annotations as measurements

Plan: `docs/redesign_plan.md` (v1 commit 490da0e; v2 after independent review, commit 8a4b316). Scripts 24-29. All numbers are from
`data/derived/` tables; acceptance checks against earlier outputs passed (listed per analysis). Date 2026-11-15.

## 1. Label provenance (script 24)
- 76 of 110 genes have a catalytic-activity record, 68 a Rhea reaction, 104 GO molecular-function terms, 49 a PDB entry.
- 80 genes have a complete EC3 label (identical to the 80 used in function transfer). Of these, 58 (73%) have at least one experimentally supported
  reaction (ECO:0000269), 51 have only experimental reactions for their EC3, 7 mixed, 14 non-experimental (curator or similarity), 8 no record.
  Among the 15 ABHD genes with a complete EC3, 10 have experimental support.
- EC3 is a chemistry class: 41 of 80 labelled genes are exactly 3.1.1. GO evidence codes are not in the download, so GO is not used for transfer.

## 2. Function transfer as a decay curve (script 25)
Pairwise identity (global alignment, trimmed domains): median 0.25, 95th percentile 0.33; only 10% of pairs exceed 30% and 1.7% exceed 40%,
so cutoffs above 40% do nothing and the informative range is about 25-40%.
- Checks passed: with no identity cutoff the accuracies reproduce script 14 exactly (for example structure 0.838 / family-excluded 0.675).
- **Below 30% identity no space beats the majority-class baseline**: accuracy 0.39 (structure), 0.55 (sequence), 0.51 (FunFam), 0.575 (elbow window)
  versus a 0.59 baseline; bootstrap (stem-cluster) differences -0.20 [-0.38, 0.01], -0.04 [-0.19, 0.11], -0.07 [-0.19, 0.05], -0.01 [-0.12, 0.10].
  Structure falls to the random-neighbour level (0.39-0.41).
- With family exclusion only (no identity cutoff) the gain over baseline is +0.08 to +0.11 and its intervals include 0
  (structure +0.09 [-0.08, 0.26]; FunFam +0.09 [-0.03, 0.21]; elbow window +0.11 [-0.01, 0.24]).
- Accuracy is high when the actual nearest labelled neighbour has identity of at least 0.30 (0.83-1.00, n = 10-12 in 0.30-0.40; 0.92, n about 50 above 0.40) and
  low when it is below 0.25 (0.38-0.50, n = 7-13). The relation is identity-driven in every space, so the spaces do not differ in a way this sample can resolve.
- Label quality (descriptive): restricting queries to experimentally evidenced EC (n = 58) changes accuracies by at most a few points (structure 0.88 versus 0.84
  with any neighbour; 0.66 versus 0.68 family-excluded), within the intervals. Requiring neighbours to be experimental too (n = 51) lowers family-excluded
  accuracy to 0.55-0.67.
- **Against the pre-registered statement** "across families neighbours beat guessing by about ten points": the point estimate (+0.08 to +0.11) is unchanged without an identity cutoff, but its
  interval includes zero, and below 30% identity transfer is no better than guessing. The earlier statement is revised to: function transfer works when a
  neighbour with at least about 30% identity exists and not otherwise.

## 3. A third universe source (script 26)
Reviewed human UniProt entries with InterPro IPR029058 (121) or any of 12 Pfam alpha/beta-hydrolase families (78): the third source contains all 110 genes,
so it is a near-superset and gives no independent recall estimate; no population-size estimate is made. It lists **11 reviewed candidates not in the 110**:
ABHD12B, AADACL2, AADACL3, AADACL4, TEX30, FAM135A, FAM135B, TMCO4, TMEM53, ARB2A (FAM172A) and thyroglobulin (TG). The pre-registered threshold (10 or more)
is met: **110 is a lower bound of the reviewed human set**. The candidates were not added or analysed (they need structure QC, an elbow call and a triad check);
AADACL2-4 and ABHD12B are the most plausible hydrolases. Fetch date recorded in `data/raw/uniprot/third_source_fetched.txt`.

## 4. Structure-derived cores in place of CATH ranges (script 27)
Core = the 200 (also 150, 180) confident residues of the full model nearest the elbow C-alpha; none lie outside the CATH range for the median gene, but cores are
spatial sets with a median of 20 chain segments, so sequence measures on cores are fragmentary. Two genes have fewer than 200 confident residues.
- Core TM distance correlates with the CATH-domain TM distance at Spearman 0.78 (K = 200; 0.73 at K = 150).
- ABHD nearest-neighbour rate (structure): core 0.85 versus CATH 0.80 (difference +0.05, paired bootstrap [0.00, 0.15]); R-precision 0.403 versus 0.400
  (equivalent within the 0.10 margin: yes). With the closest 3.4% pairs removed 0.55 versus 0.50. For the shorter-chain normalisation R-precision is +0.05 [-0.05, 0.16]
  (not shown to be equivalent). The nearest-neighbour comparison cannot be called equivalent because 20 genes give intervals of +-0.15.
- Sequence on cores is weaker (R-precision 0.42 versus 0.48), consistent with fragmentation.
- **Against the pre-registered statement**: the structure-space ABHD signal is not an artefact of CATH trimming or length (R-precision unchanged, point estimates unchanged); the
  nearest-neighbour comparison is underpowered. Cores are centred on the elbow, so this is not independent of the elbow measures.

## 5. Multi-view agreement (script 28)
Spearman correlation of pairwise distances (cross-family pairs only, stem-cluster bootstrap 95% CI): structure-sequence 0.50 [0.38, 0.60]; sequence-FunFam 0.54 [0.44, 0.61];
structure-FunFam 0.36 [0.19, 0.50]; structure versus its own core 0.77 [0.69, 0.84]; core-sequence 0.31 [0.19, 0.41]. All-pair values are within 0.03 of these.
**EC distance is unrelated to every other view** (structure -0.13 [-0.25, 0.01]; core -0.00; sequence 0.03; FunFam 0.07), partly because EC3 is coarse (41 of 80 are 3.1.1) and
because EC separates esterases from peptidases that share the fold. Sequence and FunFam, and EC with sequence/FunFam, are dependent by construction (audit table).
Across-family nearest-neighbour EC3 agreement is 0.675 (structure), 0.70 (core), 0.66 (sequence), 0.675 (FunFam). Structure agrees with EC3 where sequence does not for
LIPJ, ACOT1, ACOT2, ACOT4 and PAFAH2; sequence agrees where structure does not for PPT1, PPT2, NCEH1 and AFMID.
Pre-registered statement "structure and sequence agree moderately" is supported across families.

## 6. Nucleophile calls as three measures (script 29)
On the 31 gold genes: elbow 31/31, sequence motif 26/31, leave-one-out structural transfer from the nearest other gold gene 31/31. By annotation evidence:
experimental (n = 7) 7 / 6 / 7; similarity (14) 14 / 11 / 14; curator (5) 5 / 4 / 5; rule (3) 3 / 3 / 3; none (2) 2 / 2 / 2. The ordering of the measures does not reverse in the
experimental subset. Transfer is not independent of the elbow measure for specificity (all gold genes are positives and mostly paralogs); no latent-class model was fitted.

## What changes in earlier conclusions
1. Function transfer (report 07/00): "about ten points above guessing" becomes "not distinguishable from guessing below about 30% identity; a gain exists only where a neighbour with at least 30% identity is available".
2. The family size: 110 is a lower bound (11 further reviewed candidates).
3. EC is a weak view of relatedness across families; it should not be used as the reference for structure or sequence.
4. Unchanged: the elbow benchmark ordering, the ABHD-name conclusions (core test), structure-sequence agreement.

## Slides that would change if the deck is updated
Function-transfer slide (add the decay curve and identity statement); Question 5 divider answer; limits slide (110 as lower bound); a new slide for the multi-view agreement.

## Caveats
- Identity is a pairwise alignment proxy; 99% of pairs are below 50%, so the decay curve rests on 80 queries with nested, dependent pools.
- 20 ABHD genes limit every ABHD comparison; bootstrap intervals are wide.
- Core sequence measures are fragmentary; the third source is a near-superset of the first two.
