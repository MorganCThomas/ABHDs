# Stage 8: alignment anchored on the catalytic residues

Script: `scripts/15_anchored_alignment.py`. Outputs: `data/derived/dist_anchor_elbow.csv`, `dist_anchor_cat.csv`,
`anchor_windows.csv`, `anchored_summary.csv`, `anchored_clusters.csv`; figure `figures/stage8_anchored.html`.

## What was done
Residues were placed in fixed columns by structural role, so no alignment search was needed: a +-15-residue window around
the elbow residue (110 genes) and a concatenation of +-7 windows around nucleophile, acid and His in chain order (93 genes
with a triad). Distance = 1 - ungapped BLOSUM62 score, normalised by self scores.

## Results (revised after review)
The whole-domain sequence baseline was replaced: the original BLOSUM62 global score was clipped at zero, so 94% of pairs were
tied at distance 1.0 and the "0.40" baseline was a degenerate comparison. The baseline is now a Smith-Waterman local score.
- Rank correlation with FunFam distances: whole-domain local sequence 0.56, elbow window 0.53, catalytic windows 0.49
  (the last on 93 genes, not the same gene set).
- Agreement with FunFam (adjusted Rand), 30 clusters: catalytic windows 0.39 (87 genes), whole-domain local sequence 0.36,
  elbow window 0.34, TM-align 0.25 (103 genes). At 20 clusters: 0.25, 0.20, 0.22, 0.14.
- EC level-3 transfer, leave-family-out: elbow window 0.70, catalytic windows 0.67, local sequence 0.66, TM-align 0.68; not
  distinguishable (CI about +-10).
- ABHD label: catalytic windows give an ABHD nearest neighbour for 67% of ABHD genes (33% skipping paralogs), the elbow
  window 75% (65%); R-precision 0.41 and 0.29.

## Interpretation
A 31- or 63-residue structurally defined window is as informative about FunFam relatedness and function as whole-domain
sequence, with no clear advantage. The earlier statement that the family signal is concentrated at the catalytic elbow is not
supported.

## Caveats
- Windows are ungapped; insertions inside a window shift columns, which blurs distant pairs.
- Differences in ARI between spaces are small and were not given confidence intervals. No window control (same width at a non-catalytic conserved strand) was run, and gene sets differ between rows (93, 103, 110).
