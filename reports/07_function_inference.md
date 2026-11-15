# Stage 7: function inference

Script: `scripts/14_function_inference.py`. Outputs: `data/derived/function_loo.csv`, `function_hypotheses.csv`; figure
`figures/stage7_loo.html`.

## Leave-one-out EC transfer (revised after review)
Original run counted incomplete ECs (3.-.-, 3.1.-) as classes, so two genes with different unknown enzymes "agreed". They are
now dropped at levels they do not define (80 genes at level 3, 92 at level 1). The table gives nearest-annotated-neighbour
accuracy at EC level 3 with Wilson 95% CI. "Leave-family-out" additionally excludes neighbours with the same FunFam or the same
paralog stem (CES, DPP, NLGN, ABHD17...). Baselines: majority class 0.59 (3.1.1); random labelled neighbour 0.40.

| Space | all neighbours | leave-family-out | 95% CI (family-out) |
|---|---|---|---|
| Anchored elbow window | 0.86 | 0.70 | 0.59-0.79 |
| TM-align | 0.84 | 0.68 | 0.57-0.77 |
| Sequence (local SW) | 0.85 | 0.66 | 0.55-0.76 |
| FunFam | 0.84 (77 genes) | 0.66 | 0.55-0.76 |
| Anchored catalytic windows | 0.85 | 0.67 | 0.56-0.76 |
| TM-align, shorter-chain norm. | 0.83 | 0.60 | 0.49-0.70 |
| SOAP at elbow | 0.70 | 0.54 | 0.43-0.64 |
| Triad RMSD | 0.76 | 0.52 | 0.41-0.63 |
| ProtBert | 0.65 | 0.53 | 0.42-0.63 |

## Interpretation
Across families, EC level-3 transfer from the best spaces is only about 10 points above the majority baseline and the intervals
overlap it; most of the apparent skill comes from near-identical paralogs (CES, DPP, EPHX, LIP, NLGN). The spaces cannot be
ranked at this n. Local catalytic-site similarity is at or below baseline, consistent with it describing shared chemistry,
but EC level 3 is a chemistry label, so this does not test substrate specificity.

## Hypotheses for the 30 genes without a complete EC number
(Revised: 18 had no EC at all; 12 more carry only partial ECs such as 3.-.-.) Transfer is flagged **not applicable** for 12
genes whose elbow carries no nucleophile (ABHD15, the five neuroligins, DPP10, NDRG1-4, KANSL3): they resemble carboxylesterases
or peptidases structurally (neuroligins -> BCHE/CES/ACHE) but lack the catalytic residue, so inheriting an EC would be wrong.
Predictions for the remaining 18 come from the 3 nearest TM-align neighbours and are hypotheses only: given the family-out
accuracy above, treat any prediction whose neighbours are all one paralog family as weak. DPP6 is listed as applicable by the
rule, but I believe it is catalytically inactive (not verified here).
