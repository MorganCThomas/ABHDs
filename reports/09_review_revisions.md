# Independent review and revisions (2026-11-12)

An independent reviewer (a separate agent with no access to the authors' reasoning) read the code and reports and recomputed
the headline numbers. Every verified number matched (110 genes, 31/31 vs 26/31 benchmark, 36-row grid, 93 triads, 47/49
experimental, stage 6 statistics). It also found bugs and overstated conclusions. Response, item by item.

## Accepted and fixed
| Finding | Action |
|---|---|
| Incomplete ECs (3.-.-) counted as classes in stage 7; 15 of 92 genes; baseline 0.59 not 0.51; no FunFam baseline | `ec_set` drops incomplete ECs; added leave-family-out, random-neighbour baseline, Wilson CIs; stage 7 re-run and rewritten. Hypotheses now cover 30 genes without a complete EC. |
| Paralog leakage inflates EC transfer | Leave-family-out (same FunFam or paralog stem excluded): EC3 accuracy falls from 0.84 to 0.68 (TM), 0.60 to 0.70 across top spaces; CIs about +-10. |
| seqid clipped at 0: 94% tied at 1.0, degenerate clustering | Replaced by Smith-Waterman local score; stages 6, 6b, 7, 8 re-run. This changed several results (below). |
| "ABHD genes span 4 of 20 clusters" is evidence for, not against, compactness | Rewritten; random-set null across linkages in `review_cluster_null.csv`. |
| ABHD flag misses EPHX3/EPHX4 (= ABHD9/7) and LIPJ/K/M/N | Added `abhd_family` (26 genes) to `universe.csv`; tests repeated under both definitions. ABHD12B is absent from both sources and so from the universe (noted). |
| Stage 8 "signal concentrated at the elbow" | Withdrawn: with a fair sequence baseline, whole-domain local sequence is as good or better (Spearman 0.56 vs 0.53; ARI 0.36 vs 0.34). |
| "95 of 96", not "94 of 95" (report 03) | Corrected. |
| Pocket proxy "no clear difference" | Corrected: residues within 10 A differ (p=0.001) but track domain length (rho=0.49). |
| Elbow RMSD 0.22 A and "31/31" presented as findings | Reworded as training recall and a selection-consistency check; the gold set is mostly similarity-propagated (7/31 experimental). |
| "All 15 ABHD FunFams pure" | Reworded: 14 of 15 have a single human member; ABHD17A/B/C is the only test. |
| "Local similarity reports chemistry not substrate" | Reworded: substrate level was not tested; EC3 is a chemistry label. |
| Neuroligin Gly "register shift" | Dropped: the Gly aligns with the carboxylesterase nucleophile position. ABHD15 (SEGSG) may be shifted. |
| Domain-span length confound | Added `tmshort` (shorter-chain TM norm), length-matched and length-stratified tests. Result: the ABHD signal survives (R-precision 0.40 vs null 0.27, p=0.004). |
| Paralog-driven signal | Added exclusion of near neighbours, quantile-matched across spaces (round 2 fix; the first version used TM < 0.2 in every space): nearest-neighbour rate falls to 0.45 (TM), 0.60 (sequence), 0.71 (FunFam). |
| No CIs | Gene-bootstrap CIs on nearest-neighbour rate and R-precision; Wilson CIs on EC transfer. |

## Net effect on the conclusions
- The central claim stays "partly refuted", with a changed shape: by 20-cluster **sequence** clustering, ABHD genes look cohesive (19 of 20 in
  majority-ABHD clusters at 20 clusters); by **structure**, the signal is carried by close paralog sub-groups and weakens
  once those are excluded. The name is a modest neighbourhood signal, not a clade, and not nothing.
- Function transfer is weaker than first reported and cannot discriminate between spaces.
- The elbow is a robust description of a backbone feature; its recall figure is training recall.

## Rebutted or deferred (not done)
- **Held-out elbow calibration** was run (`review_elbow_heldout.csv`) but is uninformative: gold lies at phi 46-60, psi
  -130..-113, inside every window, so held-out recall is 1.0. A real test needs a non-ABH negative structure set.
- **Re-trimming PREP/PREPL** to the hydrolase fold, **substrate-level labels** (EC4/Rhea), a **matched SOAP null** (other
  strand-to-helix turns), **a window control for stage 8** and **block permutation by paralog group**: deferred.
- The reviewer's suggestion to centre SOAP on the nucleophile side chain is deferred.
- Clustering-level ARI CIs not computed.

## New questions raised (prioritised)
1. ABHD1/2/3/15 and ABHD4/5: active/inactive pairs in the same sub-clade; what pocket and cap features separate them
   (ABHD15 TM distance to ABHD3 is 0.15; ABHD4-ABHD5 0.07)? Needs per-residue conservation within each pair.
2. ABHD17A/B/C vs LYPLA1/2 depalmitoylases: TM distance 0.30, but ABHD17 is closer to ABHD13/14 (0.15-0.17); shared lid
   geometry or convergence?
3. 2-AG hydrolases: MGLL-ABHD6 0.25, MGLL-ABHD10 0.23, MGLL-ABHD12 0.34; compare oxyanion residues and pocket volume
   (ABHD2 I134, ABHD1 I131, ABHD3 L148, ABHD4 W172 from the exploratory oxyanion call). Directly relevant to the ABHD2 case study.
4. Do TM-bearing ABHDs (ABHD1/2/3/6/12/13/14A/16A) cluster by cap size independent of the TM helix? (`structure_qc.n_tm`.)
5. Is EPHX3/4 + LIPJ/K/M/N a bridging clade between ABHDs and epoxide hydrolases?
6. Do the no-nucleophile genes (NDRG, neuroligin, ABHD5, ABHD15) share an interface-like cap, and is the Gly/Asn conserved
   in orthologues?

Reproduce: `uv run python scripts/17_review_controls.py` (after `make study`).
