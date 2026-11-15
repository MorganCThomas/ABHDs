# ABHD superfamily: summary of the executed workflow (stages 1-8)

Design: `docs/study_design.md`. Per-stage reports: `reports/01_*.md` ... `08_*.md`. Figures: `figures/stage*.html`.
Tables: `data/derived/`. Re-run: `make study`. Run on 2026-11-11; revised 2026-11-12 after an independent review (`reports/09_review_revisions.md`).

## Headline results
1. **The UniProt family query is not the superfamily.** Merging it with the CATH FunFam human members gives 110 genes
   (20 ABHD-named): 44 in both, 59 FunFam only (including ABHD4, ABHD5, ABHD12, ABHD13, the carboxylesterases, DPP4, NDRGs,
   neuroligins, EPHX1), 7 UniProt only (including ABHD16A/B, ABHD18). The supplied FunFam list is itself incomplete
   (FF-000064 and FF-000074 contain DAGLB and ABHD16A but are not in the metadata).
2. **The elbow definition recovers known nucleophiles (training recall, not independent validation).** A phi/psi/secondary-
   structure rule, with windows set from known nucleophiles (ABHD2, ABHD14B, LPL, MGLL), recovers 31/31 UniProt-annotated
   nucleophiles; only 7 of those 31 have experimental evidence, 15 are propagated by similarity, and many are pancreatic-lipase
   paralogs. 47 of 49 genes reproduce the turn in experimental structures (ABHD6 7OTS is an S148A mutant, so the turn is a
   backbone feature). The regex recovers 26/31, but it has no Asp; on the 27 Ser/Cys/Thr gold residues it is 26/27 vs 27/27.
   The review's held-out test is uninformative (gold phi 46-60, psi -130..-113 is far inside every window). What the elbow adds
   is specificity: 95 of the 96 genes with an S/C/T/D elbow have exactly one such residue, where the regex flags about three.
   The geometry alone also passes 93 Gly and 88 Ser residues; nucleophile identity comes from the residue type filter.
3. **Catalytic categories.** 77 genes are Ser-His-Asp (18 of 20 ABHD genes plus 59 others); 8 Ser-His-Glu (cholinesterases and
   carboxylesterases), 5 + 1 Asp-based (epoxide hydrolases, MEST, DPP6), 2 Cys-His-Asp, 3 Ser with no His, and 14 with no
   nucleophile at the elbow (neuroligins, NDRGs, DPP10, ACOT6, KANSL3, ABHD5, ABHD15). The elbow turn is present everywhere;
   only the residue changes.
4. **The ABHD name carries a modest, partly length- and paralog-driven neighbourhood signal.** An ABHD gene's nearest neighbour
   is ABHD in 75-85% of cases (17% by chance; bootstrap 95% CI about 60-95%). R-precision is 0.35-0.48 (chance 0.17). The signal
   survives a length-stratified permutation (observed 0.40 vs null 0.27, p=0.004) but drops when close neighbours are
   removed: dropping the closest 3.4% of gene pairs in each space (the share a TM-align distance below 0.2 removes) takes the
   nearest-neighbour rate to 0.45 (TM-align), 0.60 (shorter-chain TM and local sequence), 0.71 (FunFam) and 0.55 (elbow window); an
   unmatched 0.2 cut in sequence space removed only 5 pairs, so the earlier "sequence stays at 0.80" was vacuous. So much of the signal is the ABHD1/2/3, ABHD4/5, ABHD17A/B/C and ABHD14A/B sub-groups. The set is
   not compact (most nearest neighbours are non-ABHD). Calling EPHX3 and EPHX4 (UniProt names "ABHD9" and "ABHD7") and
   LIPJ/K/M/N ABHD-family raises the TM nearest-neighbour rate to 0.85 (26 genes). "ABHD genes occupy 4 of 20 clusters" is
   evidence for clustering, not against it: random 20-gene sets occupy 8.8 on average (p=0.001; linkage-dependent, 2-6).
5. **Local catalytic-site similarity reflects shared chemistry; it was not tested against substrate.** Elbow backbone RMSD
   (0.22 A) is circular because the residues were selected on that geometry; triad skeleton RMSD is 0.58 A. After removing
   incomplete ECs (3.-.-) and excluding neighbours from the same paralog family or FunFam, leave-one-out EC level-3 transfer is
   0.60-0.70 for structure, sequence, FunFam and the elbow window (95% CI about +-10; majority baseline 0.59, random neighbour
   0.40) and 0.52-0.54 for SOAP and triad RMSD. With all neighbours allowed, structure/sequence/FunFam reach 0.83-0.86.
   Paralog leakage accounts for much of the apparent skill; spaces cannot be ranked at n=80. EC level 3 is itself a chemistry
   label, so substrate-level transfer remains untested.
6. **Elbow-anchored windows are as informative as whole-domain sequence, not better.** With a corrected sequence baseline
   (local Smith-Waterman instead of a clipped global score, which gave 94% tied distances), the whole-domain local score
   correlates with FunFam distance at 0.56 vs 0.53 (elbow window) and 0.49 (catalytic windows), and ARI at 30 clusters is 0.36
   vs 0.34. A claim that the signal is concentrated at the elbow is not supported.
7. **Function transfer must be gated on the catalytic category.** Neuroligins resemble carboxylesterases (nearest neighbours
   BCHE/CES/ACHE) but have no nucleophile, so an inherited EC would be wrong; 12 of the 30 genes without a complete EC were flagged.

## Follow-ups (reports/10_followup_experiments.md)
- In the AlphaFold models ABHD15 and ABHD5 lack only the nucleophile: His and an acid sit at the aligned positions and ABHD15 Ser238
  is 6.3 A from the His (a register shift is not supported by geometry). After orthologue clean-up (reports/12): ABHD5 Asn is
  434/434; ABHD15 Gly is present in most amniotes and amphibians but Ser in crocodilians and a few basal fish (polarity unresolved,
  no tree). The ABHD4/ABHD5 difference near the pocket is a conserved ABHD4-specific insertion (ABHD4 180-198), not ABHD5 remodelling.
  Like-for-like cavities: ABHD15 < ABHD3, ABHD2 (18/18 variants); ABHD5 > ABHD4 (18/18).
- Oxyanion hole: Nu+1 N is a near-universal first donor (2.8-3.7 A); the second donor residue varies (Ala, Phe, Tyr, Asn, Ile, Leu).
  Pocket volumes depend on model (trimmed vs full-length, rho 0.6-0.7) and anchor; only contrasts consistent across 18 variants are
  reported (ABHD6 and ABHD5 largest; ABHD2 < ABHD6; ABHD12 < MGLL; ABHD2 vs MGLL unordered). No selectivity claim is made.

## Changes from the plan
- Foldseek (planned) was not used: downloading and running its binary was blocked in this environment; TM-align via the
  `tmtools` Python package was used instead.
- Domains use CATH ranges from the FunFam alignment headers (union across FunFams) rather than a pLDDT heuristic, except for
  5 genes with no alignment entry.
- "Nullified" sites are reported as "no nucleophile at elbow", without orthologue conservation support (not done).
- Pocket descriptors: a composition proxy (stage 5) plus a grid-cavity calculation (follow-up) that is not robust to model, anchor and clearance.
- Oxyanion-hole calls are exploratory and unvalidated.
- The attrition ledger exists for stages 1 (444 accessions) and 2-4 (counts by flag/category), not as a single diagram.

## Limits
- Human only; 110 genes. 31 genes have an annotated nucleophile for benchmarking; EC labels exist for 92.
- Domain spans from CATH headers include non-hydrolase domains for some genes (PREP/PREPL span the beta-propeller, 520-700
  residues); TM-align distance is partly length-driven (a shorter-chain-normalised variant, `tmshort`, was added).
- Distances use different domain lengths; ProtBert is mean-pooled and weakly correlated with everything else.
- Gene-bootstrap CIs are given for the label statistics and Wilson CIs for EC transfer; ARI differences still have none.
- Permutation tests treat the 20 ABHD genes as independent; they form about 8-10 paralog groups.
- Literature statements in the reports (for example DPP6 being inactive) come from background knowledge and were not checked.

## Suggested next steps
1. Ask Christina why FF-000064/FF-000074 are missing from the metadata and whether ABHD16B, ABHD18, SERHL, SERHL2 and EPHX4
   are in CATH 3.40.50.1820 at all.
2. Orthologue conservation is done for ABHD15 and ABHD5 only; extend it to the other 12 no-nucleophile genes (see `docs/next_steps.md`).
3. ABHD15 register offset: not supported by geometry (follow-up A1/A5); neuroligin Gly aligns to the GESAG nucleophile position.
4. Review follow-ups not yet done: re-trim PREP/PREPL to the hydrolase fold, substrate-level (EC4/Rhea) labels, a matched SOAP null, a non-ABH negative set for the elbow, window controls for stage 8. See `reports/09_review_revisions.md`.
5. Link the case study: ABHD2's neighbours (ABHD3, ABHD1) share its Ser-His-Asp site and elbow; use pocket geometry, not the
   elbow, to look for differences relevant to the ligand series.
