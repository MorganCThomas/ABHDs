# Next steps (noted 2026-11-13, after follow-up experiments e214083)

Carried forward from `reports/10_followup_experiments.md` and `reports/09_review_revisions.md`.

## From the follow-up experiments
1. ~~ABHD4 180-198 conservation~~ done: conserved ABHD4-specific insertion absent from ABHD5. Original: **ABHD4 180-198 / ABHD5 equivalent segment.** Align it across vertebrate orthologues of both genes and test whether it is
   conserved in ABHD5 (a partner-binding or lipid-droplet surface?) and divergent from ABHD4. Needs the existing orthologue
   sets (`data/raw/orthologues/`); the question is whether the 9/28 shell divergence is localised to this segment.
2. **Ligand placement in the 2-AG panel.** Dock or restrain a monoacylglycerol in MGLL, ABHD6, ABHD2 and ABHD12 to test
   whether pocket shape, not residue identity, separates them; compare with inhibitor structure-activity data for ABHD2.
3. ~~ABHD15 reversions~~ partly done (crocodilians and basal fish are Ser); the triad in those models is untested. Original: **ABHD15 reversions.** Inspect the 8 vertebrate ABHD15 orthologues with Ser at the elbow (taxonomy, triad geometry).

## Still open from the first review
4. Re-trim PREP/PREPL (and other long-domain genes) to the hydrolase fold and re-run the length-sensitive tests.
5. Substrate-level labels (EC4, Rhea) to test whether local-site similarity predicts substrate, not just chemistry.
6. A matched SOAP null (other strand-to-helix turns) and a window control for stage 8.
7. A non-ABH negative structure set for the elbow definition, since the 31/31 benchmark is training recall.
8. Paralog-block permutation tests and bootstrap intervals for clustering agreement (ARI).
9. Orthologue conservation for the other 12 no-nucleophile genes (neuroligins, NDRGs, DPP10, KANSL3, ACOT6).
10. Add ABHD12B (Q7Z5M8) to the universe and ask Christina about FF-000064/FF-000074 and the five genes absent from the
    FunFam alignments.

## Added after the second review (reports/11_review2_response.md)
0. ~~Redo cavities~~ done (reports/12). Remaining: pocket-shell tests on the new cavities, a finer box for ABHD5/ABHD6, and a cavity-vs-segment test for ABHD4 180-198.
11. Run-aware shell tests with all active-active pairs in both directions; sensitivity to D (2-5 A) and S (4-8 A).
12. ~~Orthologue hygiene~~ done (reports/12). Remaining: a small phylogeny (or ancestral-state test) for ABHD15 Gly vs Ser; check the triad in Ser-bearing crocodilian/fish ABHD15 models; extend to other no-nucleophile genes.
13. Oxyanion donor test on paralog stems with FDR; handle S9 Tyr-donor genes explicitly; positive control for the alternative-nucleophile scan.
14. Leave-family-out by FunFam plus identity threshold; paired tests vs the 0.59 majority baseline.
15. Engineering: Make targets for structures/domains/cavities, shared code in `src/abhd`, derive residue numbers from outputs, add a smoke test, run `make study` end to end.

## Follow-up questions (final slide of the deck)
1. Models versus experiments: compare AlphaFold-based calls with PDB structures; quantify active-site variability within and between structures.
2. Chemistry of the triad: QM/MM or machine-learned potentials for nucleophile reactivity and oxyanion stabilisation across the family.
3. Known substrates: census of known substrates and how they map onto structural, pocket and sequence groups.
4. Predicting the unknown: substrate or activity assignment for unannotated proteins, gated by catalytic category.
5. Ligands and dynamics: docking or simulation of lipid substrates in the lipid-hydrolase panel.
6. Roles of the inactive proteins: shared surfaces or caps across the 14 no-nucleophile proteins (partner binding?); Ser-restoring mutants of ABHD15 and ABHD5.

## From the labels-as-measures analysis (report 13)
- Add the 11 candidate genes (ABHD12B, AADACL2/3/4, TEX30, FAM135A/B, TMCO4, TMEM53, ARB2A, TG) with structure QC and elbow calls; rerun the ABHD tests.
- Use experimental-only EC or Rhea reaction labels, and a finer functional label than EC3, to test transfer above and below 30% identity.
- Update the presentation (function slide, Question 5 answer, limits slide, new multi-view slide).
