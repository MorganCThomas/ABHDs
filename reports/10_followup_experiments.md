# Follow-up experiments (2026-11-12)

> **Superseded in part by `reports/11_review2_response.md` and `reports/12_cavity_orthologue_redo.md`.** The pocket-size, cavity-volume and orthologue statements below were revised there: ABHD15's cavity is smaller than ABHD3's/ABHD2's, ABHD5's larger than ABHD4's, ABHD15 Gly is not universal across vertebrates, and the ABHD4 180-198 difference is an ABHD4-specific insertion.

Design (hypotheses, controls fixed in advance): `docs/followup_design.md`. Scripts: `18_pocket_oxyanion.py`,
`19_inactive_vs_active.py`, `20_panel_pockets.py`. Tables: `data/derived/followup_*.csv`, `pocket_*.csv`, `oxyanion_donors.csv`;
figure `figures/followup_b_pockets.html`. All structures are AlphaFold apo models; no ligand or activity data were used, and
literature recollections are not used as evidence.

## Question A: what separates inactive ABHD15 and ABHD5 from their active relatives?

| Experiment | Result |
|---|---|
| A1 alternative nucleophile | No S/C/T/D at elbow -3..+3 forms an ordered triad in any of the 14 no-nucleophile genes. Control: 38 off-register S/C/T/D candidates in nucleophile genes, 0 with a triad (0%). ABHD15 Ser238 (i+1) OG is 6.3 A from His391 NE2 (ABHD3 Ser220: 2.8 A), while Gly237 CA is 3.6 A from it. **The register-shift hypothesis for ABHD15 is rejected**: the Gly occupies the nucleophile position. |
| A2 catalytic positions | TM-align places ABHD15 His391 and Asp360 on ABHD1/2/3's His and acid (CA displacement 0.5-2.0 A), and ABHD5 His327 and Glu177 on ABHD4's His320 and Asp170 (0.2-0.7 A). **Both lost only the nucleophile**; His and an acid are present. ABHD5's acid is Glu where ABHD4 has Asp. |
| A4 orthologues | 359 aligned ABHD15 and 277 ABHD5 vertebrate sequences (>=40% identity, >=70% coverage, first entry per organism name; ~64% birds, 99% unreviewed). ABHD15 Gly at the elbow is present in 351/359 (8 Ser, 6 of them Ser-Thr in basal fish and crocodilians); Ser238 347/359; His391 100%, Asp360 99.7%. ABHD5 Asn153 276/277 (one Ser). In the active paralogs the nucleophile Ser is 98-100% conserved. **The Gly/Asn is present across the sampled vertebrates, so the human model is not an outlier**; the lineage history (when the Ser was lost) is unresolved because no tree or outgroup was used. His/acid sit at the 91-93rd percentile of domain positions (strict-less-than percentiles with ties at 1.0; the elbow Gly itself is only 0.67). Orthologue sets contain some contaminants (an ABHD5 "Ser" entry is closer to ABHD4; an ABHD4 chicken entry is ABHD5) and the UniProt query was not paginated (ABHD2/3/5 truncated at 500). |
| A3 pocket-shell divergence (CA displacement > 3 A or unaligned, residues within 6 A of the reference cavity) | **ABHD5 vs ABHD4: 9/28 shell residues diverged vs 15/282 outside (odds ratio 8.4, residue-level Fisher p=1e-4), but this does not survive a run-aware test (circular-shift p=0.14; with ABHD5 as reference odds ratio 2.6, circular p=0.44)**, because it comes from one 19-residue ABHD4 segment (180-198) with 6 shell residues. **ABHD15 vs ABHD1/2/3**: the shell is, if anything, less diverged than the rest (odds ratios 0.14-0.46) but no pair is significant by the run-aware test (p 0.11-0.86) and the active-active pair ABHD1/ABHD2 gives the same odds ratio of 0 (Fisher p=0.020, circular p=0.13). Active-active baseline otherwise 1.2. A shared ABHD3-specific C-terminal segment (383-404, absent in ABHD1/ABHD2) lines the pocket (6 shell residues). |
| Cavity volume | Loss of the nucleophile did not collapse the pocket: ABHD15 162 vs ABHD3 263 (clearance 3.2 A), ABHD5 270 vs ABHD4 185 (box-truncated for ABHD5 and ABHD3, see caveats). |

Interpretation. In the models both proteins keep His and an acid at the aligned positions and lack only the nucleophile; the register-shift
alternative for ABHD15 is not supported by CA and side-chain geometry (the OG rotamer of a non-catalytic Ser is itself predicted, and the
scan has no positive control). ABHD4 segment 180-198 is displaced in the ABHD5 model and is the one localised difference worth testing, but
the statistics do not support calling the pocket remodelled. Whether the retained sites bind ligands or partners is not addressed.

## Question B: the monoacylglycerol-hydrolase panel (MGLL, ABHD6, ABHD10, ABHD12, ABHD16A, ABHD2 and relatives)

| Experiment | Result |
|---|---|
| B1 oxyanion hole | The Nu+1 backbone N is 2.77-3.73 A from the nucleophile atom (panel 2.89-3.01; 11 genes outside 2.85-3.05), a near-universal first donor; a 3.6 A cutoff finds a second donor in only 5 genes. The second donor is the residue after the upstream HG-type motif, 2.8-5.4 A away (median 4.7); 11 genes (mostly S9 peptidases, which use a Tyr donor) have none; its atoms line the cavity in 89% of the 85 that do. The stage-3 oxyanion call is therefore the second donor and is plausible by geometry; its residue identity differs by clade (below). |
| Second-donor identity | MGLL Ala51, ABHD6 Phe80, ABHD10 Tyr84, ABHD12 Asn177, ABHD16A Asn288, ABHD2 Ile134, ABHD1 Ile131, ABHD3 Leu148, ABHD4 Phe77, ABHD13 Asn123, ABHD14B Ile41, ABHD17A Asn120, LYPLA1 Leu30. Asn is the second donor in 7 of 18 ABHD-named nucleophile genes vs 3 of 67 others (Fisher p=0.0005), but the seven are ABHD17A/B/C, ABHD16A/B, ABHD12 and ABHD13 (4 stems); collapsing paralog stems gives 4/14 vs 3/53 (p=0.03), 14 residue types were screened, and the 11 genes without a donor are excluded non-randomly. **Treat the enrichment as suggestive.** ABHD1/2/3 have Ile/Ile/Leu; MGLL has Ala and ABHD6 Phe. The donor acts through its backbone NH, so a side-chain difference is not evidence of functional difference. |
| B2 pocket volume (clearance 3.6 A; 3.2 A in brackets) | ABHD6 103 (219), MGLL 93 (199), ABHD10 75 (243), ABHD4 75 (185), ABHD2 64 (103), ABHD3 44 (263), ABHD16A 43 (128), ABHD13 42 (97), ABHD1 27 (58), ABHD12 15 (52), LYPLA1 9 (27), ABHD17A 4 (15), ABHD14B 1 (4). Volume correlates with domain length (rho 0.36, p=1e-4). Rankings across clearances agree genome-wide (rho 0.87-0.97) but **panel values are not robust**: with the full-length model instead of the CATH-trimmed domain, ABHD2 at 3.6 A goes from 64 to 148 (larger than MGLL 93 and ABHD6 98), ABHD15 162/38 to 85/2, ABHD17A 15/4 to 302/114; truncating nucleophile side chains to Gly (like-for-like with ABHD15/ABHD5, which are anchored on CA) changes volumes by up to 2x. **No panel ordering of pocket size is claimed.** |
| B3 pocket-lining residues | Mapped through TM-align, MGLL's lining residues are no more conserved in any panel partner than the aligned core in general (e.g. ABHD6 0.20 vs core 0.22; ABHD12 0.07 vs 0.22), and across all 78 pairs pocket identity exceeds the core-position 95th-percentile null in only 10. **A shared pocket-lining signature of the 2-AG hydrolases was not detected, with low power** (13-49 positions per pair, TM-align RMSD up to 4.3 A). There is mild pocket enrichment overall (10 of 78 pairs vs 3.9 expected, mostly among ABHD13/14B/17A/LYPLA1 pairs). |

For the ABHD2 case study: ABHD2's second donor is Ile134, like ABHD1 (Ile131) and ABHD3 (Leu148). No statement about inhibitor
selectivity versus MGLL or ABHD6 follows from these apo-model data (the donor acts through its backbone NH and panel pocket sizes
are not robust).

## Caveats
- (Round 2 review) Cavity volumes are model-dependent: see the interpretation of B2. Cavity volumes come from a 1 A grid in a +-14 A box: cavities that reach the box edge are truncated (34 genes at 3.2 A, 17 at
  3.6 A, including ABHD6 and ABHD5 at 3.6). Volumes are relative descriptors of apo AlphaFold models, not binding-site
  volumes, and tiny values (ABHD14B, ABHD17A) may reflect closed models.
- The anchor for genes without a nucleophile (ABHD15, ABHD5) is the elbow CA, not the side chain.
- Orthologues come from a gene-name UniProt query; one sequence per species is kept, but unreviewed entries may include
  mislabelled paralogs (filtered by >=40% identity and coverage only).
- TM-align mapping of pocket positions is only as good as the structural alignment (RMSD 1.0-3.5 A across panel pairs); small
  pockets (ABHD14B, LYPLA1, 10-17 positions) give noisy identities.
- Divergence thresholds (3 A, 6 A shell) were fixed before running, but results were not tested for sensitivity to them.
- Oxyanion-hole "second donor" was defined by geometry and motif position; its function is not validated.

## Open follow-ups
0. Redo cavities with full-length models, truncated side chains, rotated grids and a larger box; reassess every pocket statement. 1. Align the ABHD4 180-198 / ABHD5 equivalent segment across orthologues and test whether it is conserved (partner-binding surface?).
2. Dock or place a monoacylglycerol in MGLL, ABHD6, ABHD2 and ABHD12 models to test whether pocket shape, not residue identity, explains selectivity; check against inhibitor SAR for ABHD2.
3. Extend the A4 orthologue analysis to the eight ABHD15 orthologues with Ser at the elbow (possible reversions).
