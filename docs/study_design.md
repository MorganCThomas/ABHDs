# ABHD Superfamily Study Design

Drafted 2026-11-11. Living copy: https://claude.ai/artifact/2CQW1MtDCtbFq2eVRXoZbS (this file is a snapshot).

## Question and central claim

**Question.** What does the α/β-hydrolase (ABH) superfamily look like from a structural and catalytic-site perspective, and how well does the "ABHD" gene name describe it?

**Central claim to test.** The ABHD name does not mark a structural or sequence clade. Many proteins outside the ABHD naming (lipases, epoxide hydrolases, neuroligins, DPP4, LCAT and others) sit in the same fold superfamily, and some ABHD-named proteins are closer to non-ABHD proteins than to other ABHDs.

**How it is tested.** Define the superfamily from two independent sources, call each member's catalytic machinery from predicted structure, then ask whether structure-based groupings agree with sequence-based ones (CATH FunFams) and with the ABHD label. The claim is supported if label purity of the ABHD genes is low in every grouping, and refuted if ABHD genes form compact clusters.

## Pipeline overview

Eight stages, each with a gate that writes an attrition table (what was dropped and why) so every outlier can be traced to a reason.

1. **Universe.** Merge UniProt-annotated ABHD entries with CATH FunFam members (3.40.50.1820), keeping source labels.
2. **Structures.** Fetch AlphaFold models, trim to the hydrolase domain, apply quality filters.
3. **Catalytic-site calls.** Detect the nucleophile elbow, then the triad and oxyanion hole.
4. **Categories.** Classify each protein by its catalytic machinery, including nullified sites.
5. **Local comparison.** Compare aligned catalytic-site environments (RMSD, SOAP, pocket descriptors).
6. **Global comparison and FunFam concordance.** Foldseek, FunFam HMM similarity and a protein language model, compared to each other and to the ABHD label.
7. **Function inference.** Transfer substrate and EC annotations between neighbours, validated by leave-one-out.
8. **Realignment.** Rebuild sequence alignments anchored on the structural elbow and catalytic residues, and re-test stage 6.

## Stage-by-stage design

Each row states the method, the measures that decide success, and what the gate reports about outliers.

| Stage | Method | Metrics and controls | Gate: outliers reported |
| --- | --- | --- | --- |
| 1. Universe | Query UniProt for ABHD entries and take human members of CATH 3.40.50.1820 FunFams; collapse to one canonical sequence per gene; tag each protein with its source(s). Optional third source: ESTHER or InterPro. | Venn of UniProt-only, FunFam-only and both; ABHD-named vs not, per source. | Isoforms, cDNA fragments and readthrough genes dropped, each with a reason code. |
| 2. Structures | AlphaFold DB models; trim signal peptides, transmembrane helices and disordered tails to the hydrolase domain. | pLDDT per domain; check geometry against experimental structures (ABHD6, ABHD14B). | Proteins with no model, truncated models, pLDDT below cutoff. |
| 3. Catalytic-site calls | Nucleophile elbow by phi/psi (about +53 / -120, preceded by beta, followed by helix); triad by distance plus sequence order (nucleophile, acid, His); oxyanion hole. Loose regex only as a sanity check. | Precision and recall against UniProt active-site annotations and M-CSA; agreement with the regex; sensitivity to thresholds. | No elbow, elbow without His or acid, multiple elbows. |
| 4. Categories | Classify: Ser-His-Asp, Ser-His-Glu, Asp-His-Asp (epoxide hydrolase type), Cys-based, nullified (for example ABHD5 Asn), no evidence. | Counts per category, split by ABHD-named vs not and by FunFam. A nullified call needs conservation support from FunFam alignments. | Separate truly lost from not detectable. |
| 5. Local comparison | Superimpose elbow, triad and oxyanion-hole atoms; compute local RMSD; SOAP on the aligned site; pocket descriptors (volume, polarity) for specificity. | Null baseline from SOAP on random Ser residues; stability to r_cut, n_max, l_max. Nucleophile similarity reports chemistry, not substrate. | Proteins whose site does not superimpose. |
| 6. Global comparison and FunFam concordance | Foldseek or TM-align all-vs-all on trimmed domains; FunFam HHalign E-values (as -log10); ProtBERT embeddings as a third view. | ARI or NMI between clusterings; rank correlation of structure vs sequence distance; ABHD label purity per cluster; list of discordant pairs. | Structurally close but sequence-distant pairs, and the reverse. |
| 7. Function inference | Transfer EC, Rhea and ChEMBL annotations across neighbours in the stage 5 and 6 spaces. | Leave-one-out: hide each known substrate or EC and score recovery against a nearest-FunFam baseline. | Predictions with no annotated neighbour within the distance cutoff. |
| 8. Realignment | Structure-guided alignment anchored on elbow and catalytic residues; rebuild trees or clusters. | Compare with the FunFam alignment: clustering agreement, branch support, control for length and pLDDT. | Members that move between groups. |

## Cross-cutting rules

Four practices apply at every stage and make the result defensible.

- **Attrition ledger.** Every dropped protein gets a reason code at the stage that dropped it. The ledger is plotted as a flow diagram from stage 1 to stage 8 and doubles as the outlier list.
- **Gold-standard benchmark.** Catalytic-site calls are scored against UniProt active-site features and M-CSA, and geometry is checked against experimental structures where they exist.
- **Sensitivity analysis.** Vary the phi/psi windows, pLDDT cutoff, triad distances and SOAP parameters, and report which category assignments change.
- **Reproducibility.** Each stage is a script run from a Makefile target, writes a versioned table, and keeps the source-label mapping back to UniProt and CATH. Notebooks only plot and interpret.

## Planned figures and case study

| Figure | Shows | Stage |
| --- | --- | --- |
| 1 | Universe: overlap of UniProt ABHD and FunFam members, with the ABHD-named share of each | 1 |
| 2 | Attrition flow from raw entries to analysed proteins, with reason codes | 1-4 |
| 3 | Elbow geometry: phi/psi of known nucleophiles, and benchmark precision and recall | 3 |
| 4 | Catalytic-site categories across the superfamily, coloured by ABHD name | 4 |
| 5 | Structure-based tree or embedding with tracks for ABHD name, FunFam and category | 5-6 |
| 6 | Concordance of structure and sequence groupings; discordant pairs highlighted | 6 |
| 7 | Function-inference accuracy from leave-one-out | 7 |
| 8 | Effect of elbow-anchored realignment on groupings | 8 |

**Case study: ABHD2.** The ABHD2 work (front and back pockets, oxyanion hole at I134, covalent inhibition by a compound series) closes the paper. It shows how the family-level picture, such as which neighbours share its catalytic chemistry and which pocket features differ, informs ligand design for one member.

## Risks and open decisions

- [ ] **Universe scope.** Human only, or add orthologs from other species? Orthologs improve subfamily inference and nullified-site calls but enlarge the universe.
- [ ] **Third source.** Include ESTHER or InterPro as an independent family definition, or keep to UniProt and CATH.
- [ ] **Gold standard.** Confirm which active-site annotations count as ground truth (UniProt features, M-CSA, literature) before stage 3 is scored.
- [ ] **Specificity vs chemistry.** Nucleophile-centred similarity reflects catalytic chemistry, so decide which pocket descriptors stand in for substrate specificity in stage 5.
- [ ] **Model limits.** AlphaFold models are apo and may misplace active-site side chains; agree how many experimental structures to use for validation.
- [ ] **Fragment handling.** Decide whether truncated isoforms are dropped at stage 1 or kept as a labelled low-confidence group.

## Execution notes (2026-11-11)

All eight stages were run; see `reports/00_summary.md`. Departures from this plan:

- The central claim was **partly refuted**: ABHD-named genes are usually each other's nearest neighbours (75-85% in structure,
  sequence and FunFam space), but they are not a compact group (R-precision 0.40 vs 0.17 by chance) and share structural
  clusters with many non-ABHD proteins. The stage 6 wording above ("supported if label purity is low in every grouping") should
  be read with that result.
- Foldseek was not available (binary download blocked); TM-align (tmtools) was used.
- Domain trimming uses CATH ranges from FunFam alignment headers.
- After an independent review (`reports/09_review_revisions.md`): the sequence distance was corrected (local alignment, not a
  clipped global score), EC transfer was re-done without incomplete ECs and with leave-family-out, and the stage 8 claim that
  the elbow concentrates the signal was withdrawn. The ABHD name is a modest, partly paralog- and length-driven neighbourhood
  signal (19 of 20 ABHD genes sit in majority-ABHD sequence clusters, but in every space the nearest-neighbour signal falls once the closest pairs are removed).
- Not done: orthologue conservation for nullified-site calls, cavity-based pocket descriptors, confidence intervals, a single
  attrition flow diagram.

- Labels as measurements (reports/13_labels_as_measures.md): CATH, FunFam, EC and UniProt annotations were re-examined as measurements; function transfer is not distinguishable from guessing below about 30% identity; 110 genes is a lower bound (11 further reviewed candidates); EC distance is unrelated to structure or sequence across families.
