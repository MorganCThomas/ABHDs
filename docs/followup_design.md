# Follow-up experiments from the review (2026-11-12)

Two questions from `reports/09_review_revisions.md`, with hypotheses, measures and controls fixed before running.

## Question A: what separates the inactive ABHD15 from the active ABHD1/2/3 (and ABHD5 from ABHD4)?
Background: ABHD15 has Gly at the elbow (SEGSG, Ser at i+1) and sits 0.15 TM-distance from ABHD3; ABHD5 (Asn) is 0.07 from ABHD4.
| Experiment | Test | Control / null |
|---|---|---|
| A1 alternative nucleophile | For every no-nucleophile gene, test S/C/T/D residues at i-3..i+3 of the elbow for a His+acid triad (chain order, distance cutoffs of stage 3). H: ABHD15 Ser i+1 is a register-shifted nucleophile (triad present) vs a genuine loss (no triad). | Same scan on the 96 genes with a nucleophile: false-positive rate of off-register triads. |
| A2 catalytic-position mapping | TM-align pairs (ABHD15 vs ABHD1/2/3; ABHD5 vs ABHD4); read off the residue the partner has at each catalytic position. H: loss is of the nucleophile only (His/acid retained) or of the whole triad. | Active-active pairs (ABHD1/2/3 against each other). |
| A3 structural divergence | Per-residue CA displacement after superposition and unaligned insertions, mapped to distance from the pocket. H: the inactive member diverges mostly in the cap near the pocket, more than active-active pairs do. | Active-active pairs; far-from-pocket residues as the background. |
| A4 orthologue conservation | Vertebrate orthologues from UniProt (>=40% identity to the human protein, >=70% coverage). H: if the elbow Gly (ABHD15) / Asn (ABHD5) is conserved at a rate above the domain-wide rate, the loss is ancient and not a human-model artefact. | Domain-wide per-position identity; same positions in active paralogs (Ser conserved?). |

## Question B: the 2-AG / monoacylglycerol hydrolase panel (MGLL, ABHD6, ABHD10, ABHD12, ABHD16A, ABHD2) and its relatives
Panel adds LYPLA1 and ABHD17A as depalmitoylase outgroups, ABHD1/3/4/13/14B as ABHD controls.
| Experiment | Test | Control / null |
|---|---|---|
| B1 oxyanion hole | Backbone-N donors near the nucleophile O (nearest non-local N and local Nu+1/Nu+2), their residue and the motif context. The stage-3 call (nearest N within 6 A, 4.4-5.2 A) may be loose; report both a 3.6 A cutoff and the motif-positioned second donor. (Result: Nu+1 N is the universal first donor; the stage-3 call is the second donor and is geometrically plausible - see reports/10.) | Distribution over all 96 nucleophile genes. |
| B2 pocket volume | Grid cavity connected to the active site (probe clearance 3.2 A, buriedness >=70%), sensitivity to clearance 3.0/3.6; lining composition. | All 110 genes; volume vs domain length (known confound) with residuals z-scored. |
| B3 pocket-position comparison | Map each panel gene's lining residues onto MGLL through the TM-align alignment and tabulate identity/class at each position. H: ABHD6 and MGLL share a pocket-lining pattern that ABHD12 and ABHD10 do not. | Lining-residue identity between random non-panel pairs at the same TM distance. |
Everything is on AlphaFold apo models; pocket volume is a relative descriptor, not a binding-site measurement. Literature statements are not used as evidence.
