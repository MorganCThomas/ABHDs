# Presentation plan v2: "What's in a name? Mapping the alpha/beta-hydrolase superfamily from structure"

Revision (2026-11-14): v1 left out about half of the analyses. v2 includes every result and conclusion, grouped into eight
overarching questions. Each question has a divider slide (question + short answer), then slides in the order
**method -> result(s)**; a section may have one method with several results, or alternate method / result pairs.
Every slide carries a small tag ("Q4 · METHOD", "Q4 · RESULT") and a title that states its message. Audience: general scientific;
short bullets or one-sentence boxes; figures over text; one message per slide; numbers read from the committed tables.

## Contents (one slide) and sections
| Q | Question | Short answer |
|---|---|---|
| 1 | What is in the human alpha/beta-hydrolase family? | 110 genes; no single database holds them all |
| 2 | Can the backbone alone locate the nucleophile? | Yes: 31/31 known sites; 14 proteins have lost it |
| 3 | How alike are catalytic sites across the family? | Very alike: shared chemistry, not shared substrate |
| 4 | Does the ABHD name mark a clade? | A neighbourhood, not a clade |
| 5 | Can neighbours predict function? | Weakly once paralogs are excluded |
| 6 | Why are ABHD15 and ABHD5 inactive? | They lost only the nucleophile |
| 7 | Do oxyanion holes and pockets separate the lipid-hydrolase panel? | Second donor varies; pocket size partly; no lining signature |
| 8 | How robust are the conclusions? | Two reviews withdrew or revised six claims |

## Slide list (M = method, R = result)
Title; Contents.
**Q1** divider; M merge sources and collapse genes; R source overlap; R attrition of 444 accessions; R data gaps; M+R structures trimmed (pLDDT).
**Q2** divider; M elbow rule and benchmark; R phi/psi plot; R recall vs motif; R robustness (36 settings, 47/49 PDB, Gly/Ser caveat);
M triad definition; R categories; R triad geometry (93 triads).
**Q3** divider; M superposition, SOAP and null; R RMSD distributions; R SOAP vs random serines; R pocket composition proxy.
**Q4** divider; M six distance spaces and nearest-neighbour test; R agreement between spaces; R name predicts neighbours; M controls;
R paralog removal; R length control and broader label; M clustering and random-set null; R cluster occupancy; R majority-ABHD
clusters by space; R not compact, FunFam purity, discordant pairs.
**Q5** divider; M leave-one-out EC transfer and family exclusion; R EC transfer; R hypotheses and gating; M anchored windows; R windows vs whole domain.
**Q6** divider; M alternative nucleophile scan and TM-align mapping; R retained His/acid; R no alternative nucleophile; M orthologue pipeline;
R orthologue residues by class; M shell-divergence test; R shell test (not significant); R ABHD4 insertion.
**Q7** divider; M oxyanion donors; R first donor; R second-donor residues by clade; M grid cavity and variants; R model dependence; R cavity ranking;
R robust contrasts; M pocket-lining mapping; R no lining signature.
**Q8** divider; R what the two reviews changed (withdrawn/revised/held); R four errors the reviews caught; Take-home; Limits and next steps.

## New figures (scripts/23_presentation_figures.py)
attrition, RMSD distributions, SOAP null, concordance heatmap, cluster-occupancy null, majority-ABHD clusters by space,
anchored windows vs whole domain, cavity model dependence, pocket-lining vs core identity.

## Rules
Claims no stronger than reports 00-12; "AlphaFold models" stated where relevant; no clade-level selectivity claim; the review slides stay.
## Delivered (v2)
Slides artifact: https://claude.ai/artifact/24tKMq7M5BujuB9swraudj (private). 63 slides: title, question, contents, eight sections
(divider, then method / result slides), review, take-home, limits. Source in `docs/presentation_deck/project/`; figures from
`scripts/23_presentation_figures.py` (figure `/_blob/<id>` references map to PNGs in `figures/presentation/`).
Layout was checked by arithmetic only (the preview pane became unavailable); page through once.
