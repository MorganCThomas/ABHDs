# Stage 4: catalytic categories

Script: `scripts/10_catalytic_sites.py` (category column of `data/derived/catalytic_sites.csv`); figure
`figures/stage4_categories.html`.

## Result (110 genes)
| Category | Genes | ABHD-named | Members |
|---|---|---|---|
| Ser-His-Asp | 77 | 18 | 18 of the 20 ABHD genes plus 59 others |
| Ser-His-Glu | 8 | 0 | CES1, CES2, CES3, CES4A, CES5A, ACHE, BCHE, SERAC1 |
| Asp at elbow, His, Asp | 5 | 0 | EPHX2, EPHX3, EPHX4, MEST, DPP6 |
| Asp at elbow, His, Glu | 1 | 0 | EPHX1 |
| Cys-His-Asp | 2 | 0 | BAAT, CMBL |
| Ser at elbow, no His in range | 3 | 0 | SPG21, SERHL, CES1P1 |
| No nucleophile, Gly at elbow | 12 | 1 | ABHD15, ACOT6, NLGN1-4Y (5), DPP10, NDRG1-4 |
| No nucleophile, Asn at elbow | 2 | 1 | ABHD5, KANSL3 |

## Interpretation
- The elbow backbone geometry is present in all 110 genes (a strained phi/psi turn after a strand), including the 14 where
  the residue is not a nucleophile type: the fold keeps the turn when the chemistry is lost. ABHD5 Asn153 (a known
  acyltransferase co-activator, not a hydrolase) and the neuroligins, NDRGs, DPP10 and KANSL3 are the pseudoenzyme-like cases.
- Two ABHD-named genes break the "ABHD = hydrolase" expectation: ABHD5 (Asn) and ABHD15 (Gly at the elbow). Conversely 59
  non-ABHD genes carry the same Ser-His-Asp site, so the name does not identify a catalytic type.
- Ser-His-Glu in cholinesterases and carboxylesterases, Cys-His-Asp in BAAT, and Asp-based sites in the epoxide hydrolases
  match what is known about these families (stated from background knowledge, not tested here).
- The residue labelled "Asp at elbow" is not proof of function: DPP6 is, to my knowledge, catalytically inactive despite
  this geometry (not verified here), so categories describe the site, not activity.

## Caveats
- Without orthologue conservation "no nucleophile" means "no nucleophile-type residue at the elbow in the AlphaFold model". "No His in range" (SPG21, SERHL, CES1P1) may be a model artefact or a real loss; no orthologue conservation check was done.
- Gly "nucleophile" calls come from the strained phi/psi at a glycine. They are weaker evidence than Ser/Asp calls; for the
  neuroligins the Ser sits one residue before the Gly in the sequence, but the Gly aligns to the GESAG nucleophile position of the carboxylesterase family, so a register shift is not the explanation (ABHD15, SEGSG with the Ser at i+1, may be shifted).
