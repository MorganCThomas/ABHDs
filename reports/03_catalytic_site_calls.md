# Stage 3: catalytic-site calls and benchmark

Scripts: `scripts/10_catalytic_sites.py`, `src/abhd/catalytic.py`. Outputs: `data/derived/catalytic_sites.csv`,
`elbow_table.csv`, `benchmark_stage3*.csv`, `sensitivity_stage3.csv`; figure `figures/stage3_elbow_phi_psi.html`.

## What was done
Per gene: the nucleophile elbow at any residue type (phi 30-90, psi -170 to -80, beta before, helix after, pLDDT >= 70),
the His / Asp-Glu pair with heavy-atom distances <= 4.5 / 4.0 A and the chain order nucleophile < acid < His, and the
nearest non-local backbone N as an oxyanion-hole candidate. Calls were scored against UniProt "Active site" features and
against the original notebook regex.

## Results
- **Benchmark.** 31 genes carry a UniProt nucleophile annotation. The elbow recovers all 31 (31/31); the legacy regex
  recovers 26/31 (it misses the Asp nucleophiles of EPHX1-4 and PPT2 Ser111).
- **Specificity.** 95 of the 96 genes with an S/C/T/D elbow have exactly one such residue (PREP has two: Ser554, Thr590).
  The regex flags a median of 3 serine/threonine/cysteine residues per domain.
- **Robustness.** 36 threshold combinations (phi, psi, pLDDT 60/70/80, with and without the strand/helix condition) all
  recover 31/31. Without the strand/helix condition 4-5 genes get more than one S/C/T/D call.
- **Triad.** 93 genes have a complete triad within the cutoffs; all 93 obey the chain order nucleophile < acid < His.
  Distances: nucleophile-His 2.74-3.88 A (mean 2.98), His-acid 2.61-3.78 A (mean 2.79). The His agrees with the
  annotated His in 28 of 28 genes where both exist; the full His + acid pair agrees in 24 of 28. The four differences are
  EPHX1-4, where UniProt's second catalytic position is not the Asp/Glu found geometrically (e.g. EPHX2 annotated 466,
  found 496).
- **Oxyanion hole.** The nearest non-local backbone N is 2.8-5.95 A from the nucleophile (median 4.7 A) and varies in
  identity; this is descriptive only and was not validated against known oxyanion-hole residues.

## Caveats
- Only 31 of 110 genes have a UniProt nucleophile annotation, so recall is measured on a modest, curated subset; the
  annotated set is biased toward well-studied enzymes.
- Precision beyond "one call per gene" cannot be measured without negatives; the experimental-structure check (stage 2
  report) covers geometry reproducibility, not correctness of the residue assignment.
