# Stage 5: local comparison of catalytic sites

Script: `scripts/11_local_site.py`. Outputs: `data/derived/site_rmsd_elbow.csv`, `site_rmsd_triad.csv`,
`soap_site_cosine.csv`, `soap_null_summary.csv`, `site_env.csv`; figure `figures/stage5_soap_null.html`.

## Results
- **Elbow backbone** (N, CA, C of residues i-2..i+2, 110 genes): median pairwise RMSD 0.22 A. This is by construction: the residues were selected
  on phi/psi and neighbouring secondary structure, so it shows consistency of the selection, not an independent finding.
- **Triad skeleton** (CA and functional group of nucleophile, His, acid; 89 genes whose triad atoms all lie in the trimmed
  domain): median pairwise RMSD 0.58 A.
- **SOAP at the elbow CA** (8 A cutoff): median cosine distance 0.0069 between elbow sites, against 0.034 between an elbow
  site and a random serine and 0.042 between random serines. Elbow sites are about five times closer to each other than to
  random sites, but within-elbow spread (5th-95th percentile 0.002-0.034) overlaps the lower tail of the random baseline
  (5th percentile 0.008), so SOAP alone discriminates elbow from non-elbow only coarsely.
- **Pocket proxy** (residues within 10 A of the nucleophile): median 36 residues for ABHD-named vs 41 others, hydrophobic
  fraction 0.46 vs 0.45 (not different, p=0.45). Residues within 10 A do differ (36 vs 41, p=0.001) and track domain length (rho=0.49), so the proxy is largely size.

## Interpretation
The catalytic core is a conserved module; differences between family members are not in the elbow itself. Local similarity
therefore reports shared chemistry, not shared substrate preference (confirmed by stage 7, where local-site distances
transfer EC no better than a majority baseline once paralog neighbours are excluded; EC level 3 is itself a chemistry label, so substrate-level specificity was not tested). The SOAP null uses random Ser CA atoms, which is not a matched null (other strand-to-helix turns would be).

## Caveats
- The pocket proxy is a crude composition count, not a cavity calculation (no pocket-detection tool was run).
- SOAP is centred on CA so genes without a nucleophile can be included; this reduces its sensitivity to side-chain chemistry.
