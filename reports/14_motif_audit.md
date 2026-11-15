# Stage 14: motif audit (HX4D and the nucleophile elbow)

Script `scripts/30_motif_audit.py` (`make data/derived/motif_audit_summary.csv`). Outputs `data/derived/motif_audit_{triads,offsets,summary,orthologues}.csv`. Date 2026-11-15.
Question: do the historical motifs (HX4D "acyltransferase"; small-X-Nu-X-small / GXSXG elbow) mean in our models what they were taken to mean?
Method: for the 93 complete triads, read the residues around the catalytic His and the nucleophile in the full-length AlphaFold sequence.

## HX4D (Asp at His+5)
- 12 of 93 triad histidines start HX4D (13%); other histidines of the same proteins 5.3%. One-sided permutation test (one random non-catalytic His per protein,
  20,000 draws): p = 0.004. This treats proteins as independent; the 12 hits lie in 10 of 64 FunFams, so the evidence is weaker than the p-value.
- The enrichment is not specific to His+5: Asp is 17% at His-1, 14% at His-3, 16% at His+2, 12% at His+8 and about 1% at His+4 and His+6. The His loop is acid-rich with a
  period-3 pattern. Asp or Glu at His+5: 28 of 93 (30%) against 12% at other histidines; 16 triads have HX4E.
- The motif Asp is never the triad acid (0 of 93; the acid precedes the His in chain order). The motif is the His loop, not a second site.
- Hits: ABHD11, ABHD4, ABHD14B, PPME1, PLA2G7, PAFAH2, EPHX2, EPHX4, OLAH, CTSA, CPVL, SCPEP1: mostly hydrolases and serine carboxypeptidases.
- EC 2.3 (acyltransferase) records: 6 proteins, 2 with HX4D (ABHD14B, PAFAH2, which also has a hydrolase EC); versus 10 of 87 others (Fisher p = 0.17; underpowered, no evidence either way).
- ABHD-named: 3 of 18 triads (ABHD11, ABHD4, ABHD14B). ABHD5 has no triad call; its His (H327, from the TM-align mapping to ABHD4 H320) starts HYVYAD, so it keeps HX4D; ABHD15 (H391, HCGFLR) does not.
- Orthologue sets (cleaned vertebrate): the +5 residue is Asp in ABHD4 (D325, 100%) and ABHD5 (D332, 100%), Glu in ABHD1/2/3 (>= 99.8%), Arg in ABHD15 (84%). The motif follows relatedness (ABHD4/ABHD5),
  and ABHD4 is a hydrolase.
- Reading: HX4D is a mildly enriched, non-specific feature of the His loop. It was inferred from ABHD5 and does not mark acyltransferase activity in this set. Whether ABHD5 is an acyltransferase is untested here.

## Elbow motif
- small (G, A, S, C, T) at Nu-2 and Nu+2 holds in 93 of 93 triads (12% at random Ser); the elbow was called from backbone geometry, so this is independent of the sequence pattern.
- Strict GXSXG (G.[SC].G): 75 of 93, and 15 of 18 ABHD triads. Exceptions among ABHDs: ABHD14A/B (SPSLS) and ABHD16A (AWSIG).

Limits: sequence read from AlphaFold models; the +5 offset is a fixed-spacing definition of HX4D; no experimental mutants of the motif exist in the cited papers.
