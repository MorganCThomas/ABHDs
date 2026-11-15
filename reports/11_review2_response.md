# Second independent review: response (2026-11-13)

A second reviewer re-ran the revised pipeline and the follow-ups. It reproduced all headline numbers in reports 00 and 09 and most of
reports/10 (EC leave-family-out, seqid fix, abhd_family, A1/A2/A5, volumes, Asn counts), and found the following.

## Accepted and fixed
| Finding | Action |
|---|---|
| Near-neighbour exclusion used TM < 0.2 in every space (removed 202 TM pairs but 5 sequence pairs), so "sequence stays 0.80" was vacuous | `scripts/17` now drops the same share (3.4%) of closest pairs in each space. Nearest-neighbour rate: TM 0.45, shorter-chain TM 0.60, sequence 0.60, FunFam 0.71, elbow window 0.55. Reports 00, 06, 09 and the design doc reworded: the signal is paralog-driven in sequence too. |
| ABHD5/ABHD4 shell test treats residues as independent | `scripts/19` adds a circular-shift (run-aware) p and both reference directions: ABHD5 vs ABHD4 p=0.14 (reference ABHD4), 0.44 (reference ABHD5). No ABHD15 contrast is significant, and the active-active ABHD1/ABHD2 pair gives the same odds ratio of 0. Report 10 reworded; "pocket remodelled" withdrawn. |
| Cavity volumes not robust to model (trimmed vs full-length), clearance and CA vs OG anchor | Reviewer's recomputation accepted (ABHD2 at 3.6 A: 64 -> 148, above MGLL/ABHD6). All panel pocket-size ordering and the ABHD2 selectivity paragraph withdrawn. Redo listed as next step 0. |
| Report 10 ranges did not match the CSVs (Nu+1 N 2.77-3.73, second donor 2.8-5.4, 11 not 10 genes without donor, RMSD to 4.26, 351/359 not 351/363) | Corrected. |
| Asn enrichment overstated (4 stems; 14 residue types tested; non-random missing donors) | Reworded as suggestive; stem-collapsed 4/14 vs 3/53 (p=0.03). |
| Orthologue contaminants and unpaginated query | Acknowledged in report 10; fix listed as next step. |
| "vertebrate-wide", "late loss", ligand-binding speculation, ABHD2 inhibitor-selectivity | Removed or reworded. |
| 00_summary inconsistencies (pocket-descriptor line; completed next steps) | Corrected. |

## Accepted, not yet done (added to docs/next_steps.md)
- Cavities with full-length models, truncated side chains, rotated grids, larger box.
- Run-aware shell tests with all active-active pairs in both directions; sensitivity to the 3 A and 6 A thresholds.
- Orthologue hygiene: pagination, best-hit against all paralogs, best identity per species, class-level counts and a small tree for ABHD15 Gly vs Ser-Thr.
- Oxyanion test on stems with FDR, S9 peptidases handled explicitly; positive control for the alternative-nucleophile scan.
- Leave-family-out by FunFam plus sequence-identity threshold (LYPLAL1 and LYPLA1/2 currently treated as different families) and paired tests against the 0.59 baseline.
- Engineering: `make` does not regenerate git-ignored structures/domains/cavity files on a clean clone and `$(FOLLOW)` has incomplete dependencies; script 20 imports script 19 by path (move shared code to `src/abhd`); hardcoded residue numbers in script 19; no tests; rerun `make study` end to end.

## Net effect
No headline of reports 00/09 changes beyond the exclusion control. The follow-up conclusions shrink to: both inactive proteins lack
only the nucleophile in the models and the Gly/Asn is common across sampled vertebrates; localised ABHD4 180-198 displacement is
unconfirmed; pocket-size and clade-selectivity claims are not supported.
