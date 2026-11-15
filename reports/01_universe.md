# Stage 1: merged universe

Script: `scripts/07_universe.py` (`make universe`). Outputs: `data/derived/universe.csv` (one row per gene),
`universe_entries.csv` (every source accession and its fate), `attrition_stage1.csv`;
figures `figures/stage1_universe_overlap.html`, `figures/stage1_attrition.html`.
Design: `docs/study_design.md`, stage 1.

## What was done
Three sources were merged and collapsed to one canonical (reviewed) entry per gene:
the frozen UniProt family query (210 entries, 2025), the same query re-run 2026-11-11 (378 entries), and the human
members of the CATH 3.40.50.1820 FunFams (103 accessions in 76 FunFams). Unreviewed entries were placed on a gene by
name, or by local alignment (>=90% identity over >=80% of the fragment). Everything else is dropped with a reason code.

## Result
- **Universe: 110 genes, 20 ABHD-named (18%).** 44 are in both UniProt and FunFam, 59 FunFam only, 7 UniProt only.
- **Attrition (444 distinct source accessions):** 110 canonical; 274 alternative entries matched to those genes by
  name (isoforms, fragments); 54 more placed on a gene by alignment; 5 unassignable fragments; 1 left over
  (`ABHD14A-ACY1` readthrough).
- **The UniProt family query is not the superfamily.** Its reviewed hits are 51 genes in both 2025 and 2026. The 59
  FunFam-only genes include ABHD4, ABHD5, ABHD12 and ABHD13, plus carboxylesterases, DPP4/FAP/PREP, neuroligins,
  NDRG1-4, ACHE/BCHE, EPHX1, ACOT1/2/4/6, FASN and others.
  Likely reason (family line of the 59, from UniProt "Sequence similarities"): they are annotated under other named
  families (type-B carboxylesterase, peptidases S9/S10/S28/S33, NDRG, thioesterase) rather than "AB hydrolase
  superfamily". I read the family text for the FunFam-only set; I did not test the query logic itself.
- **Three ABHD-named genes are not in the supplied human FunFam list:** ABHD16A, ABHD16B and ABHD18. Checked against the
  1,873 alignment files in `christina/FF_alns`: ABHD16A is in FF-000074 and DAGLB in FF-000064, but neither FunFam
  appears in `funfam_metadata.csv` (so Christina's list of 76 FunFams with human members is incomplete); ABHD16B,
  ABHD18, SERHL, SERHL2 and EPHX4 are in no FunFam alignment at all. Stage 2 uses the alignment headers for domain
  boundaries, but stage 6 FunFam distances (HHalign matrix) still lack ABHD16A and DAGLB.
- **Name vs family:** in the FunFam set 86 of 103 human proteins are not ABHD-named; only 17 are.
  (README previously said 17 of 107 genes; this count uses accessions mapped to canonical genes.)

## Caveats / open decisions
- AlphaFold models exist for only 203 of the original 210 entries, and for none of the 59 FunFam-only genes yet (stage 2).
- Pseudogene `CES1P1` and non-hydrolase FunFam members (FASN, KANSL3, NDRG) are in the universe: keep them for the
  naming question, but stage 3 should expect them to lack a canonical triad.
- Alias gene names (e.g. BAT5 = ABHD16A) resolved by sequence; not checked by hand beyond the ledger.
- Reviewed-only canonical choice means genes with only unreviewed entries (none found beyond `ABHD14A-ACY1`) would be lost.
