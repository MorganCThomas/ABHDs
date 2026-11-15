# Stage 2: structures and domain trimming

Scripts: `scripts/02_fetch_structures.py --universe`, `scripts/08_fetch_annotations.py`, `scripts/09_structure_qc.py`.
Outputs: `data/derived/structure_qc.csv`, `data/derived/domains/*.pdb`, `attrition_stage2.csv`; experimental check in
`experimental_benchmark.csv` (script 16).

## What was done
AlphaFold models were fetched for the canonical entry of every universe gene (110/110 found; 59 new, 51 already held).
Each model was trimmed to the hydrolase domain: the span of the CATH ranges recorded for that protein in the FunFam
alignment headers (union over all FunFams containing it), or for the 5 proteins outside any FunFam (ABHD16B, ABHD18,
SERHL, SERHL2, EPHX4) the longest run of residues with pLDDT >= 70
(gaps of <= 15 bridged). These five are genes with no FunFam alignment entry (ABHD16B, ABHD18, SERHL, SERHL2, EPHX4).

## Results
- No gene was dropped: every domain has mean pLDDT >= 70 (median 94.6, minimum 70.3) and length >= 177 (median 299,
  maximum 685 residues).
- 42 of 110 proteins have a signal peptide and 27 a transmembrane helix (UniProt), which the domain trimming removes.
- **Experimental check.** 49 universe genes have PDB entries (up to 3 each tested). The AlphaFold elbow (same residue type
  and 5-residue context) is found in at least one experimental structure for 47 of 49 genes (92 of 109 entries tested;
  the others are fragments or constructs without the region). The two misses are informative: ABHD6 (7OTS) is a
  Ser148Ala variant - the elbow geometry is present but the residue is Ala, so the elbow is a backbone feature that does
  not depend on the nucleophile being there; for LIPE no elbow was found in the single entry tested (not investigated).

## Caveats
- Domain boundaries are CATH's, so they inherit CATH's choices (LCAT is split across three FunFams in the alignment files
  although the supplied metadata lists it in one; the union span is used).
- The 5 pLDDT-defined domains have no CATH range, so their boundaries are rougher.
