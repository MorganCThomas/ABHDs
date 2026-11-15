# Redesign plan: CATH, FunFam, EC and UniProt annotations as measurements, not ground truth

Status: PLAN v2 (revised after independent review; v1 is commit 490da0e). Date: 2026-11-15.

## Principle
Every label we used as truth is itself a measurement with its own error and its own dependence on the others:
- **CATH/FunFam**: built from sequence (HMM) similarity; defines the universe (59 of 110 genes only come from it) and the domain boundaries.
- **EC numbers**: partly assigned by similarity (inferred), many incomplete; EC level 3 is a chemistry class, not a substrate.
- **UniProt active-site annotations** (the nucleophile gold set): 7 of 31 experimental, 15 propagated by similarity.
- **ABHD name**: a nomenclature history, also a measurement.
The redesign (a) replaces "accuracy against truth" with agreement between views where the label was the reference, (b) removes
single-definition dependence (family exclusion, domain trimming, universe), and (c) tracks label provenance.
Analyses that do not use these labels (elbow geometry, local sites, cavities, orthologues) are unchanged.

## Planned changes (each is a new script plus outputs; no existing script is modified except where stated)

### P1. Label provenance (scripts/24_label_provenance.py)
- Inputs: `data/raw/uniprot/universe_annotations.tsv` (EC, Catalytic activity with Evidence=ECO codes and Rhea IDs, GO MF, PDB).
- Parse per gene: EC numbers with evidence class (experimental = any ECO:0000269 / ECO:0000314 / ECO:0000315 / ECO:0000305-curated not counted;
  inferred = ECO:0000250 (similarity), ECO:0000255 (rule/automatic), ECO:0000256/0000259/0000501 (automatic), ECO:0000305 (curator inference)),
  Rhea reaction IDs, GO MF terms, count of PDB entries.
- Output `data/derived/label_provenance.csv` (one row per gene) with: ec_full_levels, ec_experimental (bool), n_rhea, n_go_mf, n_pdb.
- Report: how many EC labels are complete, experimental, inferred. Expected: a minority experimental.

### P2. Function transfer as a decay curve, stratified by label quality (scripts/25_function_decay.py)
- Compute true pairwise percent identity (BLOSUM-free, global alignment identity over aligned columns, on the trimmed domain) for all gene pairs -> `data/derived/identity_pairs.csv`.
- For each transfer space (tmalign, tmshort, seqid, funfam, anchor_elbow, anchor_cat, soap, triad, protbert) and each maximum allowed identity cutoff
  c in {100, 80, 60, 50, 40, 30, 25%}: leave-one-out EC3 nearest-neighbour accuracy using only neighbours with identity <= c; Wilson CI; n queries; two
  baselines (majority class, random neighbour) recomputed on the same pool. Output `function_decay.csv` and `figures/function_decay.html`.
- Stratification: (i) all complete EC3; (ii) query's EC is experimentally evidenced; (iii) query and neighbour both experimental. Output `function_decay_strata.csv`.
- Existing leave-family-out result is kept as one point on the curve and cross-checked.
- Decision rule: report the identity at which the best space's accuracy CI first overlaps the majority baseline.

### P3. Third universe source and a lower bound on family size (scripts/26_universe_sources.py)
- Query UniProt (reviewed human) for xref InterPro IPR029058 (alpha/beta hydrolase fold superfamily) and Pfam clan members via `xref:pfam-PF00561,PF12697,PF12146,PF00135,PF00326,PF01738`; write `data/raw/uniprot/third_source.tsv`.
- Compute overlaps among the three lists (UniProt 2025/2026 family query, CATH FunFam, InterPro/Pfam) with gene-level collapsing identical to stage 1; list genes in the third source but not in the 110 (candidate additions, not auto-added).
- Estimate size: Chao lower bound from the three-list capture pattern, stated with the independence caveat (lists share sequence evidence).
- Sensitivity: rerun the ABHD nearest-neighbour test (script 17 labels) after adding those candidates whose AlphaFold model passes QC IF they number <= 30; otherwise report counts only.
- Output `universe_sources.csv`, `universe_candidates.csv`.

### P4. Structure-derived, size-controlled cores instead of CATH ranges (scripts/27_core_domains.py)
- For each of 110 genes take the K=180 residues whose C-alpha are nearest to the elbow C-alpha (sequence order retained, fragments allowed); K chosen as the 10th percentile of CATH domain lengths (check, adjust).
- TM-align all-vs-all on cores (`tmtools`, same normalisations) -> `dist_tmalign_core.csv`, `dist_tmshort_core.csv`; local-sequence distance on core residues -> `dist_seqid_core.csv`.
- Rerun the ABHD nearest-neighbour test, R-precision, length-stratified permutation (now unnecessary but reported) and the quantile-matched close-pair exclusion on core distances; compare with CATH-domain results in `core_vs_cath_label_tests.csv`.
- Decision rule: the CATH-dependence is a problem if the nearest-neighbour rate or R-precision changes by more than the bootstrap CI half-width.

### P5. Multi-view agreement framework (scripts/28_multiview.py)
- Views, each as a gene x gene distance matrix on the same 110 genes (NaN where undefined): structure (TM), structure short-norm, local sequence, FunFam HMM, protein-language-model, SOAP, triad RMSD,
  elbow window, catalytic windows, core structure (P4), EC (1 - depth of shared EC prefix / 4, over complete levels; undefined if no complete EC), GO MF (Jaccard distance on terms),
  Rhea (Jaccard on reaction sets), catalytic category (0/1), ABHD label (0 within the labelled group, 1 across; as an indicator view), name stem family (0 same stem / 1).
- Association between every pair of views: Spearman correlation of upper-triangle distances (on genes defined in both) plus distance correlation; permutation p (5000) by permuting gene labels of one matrix (Mantel); gene-block bootstrap CI (500) to account for dependence of pairs sharing a gene.
- Independence audit: report which pairs are dependent by construction (FunFam-sequence; EC-FunFam; Rhea-EC; GO-EC) and mark them in the figure.
- Consensus: mean of rank-normalised matrices from structure, sequence, FunFam, windows; for each view the correlation with the leave-one-out consensus (a measure of view reliability); for each gene the discordance = 1 - Spearman between its row in structure space and its row in function views (EC, GO, Rhea); list the 15 most discordant genes.
- Outputs: `multiview_association.csv`, `multiview_view_reliability.csv`, `multiview_gene_discordance.csv`, `figures/multiview_heatmap.png` plus an HTML copy.

### P6. Nucleophile calls as three measures (scripts/29_nucleophile_measures.py)
- Measures on all 110 genes for "a nucleophile-type residue (S/C/T/D) exists at the elbow position": (a) elbow call; (b) legacy sequence motif matches at that residue; (c) structural transfer: the residue aligned (TM-align on cores) to the nucleophile position of the nearest of the 31 annotated genes.
- Agreement table (Cohen kappa pairwise) and a Dawid-Skene latent-class estimate of each measure's sensitivity and specificity (EM, 3 raters, binary; reported only if the EM converges and the three measures are not perfectly collinear).
- Gold subset (31 UniProt nucleophiles): evaluate each measure; stratify by evidence (experimental vs inferred).
- Output `nucleophile_measures.csv`, `nucleophile_latent.csv`.

### P7. Reports and plumbing
- `reports/13_labels_as_measures.md`: methods, results, what changed in earlier conclusions (stated explicitly), caveats.
- `Makefile`: targets for scripts 24-29 (dependencies on prior outputs); `README.md` stage table rows 24-29; `docs/study_design.md` execution note; `docs/next_steps.md`.
- Deck: not edited in this pass; a list of slides that would change is written to the report.

## Explicit non-goals / risks
- No change to the universe (110) in this pass beyond reporting candidates.
- Percent identity from pairwise global alignment is a proxy for "family" similarity; thresholds are exploratory.
- GO, EC and Rhea are not independent of each other; the multi-view association is not evidence of independent confirmation.
- Latent-class estimation with 3 raters is weakly identified; treat as descriptive.
- Network access to UniProt is required for P1 (already cached) and P3.

## Acceptance criteria
Each script runs from the Makefile on cached data (P3 fetches once and caches), writes its tables and a log line with counts; report numbers are read from tables; existing outputs are not overwritten.


---
# v2 amendments (supersede the conflicting parts of P1-P6 above)

Independent review findings accepted: ECO parsing was mis-specified; identity cutoffs above 40% are vacuous; the "best space" decision rule is a
winner's-curse; Chao and Dawid-Skene are not defensible; Mantel p-values are trivially small with shared paralog structure; K=180 does not match the
domain length distribution (5th percentile 210, 10th 221).

## Pre-registered conclusions each analysis could overturn
| Analysis | Conclusion in reports 00-12 it tests | Overturned if |
|---|---|---|
| P2 decay | "Across families, neighbours beat guessing by about ten points" (EC3 0.60-0.70 vs 0.59) | at identity <= 30% the paired difference to the majority baseline has a stem-cluster bootstrap CI that includes 0 for all prespecified spaces (supports) or excludes 0 (overturns the "weak" label) |
| P2 strata | EC transfer is not driven by inferred labels | experimental-only accuracy differs from all-label accuracy by more than the CI |
| P4 cores | ABHD nearest-neighbour signal is not a length/CATH-trimming artefact | nearest-neighbour rate on cores differs from CATH-domain rate by more than the paired-bootstrap equivalence margin (0.10) |
| P5 agreement | structure and sequence agree moderately; FunFam is a sequence-derived view | cross-family structure-sequence correlation is near 0 (overturns "agree moderately" within families only) |
| P3 sources | 110 is the family | a third source adds >= 10 structurally valid genes |
| P6 table | elbow beats the motif on the gold set | experimental-only subset reverses the ordering |

## P1 (revised): scripts/24_label_provenance.py
- Parse each `Catalytic activity` reaction block: EC (`EC=`), Rhea ID(s), Evidence ECO set. Per reaction class: experimental if ECO:0000269 present;
  else curator (ECO:0000305) ; else similarity (ECO:0000250); else rule (ECO:0000255); else none. (Codes present in the file: 269, 305, 250, 255; active-site also 303, 7744.)
- Per gene: n_reactions, n_rhea, n_go_mf, n_pdb, ec_complete_levels (0-4), ec3_set, evidence_any_experimental, evidence_all_experimental, evidence_class (experimental-only / mixed / non-experimental / no_record).
- GO evidence codes are not available in the download; GO is reported but not used for transfer.
- Acceptance: genes with catalytic blocks = 76; genes with Rhea = 68; complete EC3 genes = 80 and identical to the set in function_loo.csv; some-experimental ~58 and all-experimental ~39 among them (explain any difference).

## P2 (revised): scripts/25_function_decay.py
- Identity: global Needleman-Wunsch (BLOSUM62, gap -11/-1) on trimmed domains; identity = identical columns / aligned (non-gap) columns; also coverage; save `identity_pairs.csv`. (Local SW aligned-column identity is not used.)
- Cutoff grid on maximum neighbour identity: {100, 40, 35, 30, 25, 20}% and, alternatively, bins of the identity of each query's actual nearest labelled neighbour; log the pool size per cutoff.
- Prespecified spaces only: tmalign, seqid, funfam, anchor_elbow (others in a supplementary table). Neighbours must satisfy the identity cutoff AND, in the family-exclusion mode, the FunFam/stem exclusion of script 14.
- Statistics: per space, difference to the majority-class baseline (computed on the same queries) with a stem-cluster bootstrap (resample paralog stems, 2000 draws, drop duplicate genes); no best-space selection; nested pools are dependent and reported as a curve, not independent tests. McNemar on the 100% point vs script 14.
- Strata (descriptive): all complete-EC3; evidence_any_experimental; query and neighbour both all-experimental.
- Acceptance: c=100 reproduces function_loo.csv "all" rows (tmalign 0.838, n=80, baseline 0.588); family-exclusion at c=100 reproduces leave_family_out (tmalign 0.675); pools monotone in c; identity matrix symmetric with unit diagonal.

## P3 (revised): scripts/26_universe_sources.py
- Third source: UniProt reviewed human with xref InterPro IPR029058 (alpha/beta hydrolase fold superfamily) and, secondarily, Pfam families PF00561, PF12697, PF12146, PF00135, PF00326, PF01738, PF07859, PF02230, PF03096, PF06441, PF02450, PF00975 (record the fetch date; cache).
- Report a three-list overlap table (collapsed to genes as in stage 1) and the candidate list (third source not in the 110), with annotation columns. No Chao estimate. No rerun of analyses (the candidates are listed for structure QC in follow-up work); if >= 10 candidates, flag in the report.
- Acceptance: the first two lists reproduce 44 / 59 / 7 and the 110 genes; the overlap table sums to the union.

## P4 (revised): scripts/27_core_domains.py
- Core = the K residues whose C-alpha are nearest to the elbow C-alpha (3D), original order retained, for K in {150, 180, 200}; primary K=200 (about 0.9 x the 5th percentile of CATH domain length, 210); genes with fewer than K residues in the trimmed model use the whole model and are listed.
- Record per gene: fraction of core residues outside the CATH range, number of chain breaks. TM-align on cores (tmalign normalisation and shorter-chain), local-sequence distance on core residues.
- ABHD tests on cores: nearest-neighbour rate, R-precision, quantile-matched close-pair exclusion; paired bootstrap over ABHD genes of (core - CATH) differences with an equivalence margin of 0.10.
- The core is centred on the elbow, so it is NOT an independent check of the elbow measures; state this. PREP/PREPL cores include propeller residues.
- Acceptance: core TM distance correlates with dist_tmalign at Spearman > 0.7; CATH-domain ABHD nearest-neighbour rate reproduces 0.80 (tmalign, n=20); K=200 core length constant.

## P5 (revised, reduced): scripts/28_multiview.py
- Views (5): structure (TM-align, CATH domain), structure (core, K=200), local sequence, FunFam HMM, EC (complete EC, hierarchical distance 1 - shared prefix depth/4, defined for genes with complete EC3).
- For each pair of views: Spearman on upper-triangle pairs defined in both (n pairs), for ALL pairs and for CROSS-FAMILY pairs only (different paralog stem and no shared FunFam); stem-cluster bootstrap CI (resample stems, dedupe genes, 1000 draws). No Mantel permutation p.
- Dependence audit table: which pairs are dependent by construction.
- Per-gene discordance: EC3 of the nearest cross-family neighbour in structure, sequence and FunFam space; list genes where structure agrees with EC3 but sequence does not, and the reverse. No leave-one-out "consensus reliability".
- Acceptance: n pairs and NaN counts logged; the all-pair structure-sequence Spearman reproduces global_concordance.csv (0.52 for tmalign-seqid).

## P6 (revised, reduced): scripts/29_nucleophile_measures.py
- Descriptive per-residue table on the 31 gold genes only: elbow call, regex hit at the gold residue, structural-transfer call (leave-one-out: residue aligned to the gold nucleophile of the nearest other gold gene by TM-align on cores); by evidence of the gold annotation (experimental = ECO:0000269 in the Active site evidence; else inferred). Plus per-residue regex hit counts for the 14 no-nucleophile genes. No latent-class model.
- Acceptance: reproduces 31/31 (elbow) and 26/31 (regex); transfer is leave-one-out.

## Plumbing
Makefile targets for 24-29; README rows; `docs/study_design.md` note; `reports/13_labels_as_measures.md` stating results against the pre-registered table above.
