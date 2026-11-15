# ABHD superfamily: structure vs. name

**Question.** The α/β-hydrolase (ABH) superfamily is much larger than the genes named `ABHD*`. This repo maps that
(human UniProt entries, CATH FunFams, AlphaFold models) and asks whether sequence, FunFam and active-site
geometry agree on who is related to whom.

## Run

```bash
make env        # once: uv sync (Python 3.12, incl. torch for the optional ProtBert stage)
make            # data -> derived tables -> executed notebooks (only stale steps re-run)
make help       # list targets;  make clean-derived  to rebuild from raw
```

## Protein tools (case studies)

Everything the study computed about one protein, or the relationship between two, in one place (all from the derived tables; nothing is recomputed):

```bash
uv run python scripts/protein_report.py ABHD2               # HTML report + superposed 3D view -> out/protein_reports/
uv run python scripts/protein_report.py ABHD2 --orthologues # adds orthologue conservation (fetches UniProt the first time)
uv run python scripts/protein_report.py --all --no-structure # every gene + summary.csv
uv run python scripts/protein_compare.py ABHD2 MGLL         # similarity in every view, aligned catalytic and pocket residues
make report GENE=ABHD2                                      # same as the first line
```

The report covers identity and domain, catalytic site (elbow, triad, oxyanion donors, cavity, pocket residues), neighbours in five views with identity,
function and label quality (with the identity tiers found in report 13), position in the family, and a guardrails box. Examples: `docs/examples/`.
Python: `from abhd.atlas import Atlas` (data access), `abhd.protein_report`, `abhd.compare`, `abhd.orthologs`, `abhd.viewer`. Template notebook:
`notebooks/templates/protein_report.ipynb`; smoke test: `uv run python tests/smoke.py`. Genes outside the 110-gene universe are not supported yet.

**Publishing.** `scripts/publish_public.sh` rebuilds the local `public` branch from `main` for a free public GitHub repo (no raw data, no large files, redactions applied, safety checks); it never pushes: `git push origin public:main`.

**Data.** Raw inputs are not in the public repository; see [`data/README.md`](data/README.md) for how to obtain or regenerate them (UniProt, AlphaFold and PDB are fetched by the scripts; the CATH FunFam inputs must be supplied).

Raw downloads are never re-fetched by `make` (UniProt drifts; see *Gotchas*). Needs [uv](https://docs.astral.sh/uv/).

## Layout

| Path | What |
|---|---|
| `data/raw/` | Inputs, never edited: `uniprot/` (frozen `accessions.txt` + entries/FASTA), `structures/` (AlphaFold PDBs), `christina/` (CATH FunFam alignments + HHalign E-value matrix) |
| `data/derived/` | Generated tables: `protbert.npz`, `nucleophile_elbow.csv`, `soap.npz`, `*_cosine.csv`, … |
| `scripts/01…06_*.py` | One script per stage (see below) |
| `src/abhd/` | Shared code. **All paths and cutoffs live in `config.py`** |
| `notebooks/` | Interpretation only; they read `data/derived/` and are executed by `make` |
| `figures/` | Interactive Plotly HTML (hover for gene/accession) |
| `archive/`, `docs/`, `pymol/` | Original notebook and prototypes, slides, PyMOL session |

## Stages

| # | Script | Output |
|---|---|---|
| 1 | `01_fetch_uniprot.py` | entries + sequences for the 210 frozen accessions |
| 2 | `02_fetch_structures.py` | AlphaFold PDBs (203/210 exist) |
| 3 | `03_embed_protbert.py` | ProtBert mean embeddings (reuses existing `.pt`; model only loaded if some are missing) |
| 4 | `04_nucleophile_elbow.py` | `nucleophile_elbow.csv`: elbow call, legacy-regex call, His/Asp triad per candidate |
| 5 | `05_soap.py` | SOAP at the elbow nucleophile; ProtBert and SOAP cosine-distance tables |
| 6 | `06_funfam_names.py` | FunFam UniProt → gene map |
| 7 | `07_universe.py` | study stage 1: merged gene-level universe + attrition ledger (`make universe`) |
| 8 | `08_fetch_annotations.py` | UniProt active-site, EC, PDB annotations for the universe |
| 9 | `09_structure_qc.py` | stage 2: QC + hydrolase-domain trimming |
| 10 | `10_catalytic_sites.py` | stages 3-4: elbow, triad, benchmark, categories |
| 11 | `11_local_site.py` | stage 5: site RMSD, SOAP, null baseline |
| 12 | `12_global_comparison.py` | stage 6: TM-align, sequence, ProtBert, FunFam spaces |
| 13 | `13_abhd_label_analysis.py` | stage 6b: does the ABHD name mark a clade? |
| 14 | `14_function_inference.py` | stage 7: leave-one-out EC transfer |
| 15 | `15_anchored_alignment.py` | stage 8: alignment anchored on catalytic residues |
| 16 | `16_experimental_benchmark.py` | AlphaFold vs PDB elbow check |
| 17 | `17_review_controls.py` | Controls added after the independent review |
| 18 | `18_pocket_oxyanion.py` | Grid-cavity pocket descriptors and oxyanion donors |
| 19 | `19_inactive_vs_active.py` | ABHD15/ABHD5 vs active relatives: alt nucleophile, TM-align mapping, orthologues |
| 20 | `20_panel_pockets.py` | MGLL/ABHD6/ABHD2 panel: second donor, aligned pocket residues |
| 21 | `21_cavity_robustness.py` | Cavity volumes across model, anchor, orientation and clearance |
| 22 | `22_orthologues.py` | Cleaned vertebrate orthologue sets, conservation by taxon class |
| 24 | `24_label_provenance.py` | Evidence classes of EC / Rhea / GO labels |
| 25 | `25_function_decay.py` | Function transfer vs allowed neighbour identity, label-quality strata |
| 26 | `26_universe_sources.py` | Third universe source (InterPro/Pfam) and candidate genes |
| 27 | `27_core_domains.py` | Elbow-centred cores versus CATH domains; ABHD tests |
| 28 | `28_multiview.py` | Agreement between structure, sequence, FunFam and EC views |
| 29 | `29_nucleophile_measures.py` | Elbow vs motif vs structural transfer on the gold set |

`make study` runs stages 1-8 (needs network for the first fetches and `--extra embed` for ProtBert).

## Reading the results

1. `01_superfamily_overview`: only 17 of 107 human genes in the 76 FunFams are `ABHD*`; 15 of 76 FunFams contain one.
2. `02_sequence_embeddings`: ProtBert neighbours/UMAP. Distances are compressed; trust ranks, not values.
3. `03_nucleophile_elbow`: the elbow is the nucleophile at the β5→αC turn (φ≈+52°, ψ≈−120°, β before, helix after).
   It finds exactly one nucleophile in 160 of 203 models (the rest are mostly truncated isoforms) and also recovers
   cases the Gly-motif regex misses (ABHD16A, PLA2G15, Asp-nucleophile epoxide hydrolases).
4. `04_structural_similarity`: SOAP around the nucleophile vs ProtBert vs FunFam E-values. They agree on the closest
   paralog and disagree beyond it.

## Study progress

Design: `docs/study_design.md`. Per-stage observations and figures are in `reports/` (start with `00_summary.md`).

## Gotchas

- **UniProt drifts.** The family query returned 210 human entries in 2025 and 378 in Nov 2026; four of the 210 have
  since been deleted and one sequence (ABHD18, Q0P651) changed length. The study set is frozen in
  `data/raw/uniprot/accessions.txt`; `--from-query` re-runs the query (the Nov 2026 result is kept in
  `archive/uniprot_family_query_2026-11-11/`).
- "ABHD-named" = gene name starts with `ABHD` (36 of 210 entries now; the count was 37 in the original run).
- Residue numbers are UniProt/AlphaFold numbering (the slides' ABHD6 S147 is S148 here).
- Many entries are isoform/fragment models; analyses mix them with canonical entries (reviewed-only views are shown where it matters).
- ProtBert truncates at 510 residues; SOAP/elbow depend on AlphaFold model quality (pLDDT ≥ 70 required at the nucleophile).
- `christina/` is read-only input; its paths inside `funfam_metadata.csv` refer to her machine.
