# ABHD pipeline. `make` rebuilds only what is stale. Needs uv (https://docs.astral.sh/uv/).
PY  := uv run python
PYE := uv run --extra embed python
SRC := $(wildcard src/abhd/*.py)

ENTRIES   := data/raw/uniprot/abhd_entries.json
FASTA     := data/raw/uniprot/abhd_sequences.fasta
STRUCTS   := data/raw/structures/fetch_log.txt
PROTBERT  := data/derived/protbert.npz
ELBOW     := data/derived/nucleophile_elbow.csv
SUMMARY   := data/derived/structure_summary.csv
SOAP      := data/derived/soap.npz
UNIVERSE  := data/derived/universe.csv
FFMAP     := data/derived/ff_uniprot_map.json
NBS       := $(patsubst notebooks/%.ipynb,notebooks/.done/%,$(wildcard notebooks/*.ipynb))

.PHONY: all data derived notebooks env clean-derived help universe study
all: notebooks                        ## everything: data -> derived tables -> executed notebooks
help:                                 ## list targets
	@grep -E '^[a-z-]+:.*##' Makefile | sed 's/:.*##/ -/'
env:                                  ## create .venv (add --extra embed yourself to recompute ProtBert)
	uv sync --extra embed
data: $(ENTRIES) $(STRUCTS) $(FFMAP) ## stages 1,2,6: network downloads (skipped if files exist)
derived: $(PROTBERT) $(ELBOW) $(SUMMARY) $(SOAP) ## stages 3-5: embeddings, elbow/triads, SOAP + distance tables
notebooks: $(NBS)                     ## execute notebooks in place

# Raw downloads are only fetched when missing (order-only prerequisites): editing code never re-downloads.
# To refresh on purpose: delete the file and `make`, or run scripts/01_fetch_uniprot.py --from-query.
$(ENTRIES) $(FASTA): | data/raw/uniprot/accessions.txt
	$(PY) scripts/01_fetch_uniprot.py
$(STRUCTS): | $(ENTRIES)
	$(PY) scripts/02_fetch_structures.py
$(PROTBERT): $(FASTA) scripts/03_embed_protbert.py src/abhd/config.py
	$(PYE) scripts/03_embed_protbert.py
$(ELBOW): $(STRUCTS) scripts/04_nucleophile_elbow.py $(SRC)
	$(PY) scripts/04_nucleophile_elbow.py
$(SUMMARY): $(ELBOW)   # written by the same script as $(ELBOW)
$(SOAP): $(ELBOW) $(PROTBERT) scripts/05_soap.py $(SRC)
	$(PY) scripts/05_soap.py
universe: $(UNIVERSE)                 ## study stage 1: merged UniProt + FunFam universe, attrition ledger, figures
$(UNIVERSE): $(ENTRIES) $(FFMAP) scripts/07_universe.py $(SRC)
	$(PY) scripts/07_universe.py
$(FFMAP): data/raw/christina/funfam_metadata.csv scripts/06_funfam_names.py
	$(PY) scripts/06_funfam_names.py

# ---- study stages 2-8 (see docs/study_design.md, reports/) ------------------------------------------------------
UFETCH   := data/raw/structures/fetch_log_universe.txt
UANN     := data/raw/uniprot/universe_annotations.tsv
QC       := data/derived/structure_qc.csv
SITES    := data/derived/catalytic_sites.csv
LOCAL    := data/derived/soap_site_cosine.csv
GLOBAL   := data/derived/dist_tmalign.csv
LABEL    := data/derived/abhd_label_summary.csv
FUNC     := data/derived/function_loo.csv
ANCHOR   := data/derived/dist_anchor_elbow.csv
REVIEW   := data/derived/review_label_robustness.csv
FOLLOW   := data/derived/followup_b3_pocket_pair_identity.csv
REDO     := data/derived/cavity_contrasts.csv data/derived/orth_segment.csv
MEAS     := data/derived/nucleophile_measures.csv data/derived/multiview_association.csv
EXPB     := data/derived/experimental_benchmark.csv
study: $(EXPB) $(LABEL) $(FUNC) $(ANCHOR) $(REVIEW) $(FOLLOW) $(REDO) $(MEAS)   ## stages 1-8: universe -> structures -> sites -> comparison -> function -> realignment
$(UFETCH): $(UNIVERSE)
	$(PY) scripts/02_fetch_structures.py --universe
$(UANN): | $(UNIVERSE)
	$(PY) scripts/08_fetch_annotations.py
$(QC): $(UFETCH) $(UANN) scripts/09_structure_qc.py $(SRC)
	$(PY) scripts/09_structure_qc.py
$(SITES): $(QC) scripts/10_catalytic_sites.py $(SRC)
	$(PY) scripts/10_catalytic_sites.py
$(EXPB): $(SITES) scripts/16_experimental_benchmark.py
	$(PY) scripts/16_experimental_benchmark.py
$(LOCAL): $(SITES) scripts/11_local_site.py
	$(PY) scripts/11_local_site.py
$(GLOBAL): $(LOCAL) scripts/12_global_comparison.py
	$(PYE) scripts/12_global_comparison.py
$(ANCHOR): $(SITES) scripts/15_anchored_alignment.py
	$(PY) scripts/15_anchored_alignment.py
$(LABEL): $(GLOBAL) $(ANCHOR) scripts/13_abhd_label_analysis.py
	$(PY) scripts/13_abhd_label_analysis.py
$(FUNC): $(GLOBAL) $(ANCHOR) scripts/14_function_inference.py
	$(PY) scripts/14_function_inference.py

$(REVIEW): $(LABEL) $(FUNC) scripts/17_review_controls.py
	$(PY) scripts/17_review_controls.py

$(FOLLOW): $(REVIEW) scripts/18_pocket_oxyanion.py scripts/19_inactive_vs_active.py scripts/20_panel_pockets.py
	$(PY) scripts/18_pocket_oxyanion.py
	$(PY) scripts/19_inactive_vs_active.py
	$(PY) scripts/20_panel_pockets.py

data/derived/cavity_contrasts.csv: $(FOLLOW) scripts/21_cavity_robustness.py src/abhd/structops.py
	$(PY) scripts/21_cavity_robustness.py

data/derived/orth_segment.csv: $(FOLLOW) scripts/22_orthologues.py src/abhd/structops.py
	$(PY) scripts/22_orthologues.py

data/derived/label_provenance.csv: $(LABEL) scripts/24_label_provenance.py
	$(PY) scripts/24_label_provenance.py

data/derived/function_decay.csv: data/derived/label_provenance.csv $(FUNC) scripts/25_function_decay.py
	$(PY) scripts/25_function_decay.py

data/derived/universe_sources.csv: $(UNIVERSE) scripts/26_universe_sources.py
	$(PY) scripts/26_universe_sources.py

data/derived/core_vs_cath_label_tests.csv: $(GLOBAL) scripts/27_core_domains.py
	$(PY) scripts/27_core_domains.py

data/derived/multiview_association.csv: data/derived/core_vs_cath_label_tests.csv data/derived/function_decay.csv data/derived/universe_sources.csv scripts/28_multiview.py
	$(PY) scripts/28_multiview.py

data/derived/nucleophile_measures.csv: data/derived/core_vs_cath_label_tests.csv scripts/29_nucleophile_measures.py
	$(PY) scripts/29_nucleophile_measures.py

notebooks/.done/%: notebooks/%.ipynb $(PROTBERT) $(ELBOW) $(SUMMARY) $(SOAP) $(FFMAP) $(SRC)
	@mkdir -p notebooks/.done
	uv run jupyter nbconvert --to notebook --execute --inplace --ExecutePreprocessor.timeout=900 $<
	@touch $@

clean-derived:                        ## delete derived tables (keeps raw data and per-protein ProtBert tensors)
	rm -f $(ELBOW) $(SOAP) $(PROTBERT) data/derived/*_cosine.csv figures/*.html
	rm -rf notebooks/.done
