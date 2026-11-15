#!/usr/bin/env python3
"""
Study stage 1 - build the merged protein universe and its attrition ledger.

Sources (each protein keeps a record of where it came from):
  uniprot_2025  frozen UniProt family query (data/raw/uniprot/accessions.txt, 210 entries)
  uniprot_2026  same query re-run on 2026-11-11 (archive/uniprot_family_query_2026-11-11, 378 entries)
  funfam        human members of CATH 3.40.50.1820 FunFams (data/raw/christina/funfam_metadata.csv)

Collapsing rule (one row per gene):
  * a reviewed (Swiss-Prot) entry with a gene name is the canonical entry for that gene;
  * unreviewed entries with the same gene name are alternative entries (reason: `alt_entry_of_<gene>`);
  * unreviewed entries without a gene name are aligned to the canonical sequences; if >= FRAG_MIN_IDENTITY identity
    over >= FRAG_MIN_COVERAGE of the fragment they are assigned to that gene (reason: `fragment_of_<gene>`),
    otherwise they are dropped as `unassigned_fragment`;
  * unreviewed entries with a gene name that has no reviewed entry are dropped as `unreviewed_only_gene`
    (listed in the ledger so they can be reviewed by hand).

Outputs:
  data/raw/uniprot/funfam_only_entries.json   metadata for FunFam accessions missing from the UniProt query (cached)
  data/derived/universe_entries.csv           one row per source accession with its fate (the attrition ledger)
  data/derived/universe.csv                   one row per gene: canonical accession and source flags
  data/derived/attrition_stage1.csv           counts per fate
"""
import json
import sys
from pathlib import Path

import pandas as pd
import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd import data as D

ARCHIVE = C.ROOT / "archive" / "uniprot_family_query_2026-11-11"
FUNFAM_ONLY_JSON = C.RAW / "uniprot" / "funfam_only_entries.json"
FUNFAM_ONLY_FASTA = C.RAW / "uniprot" / "funfam_only_sequences.fasta"
SEARCH = "https://rest.uniprot.org/uniprotkb/search"

FRAG_MIN_IDENTITY = 0.90
FRAG_MIN_COVERAGE = 0.80


def fetch_funfam_only(accs):
    """Entry metadata + sequences for FunFam accessions that the UniProt family query does not return."""
    if FUNFAM_ONLY_JSON.exists() and FUNFAM_ONLY_FASTA.exists():
        return json.load(open(FUNFAM_ONLY_JSON))
    q = "accession:(" + " OR ".join(sorted(accs)) + ")"
    r = requests.get(SEARCH, timeout=60, params={
        "query": q, "format": "tsv", "size": 500,
        "fields": "accession,id,protein_name,gene_primary,cc_catalytic_activity,cc_function,ec,reviewed,length"})
    r.raise_for_status()
    rows = []
    for line in r.text.rstrip("\n").splitlines()[1:]:
        p = line.split("\t") + [""] * 9
        rows.append(dict(accession=p[0], entry_name=p[1], name=p[2], gene_name=p[3], catalytic_activity=p[4],
                         function=p[5], ec=p[6], reviewed=p[7], length=int(p[8] or 0)))
    f = requests.get(SEARCH, timeout=60, params={"query": q, "format": "fasta", "size": 500})
    f.raise_for_status()
    FUNFAM_ONLY_JSON.write_text(json.dumps(rows, indent=1))
    FUNFAM_ONLY_FASTA.write_text(f.text)
    return rows


def read_fasta(path):
    return {rec.id.split("|")[1]: str(rec.seq) for rec in SeqIO.parse(path, "fasta")}


def best_gene_for_fragment(seq, canon_seqs, aligner):
    best = (0.0, 0.0, None)
    for gene, ref in canon_seqs.items():
        aln = aligner.align(seq, ref)[0]
        ident = sum(a == b for (s0, e0), (s1, e1) in zip(*aln.aligned) for a, b in zip(seq[s0:e0], ref[s1:e1]))
        aligned = sum(e0 - s0 for s0, e0 in aln.aligned[0])
        identity = ident / aligned if aligned else 0.0
        coverage = aligned / len(seq)
        if identity * coverage > best[0] * best[1]:
            best = (identity, coverage, gene)
    return best


def make_figures(uni, att):
    """figures/stage1_*.html: source overlap by ABHD name, and the entry-level attrition ledger."""
    import plotly.graph_objects as go
    C.FIGURES.mkdir(exist_ok=True)
    order = ["uniprot_only", "both", "funfam_only"]
    labels = {"uniprot_only": "UniProt query only", "both": "UniProt + FunFam", "funfam_only": "FunFam only"}
    fig = go.Figure()
    for named, name, col in ((True, "ABHD-named", "crimson"), (False, "other name", "lightgrey")):
        sub = uni[uni.abhd_named == named]
        fig.add_bar(y=[labels[o] for o in order], x=[(sub.source == o).sum() for o in order], orientation="h",
                    name=name, marker_color=col,
                    hovertext=[", ".join(sorted(sub[sub.source == o].gene)) for o in order], hoverinfo="text+x")
    fig.update_layout(barmode="stack", template="simple_white", width=800, height=330,
                      title="Stage 1: genes in the merged universe by source", xaxis_title="genes (canonical reviewed entry)")
    fig.write_html(C.FIGURES / "stage1_universe_overlap.html")
    a = att.sort_values("n_entries")
    fig = go.Figure(go.Bar(y=a.fate_class, x=a.n_entries, orientation="h", marker_color="steelblue"))
    fig.update_layout(template="simple_white", width=800, height=330, title="Stage 1: fate of every source accession",
                      xaxis_title="accessions (UniProt 2025 + 2026 + FunFam, deduplicated)")
    fig.write_html(C.FIGURES / "stage1_attrition.html")


def main():
    C.DERIVED.mkdir(parents=True, exist_ok=True)
    e25 = {e["accession"]: e for e in json.load(open(C.ENTRIES_JSON))}
    e26 = {e["accession"]: e for e in json.load(open(ARCHIVE / "abhd_entries.json"))}
    ff, _ = D.load_funfam()
    ff_acc = {}                                     # accession -> FunFam ids
    for fid, s in zip(ff.funfam_id, ff.human_uniprots):
        for a in s.split(";"):
            ff_acc.setdefault(a, []).append(fid)

    missing = set(ff_acc) - set(e25) - set(e26)
    e_ff = {e["accession"]: e for e in fetch_funfam_only(missing)}
    seqs = {**read_fasta(ARCHIVE / "abhd_sequences.fasta"), **read_fasta(C.SEQUENCES_FASTA),
            **read_fasta(FUNFAM_ONLY_FASTA)}

    allent = {**e26, **e25, **e_ff}                  # metadata, 2025 wins over 2026 where both exist
    rows = []
    for acc, e in allent.items():
        rows.append(dict(
            accession=acc, entry_name=e["entry_name"], gene_name=e["gene_name"], reviewed=e["reviewed"] == "reviewed",
            name=e["name"], ec=e["ec"], length=len(seqs.get(acc, "")),
            in_uniprot_2025=acc in e25, in_uniprot_2026=acc in e26, in_funfam=acc in ff_acc,
            funfam_ids=";".join(ff_acc.get(acc, []))))
    ent = pd.DataFrame(rows).set_index("accession", drop=False)

    # canonical = reviewed entry with a gene name
    canon = ent[ent.reviewed & ent.gene_name.ne("")]
    assert not canon.gene_name.duplicated().any(), "two reviewed entries share a gene name"
    canon_by_gene = canon.set_index("gene_name").accession.to_dict()
    canon_seqs = {g: seqs[a] for g, a in canon_by_gene.items() if a in seqs}

    aligner = PairwiseAligner(mode="local", match_score=2, mismatch_score=-1, open_gap_score=-4, extend_gap_score=-1)
    fate, gene_of = {}, {}
    for acc, r in ent.iterrows():
        if acc in set(canon.accession):
            fate[acc], gene_of[acc] = "canonical", r.gene_name
        elif r.gene_name and r.gene_name in canon_by_gene:
            fate[acc], gene_of[acc] = f"alt_entry_of_{r.gene_name}", r.gene_name
        else:
            # no canonical entry carries this entry's gene name (blank, alias such as BAT5, or no reviewed entry):
            # try to place it by sequence before giving up
            ident, cov, g = best_gene_for_fragment(seqs[acc], canon_seqs, aligner) if acc in seqs else (0, 0, None)
            if g and ident >= FRAG_MIN_IDENTITY and cov >= FRAG_MIN_COVERAGE:
                fate[acc], gene_of[acc] = f"fragment_of_{g}", g
            elif r.gene_name:
                fate[acc], gene_of[acc] = "unreviewed_only_gene", r.gene_name
            else:
                fate[acc], gene_of[acc] = "unassigned_fragment", ""
    ent["fate"] = pd.Series(fate)
    ent["assigned_gene"] = pd.Series(gene_of)
    ent.to_csv(C.DERIVED / "universe_entries.csv", index=False)

    # gene-level universe: canonical entries, with source flags pooled over all entries assigned to the gene
    g = ent[ent.assigned_gene.ne("")].groupby("assigned_gene")
    uni = canon.set_index("gene_name")[["accession", "entry_name", "name", "ec", "length", "funfam_ids"]].copy()
    uni.index.name = "gene"
    for flag in ("in_uniprot_2025", "in_uniprot_2026", "in_funfam"):
        uni[flag] = g[flag].any().reindex(uni.index).fillna(False)
    uni["n_entries"] = g.size().reindex(uni.index)
    uni["abhd_named"] = uni.index.str.startswith("ABHD")          # strict: gene symbol starts with ABHD
    # broad: UniProt protein name carries "abhydrolase domain-containing protein N" (adds EPHX3=ABHD9, EPHX4=ABHD7 and the
    # LIPJ/K/M/N "lipase-like abhydrolase domain-containing" genes)
    uni["abhd_family"] = uni.name.fillna("").str.contains("abhydrolase domain-containing", case=False)
    uni["source"] = uni.apply(lambda r: ("both" if r.in_funfam and (r.in_uniprot_2025 or r.in_uniprot_2026) else
                                         "funfam_only" if r.in_funfam else "uniprot_only"), axis=1)
    uni.reset_index().to_csv(C.DERIVED / "universe.csv", index=False)

    led = ent.assign(fate_class=ent.fate.str.replace(r"_of_.*|_entry_of_.*", "", regex=True).replace(
        {"alt": "alt_entry", "alt_entry": "alt_entry"}))
    led["fate_class"] = ent.fate.apply(lambda f: f.split("_of_")[0] if "_of_" in f else f)
    att = led.groupby("fate_class").size().rename("n_entries").reset_index()
    att.to_csv(C.DERIVED / "attrition_stage1.csv", index=False)
    print(att.to_string(index=False))
    make_figures(uni.reset_index(), att)
    print(f"\nuniverse: {len(uni)} genes; ABHD-named {int(uni.abhd_named.sum())}")
    print(pd.crosstab(uni.source, uni.abhd_named, margins=True).to_string())


if __name__ == "__main__":
    main()
