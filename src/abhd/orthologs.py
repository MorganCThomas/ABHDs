"""Vertebrate orthologue conservation for any gene in the universe (the pipeline of script 22, for one gene).

A UniProt entry is kept as an orthologue only if its best BLOSUM62 local hit among all human universe sequences is the target gene (10% margin over the
runner-up) with >= 40% identity and >= 70% coverage of the human domain; one sequence (the best scoring) per organism, human excluded.
Needs network access the first time (UniProt REST); results are cached in data/raw/orthologues/<GENE>_v2.tsv.
"""
from multiprocessing import Pool

import numpy as np
import pandas as pd
import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from scipy.stats import rankdata

from . import config as C
from .structops import load

ORTH = C.RAW / "orthologues"
CLASSES = [("Mammalia", "mammals"), ("Aves", "birds"), ("Crocodylia", "crocodilians"), ("Testudines", "turtles"), ("Lepidosauria", "lepidosaurs"),
           ("Amphibia", "amphibians"), ("Chondrichthyes", "cartilaginous fish"), ("Coelacanthimorpha", "lobe-finned fish"), ("Dipnoi", "lobe-finned fish"),
           ("Actinopterygii", "ray-finned fish"), ("Hyperoartia", "jawless fish"), ("Cyclostomata", "jawless fish")]
_AL = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-10, extend_gap_score=-0.5)
_HUMAN = {}


def _init(h):
    _HUMAN.update(h)


def _score_all(seq):
    return {g: _AL.score(seq, s) for g, s in _HUMAN.items()}


def taxon_class(lineage):
    for tag, name in CLASSES:
        if tag in str(lineage):
            return name
    return "other vertebrate"


def human_sequences(atlas):
    seqs = {}
    for f in (C.ROOT / "archive" / "uniprot_family_query_2026-11-11" / "abhd_sequences.fasta", C.SEQUENCES_FASTA, C.RAW / "uniprot" / "funfam_only_sequences.fasta"):
        if f.exists():
            for rec in SeqIO.parse(f, "fasta"):
                seqs[rec.id.split("|")[1]] = str(rec.seq)
    return {g: seqs[a] for g, a in atlas.acc_of.items() if a in seqs}


def fetch(gene):
    ORTH.mkdir(parents=True, exist_ok=True)
    p = ORTH / f"{gene}_v2.tsv"
    if not p.exists():
        url = "https://rest.uniprot.org/uniprotkb/search"
        params = {"query": f"gene:{gene} AND taxonomy_id:7742", "format": "tsv", "size": 500, "fields": "accession,organism_name,length,reviewed,sequence,lineage"}
        parts, r = [], requests.get(url, params=params, timeout=90)
        while True:
            r.raise_for_status()
            parts.append(r.text if not parts else r.text.split("\n", 1)[1])
            nxt = r.links.get("next", {}).get("url")
            if not nxt:
                break
            r = requests.get(nxt, timeout=90)
        p.write_text("".join(parts))
    return pd.read_csv(p, sep="\t").rename(columns={"Taxonomic lineage": "Lineage"})


def conservation(atlas, gene):
    """Per-position conservation across cleaned vertebrate orthologues. Returns a dict (or {'error': ...})."""
    a = atlas.acc(gene); g = atlas.gene(a)
    try:
        df = fetch(g)
    except Exception as e:                                   # no network
        return dict(error=f"orthologues could not be fetched ({type(e).__name__}); cached file missing")
    human = human_sequences(atlas)
    df = df[(df.Entry != a) & (df.Organism != "Homo sapiens (Human)") & df.Sequence.notna()].copy()
    with Pool(initializer=_init, initargs=(human,)) as p:
        sc = p.map(_score_all, list(df.Sequence), chunksize=8)
    keep = []
    for row, s in zip(df.itertuples(), sc):
        ranked = sorted(s.items(), key=lambda kv: -kv[1])
        if ranked[0][0] == g and (len(ranked) < 2 or ranked[0][1] >= 1.1 * ranked[1][1]):
            keep.append(dict(accession=row.Entry, organism=row.Organism, score=ranked[0][1], sequence=row.Sequence, lineage=row.Lineage))
    K = pd.DataFrame(keep)
    if K.empty:
        return dict(error="no orthologues passed the paralog filter", n_candidates=len(df))
    K = K.sort_values("score", ascending=False).drop_duplicates("organism")
    ids, _, hseq = load(a)
    cols, meta = [], []
    for r in K.itertuples():
        al = _AL.align(hseq, r.sequence)[0]
        m = ["-"] * len(hseq); nal = nid = 0
        for (s0, e0), (s1, e1) in zip(*al.aligned):
            for k in range(e0 - s0):
                m[s0 + k] = r.sequence[s1 + k]; nal += 1; nid += hseq[s0 + k] == r.sequence[s1 + k]
        if nal / len(hseq) >= 0.7 and nid / max(nal, 1) >= 0.4:
            cols.append("".join(m)); meta.append(dict(organism=r.organism, taxon_class=taxon_class(r.lineage), identity=nid / nal))
    if len(cols) < 6:
        return dict(error=f"only {len(cols)} orthologues passed the filters (need 6)", n_candidates=len(df))
    M = np.array([list(c) for c in cols]); meta = pd.DataFrame(meta)
    aligned = M != "-"
    cons = np.array([np.mean(M[:, k][aligned[:, k]] == hseq[k]) if aligned[:, k].any() else np.nan for k in range(len(hseq))])
    ok = ~np.isnan(cons); pct = np.full(len(cons), np.nan); pct[ok] = (rankdata(cons[ok]) - 0.5) / ok.sum()
    per_pos = pd.DataFrame(dict(resnum=ids, aa=list(hseq), conservation=cons, n_aligned=aligned.sum(0), percentile=pct))
    # key positions: nucleophile (elbow), His, acid; residues by taxon class
    c = atlas.cs.loc[a]
    roles = {"nucleophile position": c.elbow_resnum, "His": c.his, "acid": c.acid}
    rows = []
    for role, num in roles.items():
        if pd.isna(num) or int(num) not in ids:
            continue
        k = int(np.where(ids == int(num))[0][0])
        by = pd.Series(M[:, k]).groupby(meta.taxon_class.values).agg(lambda s: ";".join(f"{x}:{n}" for x, n in s.value_counts().items()))
        rows.append(dict(role=role, human=f"{hseq[k]}{int(num)}", conservation=float(cons[k]), midrank_percentile=float(pct[k]),
                         residues_by_class="; ".join(f"{cls} {v}" for cls, v in by.items())))
    return dict(per_position=per_pos, key_positions=pd.DataFrame(rows), n_orthologues=len(meta), n_candidates=len(df),
                classes=meta.taxon_class.value_counts().rename_axis("taxon_class").reset_index(name="n"))
