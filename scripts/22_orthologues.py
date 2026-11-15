#!/usr/bin/env python3
"""
Orthologue clean-up (reports/12): vertebrate orthologues of ABHD1/2/3/4/5/15 with
  * paginated UniProt retrieval (all pages, cached in data/raw/orthologues/<gene>_v2.tsv),
  * assignment by best hit against all 110 human universe sequences (a candidate is kept only if its best human hit is the target
    gene with a >=10% score margin over the runner-up, and >=40% identity / >=70% coverage of the human domain),
  * the single highest-scoring sequence per organism, human excluded,
  * taxonomic class from the UniProt lineage, class-level counts and conservation,
  * mid-rank conservation percentiles (ties are not penalised), and the ABHD4 180-198 / ABHD5-equivalent segment.
Outputs: data/derived/orthologue_sets.csv, orth_conservation.csv, orth_elbow_by_class.csv, orth_segment.csv, orth_rejected.csv
"""
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from scipy.stats import rankdata

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.structops import align, load

ORTH = C.RAW / "orthologues"
GENES = ("ABHD15", "ABHD5", "ABHD4", "ABHD3", "ABHD2", "ABHD1")
CLASSES = [("Mammalia", "mammals"), ("Aves", "birds"), ("Crocodylia", "crocodilians"), ("Testudines", "turtles"),
           ("Lepidosauria", "lepidosaurs"), ("Amphibia", "amphibians"), ("Chondrichthyes", "cartilaginous fish"),
           ("Coelacanthimorpha", "lobe-finned fish"), ("Dipnoi", "lobe-finned fish"), ("Actinopterygii", "ray-finned fish"),
           ("Hyperoartia", "jawless fish"), ("Cyclostomata", "jawless fish")]
_AL = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-10, extend_gap_score=-0.5)
HUMAN = {}


def klass(lineage):
    for tag, name in CLASSES:
        if tag in str(lineage):
            return name
    return "other vertebrate"


def fetch(gene):
    ORTH.mkdir(parents=True, exist_ok=True)
    p = ORTH / f"{gene}_v2.tsv"
    if p.exists():
        return pd.read_csv(p, sep="\t").rename(columns={"Taxonomic lineage": "Lineage"})
    url = "https://rest.uniprot.org/uniprotkb/search"
    params = {"query": f"gene:{gene} AND taxonomy_id:7742", "format": "tsv", "size": 500,
              "fields": "accession,organism_name,length,reviewed,sequence,lineage"}
    parts, r = [], requests.get(url, params=params, timeout=90)
    while True:
        r.raise_for_status(); parts.append(r.text if not parts else r.text.split("\n", 1)[1])
        nxt = r.links.get("next", {}).get("url")
        if not nxt:
            break
        r = requests.get(nxt, timeout=90)
    p.write_text("".join(parts))
    return pd.read_csv(p, sep="\t").rename(columns={"Taxonomic lineage": "Lineage"})


def _init(h):
    HUMAN.update(h)


def score_all(seq):
    return {g: _AL.score(seq, s) for g, s in HUMAN.items()}


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv")
    seqs = {}
    for f in (C.ROOT / "archive" / "uniprot_family_query_2026-11-11" / "abhd_sequences.fasta", C.SEQUENCES_FASTA,
              C.RAW / "uniprot" / "funfam_only_sequences.fasta"):
        for rec in SeqIO.parse(f, "fasta"):
            seqs[rec.id.split("|")[1]] = str(rec.seq)
    for g, a in zip(uni.gene, uni.accession):
        HUMAN[g] = seqs[a]
    gacc = dict(zip(uni.gene, uni.accession))
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv").set_index("gene")
    sets, rej, cons_rows, elbow_rows = [], [], [], []
    store = {}
    for gene in GENES:
        df = fetch(gene)
        df = df[(df.Entry != gacc[gene]) & (df.Organism != "Homo sapiens (Human)") & df.Sequence.notna()].copy()
        with Pool(initializer=_init, initargs=(HUMAN,)) as p:
            sc = p.map(score_all, list(df.Sequence), chunksize=8)
        keep = []
        for row, s in zip(df.itertuples(), sc):
            ranked = sorted(s.items(), key=lambda kv: -kv[1])
            best, second = ranked[0], ranked[1]
            ok = best[0] == gene and best[1] >= 1.1 * second[1]
            rec = dict(target=gene, accession=row.Entry, organism=row.Organism, best_human_hit=best[0], best=round(best[1], 1),
                       runner_up=ranked[1][0], runner_score=round(second[1], 1), target_score=round(s[gene], 1))
            (keep if ok else rej).append(rec | dict(Sequence=row.Sequence, lineage=row.Lineage, reviewed=row.Reviewed))
        K = pd.DataFrame(keep)
        if K.empty:
            continue
        K = K.sort_values("best", ascending=False).drop_duplicates("organism")          # best-scoring sequence per organism
        ids, _, hseq = load(gacc[gene])
        cols, meta = [], []
        for r in K.itertuples():
            a = _AL.align(hseq, r.Sequence)[0]
            m = ["-"] * len(hseq); nal = nid = 0
            for (s0, e0), (s1, e1) in zip(*a.aligned):
                for k in range(e0 - s0):
                    m[s0 + k] = r.Sequence[s1 + k]; nal += 1; nid += hseq[s0 + k] == r.Sequence[s1 + k]
            if nal / len(hseq) < 0.7 or nid / max(nal, 1) < 0.4:
                rej.append(dict(target=gene, accession=r.accession, organism=r.organism, best_human_hit="(identity/coverage)", best=r.best,
                                runner_up="", runner_score=0, target_score=r.target_score, Sequence="", lineage=r.lineage, reviewed=r.reviewed))
                continue
            cols.append("".join(m)); meta.append(dict(target=gene, accession=r.accession, organism=r.organism, taxon_class=klass(r.lineage),
                                                      reviewed=r.reviewed, identity=round(nid / nal, 3), coverage=round(nal / len(hseq), 3)))
        M = np.array([list(c) for c in cols]); meta = pd.DataFrame(meta)
        sets.append(meta); store[gene] = (ids, hseq, M, meta)
        print(f"{gene}: candidates {len(df)}, kept {len(meta)} (rejected as other paralog / low score: {len(df) - len(K)}); classes:",
              meta.taxon_class.value_counts().to_dict(), flush=True)
    pd.concat(sets).to_csv(C.DERIVED / "orthologue_sets.csv", index=False)
    pd.DataFrame(rej).drop(columns=["Sequence"], errors="ignore").to_csv(C.DERIVED / "orth_rejected.csv", index=False)

    for gene, (ids, hseq, M, meta) in store.items():
        aligned = M != "-"
        cons = np.array([np.mean(M[:, k][aligned[:, k]] == hseq[k]) if aligned[:, k].any() else np.nan for k in range(len(hseq))])
        mid = (rankdata(cons[~np.isnan(cons)]) - 0.5) / (~np.isnan(cons)).sum()
        pct = np.full(len(cons), np.nan); pct[~np.isnan(cons)] = mid
        el = int(cs.loc[gene, "elbow_resnum"]); ie = int(np.where(ids == el)[0][0])
        keys = {f"elbow{o:+d}": ie + o for o in range(-2, 3)}
        extra = {"ABHD15": {"his_aligned": 391, "acid_aligned": 360}, "ABHD5": {"his_aligned": 327, "acid_aligned": 177}}.get(gene, {})
        for role in ("his", "acid"):
            v = cs.loc[gene, role]
            if pd.notna(v) and int(v) in ids:
                keys[role] = int(np.where(ids == int(v))[0][0])
        for role, num in extra.items():
            keys[role] = int(np.where(ids == num)[0][0])
        for name, k in keys.items():
            cons_rows.append(dict(gene=gene, n=len(meta), position=name, human=f"{hseq[k]}{int(ids[k])}", n_aligned=int(aligned[:, k].sum()),
                                  conservation=round(float(cons[k]), 3), domain_median=round(float(np.nanmedian(cons)), 3),
                                  midrank_percentile=round(float(pct[k]), 3)))
        by = meta.assign(res=M[:, ie])
        for cls, g in by.groupby("taxon_class"):
            elbow_rows.append(dict(gene=gene, taxon_class=cls, n=len(g), human_residue=hseq[ie],
                                   residues=";".join(f"{a}:{c}" for a, c in g.res.value_counts().items())))
    pd.DataFrame(cons_rows).to_csv(C.DERIVED / "orth_conservation.csv", index=False)
    E = pd.DataFrame(elbow_rows); E.to_csv(C.DERIVED / "orth_elbow_by_class.csv", index=False)
    pd.set_option("display.width", 220)
    print(pd.DataFrame(cons_rows).to_string(index=False)); print(E.to_string(index=False))

    # sequences at the elbow that differ from the human residue, with organism
    odd = []
    for gene in ("ABHD15", "ABHD5"):
        ids, hseq, M, meta = store[gene]; ie = int(np.where(ids == int(cs.loc[gene, "elbow_resnum"]))[0][0])
        for r, a in zip(meta.itertuples(), M[:, ie]):
            if a != hseq[ie]:
                odd.append(dict(gene=gene, organism=r.organism, taxon_class=r.taxon_class, accession=r.accession, residue=a,
                                next=M[meta.index.get_loc(r.Index), ie + 1], identity=r.identity))
    O = pd.DataFrame(odd); O.to_csv(C.DERIVED / "orth_elbow_exceptions.csv", index=False); print(O.to_string(index=False))

    # ABHD4 180-198 segment vs the TM-aligned ABHD5 segment
    (i5, c5, s5), (i4, c4, s4), X, pairs, r = align(gacc["ABHD5"], gacc["ABHD4"])
    map45 = {int(i4[j]): int(i5[i]) for i, j in pairs}                         # ABHD4 residue -> ABHD5 residue
    ids4, h4, M4, m4 = store["ABHD4"]; ids5, h5, M5, m5 = store["ABHD5"]
    rows = []
    for pos4 in range(180, 199):
        k4 = int(np.where(ids4 == pos4)[0][0]); pos5 = map45.get(pos4)
        col4 = M4[:, k4][M4[:, k4] != "-"]
        row = dict(abhd4_pos=f"{h4[k4]}{pos4}", abhd4_conservation=round(float(np.mean(col4 == h4[k4])), 3) if len(col4) else np.nan,
                   abhd4_consensus=pd.Series(col4).mode().iat[0] if len(col4) else "-")
        if pos5 is not None:
            k5 = int(np.where(ids5 == pos5)[0][0]); col5 = M5[:, k5][M5[:, k5] != "-"]
            row.update(abhd5_pos=f"{h5[k5]}{pos5}", abhd5_conservation=round(float(np.mean(col5 == h5[k5])), 3) if len(col5) else np.nan,
                       abhd5_consensus=pd.Series(col5).mode().iat[0] if len(col5) else "-")
        else:
            row.update(abhd5_pos="unaligned")
        rows.append(row)
    G = pd.DataFrame(rows); G.to_csv(C.DERIVED / "orth_segment.csv", index=False); print(G.to_string(index=False))
    d4 = np.nanmedian([r for r in G.abhd4_conservation]); d5 = np.nanmedian([r for r in G.abhd5_conservation.dropna()])
    print(f"segment median conservation: ABHD4 {d4:.2f}, ABHD5 (aligned positions) {d5:.2f}; unaligned in ABHD5: {(G.abhd5_pos == 'unaligned').sum()}/19")


if __name__ == "__main__":
    main()
