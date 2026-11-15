#!/usr/bin/env python3
"""
Study stage 8 - alignment anchored on the structurally identified catalytic residues.

Instead of aligning whole sequences, residues are placed in columns by their structural role, so no alignment search is
needed:
  anchor_elbow  window of +-WIN residues around the elbow residue (all genes with an elbow call)
  anchor_cat    three windows (+-WIN2) around nucleophile, acid and His, concatenated in chain order
                (genes with a complete triad)
Distance = 1 - BLOSUM62 ungapped score / sqrt(self scores). Compared with whole-domain sequence identity, FunFam and
structure (TM-align) distances: rank correlation, label tests, clustering agreement and EC transfer.
Outputs: data/derived/dist_anchor_elbow.csv, dist_anchor_cat.csv, anchored_summary.csv, anchored_clusters.csv;
data/derived/anchor_windows.csv; figures/stage8_anchored.html
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.elbow import THREE_TO_ONE

WIN, WIN2 = 15, 7
B62 = substitution_matrices.load("BLOSUM62")


def window(res, centre, w):
    out = []
    for k in range(centre - w, centre + w + 1):
        out.append(THREE_TO_ONE.get(res[k].get_resname(), "X") if k in res else "-")
    return "".join(out)


def score(a, b):
    return sum(B62[x, y] for x, y in zip(a, b) if x != "-" and y != "-" and x != "X" and y != "X")


def dist_matrix(seqs):
    accs = list(seqs); n = len(accs); S = np.zeros((n, n)); self_ = [max(score(seqs[a], seqs[a]), 1) for a in accs]
    for i in range(n):
        for j in range(i, n):
            S[i, j] = S[j, i] = score(seqs[accs[i]], seqs[accs[j]])
    sim = np.clip(S / np.sqrt(np.outer(self_, self_)), 0, 1)
    return pd.DataFrame(1 - sim, index=accs, columns=accs)


def main():
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv")
    cs = cs[cs.elbow_resnum.notna()]
    uni = pd.read_csv(C.DERIVED / "universe.csv").set_index("accession")
    parser = PDBParser(QUIET=True)
    elb, cat, rows = {}, {}, []
    for r in cs.itertuples():
        s = parser.get_structure(r.accession, str(C.STRUCTURES / f"{r.accession}.pdb"))
        res = {x.get_id()[1]: x for x in s.get_residues()}
        n = int(r.elbow_resnum)
        elb[r.accession] = window(res, n, WIN)
        if pd.notna(r.his) and pd.notna(r.acid):
            cat[r.accession] = "".join(window(res, int(p), WIN2) for p in sorted([n, int(r.acid), int(r.his)]))
        rows.append(dict(accession=r.accession, gene=r.gene, elbow=elb[r.accession], catalytic=cat.get(r.accession, "")))
    pd.DataFrame(rows).to_csv(C.DERIVED / "anchor_windows.csv", index=False)
    mats = {"anchor_elbow": dist_matrix(elb), "anchor_cat": dist_matrix(cat)}
    for k, M in mats.items():
        M.to_csv(C.DERIVED / f"dist_{k}.csv")

    # concordance with existing spaces (on shared genes)
    out = []
    for k, M in mats.items():
        for other in ("seqid", "tmalign", "funfam", "triad", "soap"):
            O = pd.read_csv(C.DERIVED / f"dist_{other}.csv", index_col=0)
            idx = [a for a in M.index if a in O.index]
            iu = np.triu_indices(len(idx), 1)
            x, y = M.loc[idx, idx].values[iu], O.loc[idx, idx].values[iu]; ok = ~np.isnan(y)
            out.append(dict(anchored=k, vs=other, spearman=round(spearmanr(x[ok], y[ok])[0], 3), n_genes=len(idx)))
    conc = pd.DataFrame(out); print(conc.to_string(index=False))

    # clustering agreement with FunFam (k=30) for anchored vs whole-sequence spaces
    ffl = uni.funfam_ids.fillna("").str.split(";").str[0]
    cl = []
    for k in ("anchor_elbow", "anchor_cat", "seqid", "tmalign"):
        M = mats[k] if k in mats else pd.read_csv(C.DERIVED / f"dist_{k}.csv", index_col=0)
        accs = [a for a in M.index if ffl.get(a, "")]
        v = M.loc[accs, accs].values; v = (v + v.T) / 2; np.fill_diagonal(v, 0)
        Z = linkage(squareform(v, checks=False), method="average")
        for nc in (20, 30):
            c = fcluster(Z, nc, "maxclust")
            cl.append(dict(space=k, n_clusters=nc, n_genes=len(accs), ARI_vs_funfam=round(adjusted_rand_score(ffl[accs], c), 3),
                           NMI_vs_funfam=round(normalized_mutual_info_score(ffl[accs], c), 3)))
    cl = pd.DataFrame(cl); print(cl.to_string(index=False))
    conc.to_csv(C.DERIVED / "anchored_summary.csv", index=False); cl.to_csv(C.DERIVED / "anchored_clusters.csv", index=False)

    fig = go.Figure()
    for sp in ("anchor_elbow", "anchor_cat", "seqid", "tmalign"):
        d = cl[(cl.space == sp) & (cl.n_clusters == 30)]
        fig.add_bar(x=[sp], y=d.ARI_vs_funfam, name=sp)
    fig.update_layout(template="simple_white", width=700, height=400, showlegend=False,
                      title="Stage 8: agreement with FunFam (ARI, 30 clusters)", yaxis_title="ARI")
    fig.write_html(C.FIGURES / "stage8_anchored.html")


if __name__ == "__main__":
    main()
