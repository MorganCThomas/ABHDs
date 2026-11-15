#!/usr/bin/env python3
"""
Study stage 6 - global comparison of the universe in several spaces, and concordance with FunFam and the ABHD name.

Spaces (gene x gene distances, all on the stage-2 hydrolase domain):
  tmalign   structure: 1 - mean TM-score of both normalisations (tmtools / TM-align)
  seqid     sequence: 1 - BLOSUM62 Smith-Waterman local score normalised by sqrt(self scores)
  tmshort   structure: 1 - TM-score normalised by the shorter chain (length-robust)
  protbert  ProtBert mean embedding, cosine distance (domain sequence, <=510 residues)
  funfam    HHalign E-values between FunFam HMMs (-log10, min over the FunFams a gene belongs to); NaN if no FunFam
  soap      stage-5 SOAP at the elbow CA, cosine distance
  triad     stage-5 triad skeleton RMSD (genes with a complete triad)
Outputs (data/derived/): dist_<space>.csv, global_label_tests.csv, global_concordance.csv, global_clusters.csv,
discordant_pairs.csv; figures/stage6_*.html
"""
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from Bio import SeqIO
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser, PPBuilder
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform
from scipy.stats import spearmanr
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from tmtools import tm_align

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd import data as D

DOMAINS = C.DERIVED / "domains"
K_NN = 3
N_PERM = 5000
rng = np.random.default_rng(0)


def load_domain(acc):
    s = PDBParser(QUIET=True).get_structure(acc, str(DOMAINS / f"{acc}.pdb"))
    ca, seq = [], []
    from abhd.elbow import THREE_TO_ONE
    for r in s.get_residues():
        if "CA" in r:
            ca.append(r["CA"].coord); seq.append(THREE_TO_ONE.get(r.get_resname(), "X"))
    return np.array(ca), "".join(seq)


def _tm(args):
    (a, ca1, s1), (b, ca2, s2) = args
    r = tm_align(ca1, ca2, s1, s2)
    return a, b, (r.tm_norm_chain1 + r.tm_norm_chain2) / 2, max(r.tm_norm_chain1, r.tm_norm_chain2)


def tm_matrix(data):
    accs = list(data)
    jobs = [((a, *data[a]), (b, *data[b])) for i, a in enumerate(accs) for b in accs[i + 1:]]
    pos = {a: i for i, a in enumerate(accs)}
    M = np.eye(len(accs)); Ms = np.eye(len(accs))
    with Pool() as p:
        for a, b, t, ts in p.imap_unordered(_tm, jobs, chunksize=50):
            M[pos[a], pos[b]] = M[pos[b], pos[a]] = t
            Ms[pos[a], pos[b]] = Ms[pos[b], pos[a]] = ts     # normalised by the shorter chain: length-robust
    return pd.DataFrame(1 - M, index=accs, columns=accs), pd.DataFrame(1 - Ms, index=accs, columns=accs)


def seq_matrix(data):
    al = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"),
                         open_gap_score=-11, extend_gap_score=-1)
    accs = list(data); seqs = {a: data[a][1].replace("X", "A") for a in accs}
    self_ = {a: al.score(seqs[a], seqs[a]) for a in accs}
    M = pd.DataFrame(0.0, index=accs, columns=accs)
    for i, a in enumerate(accs):
        for b in accs[i + 1:]:
            M.loc[a, b] = M.loc[b, a] = 1 - al.score(seqs[a], seqs[b]) / np.sqrt(self_[a] * self_[b])
    return M


def protbert_matrix(data):
    import torch
    from transformers import BertModel, BertTokenizer
    tok = BertTokenizer.from_pretrained(C.PROTBERT_MODEL, do_lower_case=False)
    model = BertModel.from_pretrained(C.PROTBERT_MODEL).eval()
    embs = {}
    for a, (_, seq) in data.items():
        text = " ".join(seq[:510].replace("U", "X").replace("Z", "X").replace("O", "X").replace("B", "X"))
        enc = tok(text, return_tensors="pt")
        with torch.no_grad():
            h = model(**enc).last_hidden_state[0, 1:-1].mean(0).numpy()
        embs[a] = h
    X = np.array([embs[a] for a in data]); X = X / np.linalg.norm(X, axis=1, keepdims=True)
    return pd.DataFrame(1 - X @ X.T, index=list(data), columns=list(data))


def funfam_matrix(uni):
    ff, ev = D.load_funfam()
    d = D.evalue_to_distance(ev)
    accs = list(uni.accession)
    sets = {a: [f for f in str(uni.set_index("accession").loc[a, "funfam_ids"]).split(";") if f in d.index] for a in accs}
    M = pd.DataFrame(np.nan, index=accs, columns=accs)
    for a in accs:
        for b in accs:
            if sets[a] and sets[b]:
                M.loc[a, b] = np.min([d.loc[x, y] for x in sets[a] for y in sets[b]])
    return M


def sym(M):
    v = M.values; return pd.DataFrame(np.fmin(v, v.T) if np.isnan(v).any() else (v + v.T) / 2, index=M.index, columns=M.columns)


def perm_test(Dm, labels):
    """mean distance within the labelled group minus labelled-vs-rest; permutation p for 'closer than chance'."""
    v = Dm.values.copy(); ok = ~np.isnan(v); np.fill_diagonal(ok, False)
    lab = np.array(labels, bool)
    def stat(l):
        w = np.outer(l, l) & ok; x = np.outer(l, ~l) & ok
        if w.sum() == 0 or x.sum() == 0: return np.nan
        return v[w].mean() - v[x].mean()
    obs = stat(lab)
    null = np.array([stat(rng.permutation(lab)) for _ in range(N_PERM)])
    return obs, float((np.sum(null <= obs) + 1) / (N_PERM + 1)), float(np.nanmean(null))


def knn_purity(Dm, labels, k=K_NN):
    v = Dm.values.copy(); np.fill_diagonal(v, np.inf); v[np.isnan(v)] = np.inf
    lab = np.array(labels, bool); out = []
    for i in np.where(lab)[0]:
        order = np.argsort(v[i])[:k]; order = [j for j in order if np.isfinite(v[i, j])]
        if order: out.append(np.mean(lab[order]))
    return float(np.mean(out)) if out else np.nan


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv")
    qc = pd.read_csv(C.DERIVED / "structure_qc.csv")
    gene = dict(zip(uni.accession, uni.gene)); abhd = dict(zip(uni.accession, uni.abhd_named))
    accs = [a for a in uni.accession if (DOMAINS / f"{a}.pdb").exists()]
    data = {a: load_domain(a) for a in accs}
    mats = {}
    mats["tmalign"], mats["tmshort"] = tm_matrix(data); print("tmalign done", flush=True)
    mats["seqid"] = seq_matrix(data); print("seqid done", flush=True)
    mats["protbert"] = protbert_matrix(data); print("protbert done", flush=True)
    mats["funfam"] = funfam_matrix(uni[uni.accession.isin(accs)])
    soap = pd.read_csv(C.DERIVED / "soap_site_cosine.csv", index_col=0); mats["soap"] = soap.reindex(index=accs, columns=accs)
    tri = pd.read_csv(C.DERIVED / "site_rmsd_triad.csv", index_col=0); mats["triad"] = tri.reindex(index=accs, columns=accs)
    for k, M in mats.items():
        M = M.reindex(index=accs, columns=accs); mats[k] = M; M.to_csv(C.DERIVED / f"dist_{k}.csv")
    labels = np.array([abhd[a] for a in accs])

    # ---- label tests: are ABHD-named genes closer to each other than to others? ----
    rows = []
    n_abhd, n = labels.sum(), len(labels)
    for k, M in mats.items():
        sel = np.where(~np.isnan(M.values).all(1))[0]          # genes with any distance in this space
        lab = labels[sel]; Ms = M.iloc[sel, sel]
        obs, p, nullmean = perm_test(Ms, lab)
        pur = knn_purity(Ms, lab)
        rows.append(dict(space=k, n_genes=len(sel), n_abhd=int(lab.sum()), mean_within_minus_between=round(obs, 4),
                         perm_p_closer=round(p, 4), knn_abhd_purity=round(pur, 3),
                         knn_expected=round((lab.sum() - 1) / (len(lab) - 1), 3)))
    lt = pd.DataFrame(rows); lt.to_csv(C.DERIVED / "global_label_tests.csv", index=False)
    print(lt.to_string(index=False))

    # ---- clustering agreement ----
    cl_rows, clusters = [], pd.DataFrame({"accession": accs, "gene": [gene[a] for a in accs], "abhd": labels})
    ffl = uni.set_index("accession").funfam_ids.reindex(accs).fillna("").str.split(";").str[0]
    for k, M in mats.items():
        v = M.values.copy()
        if np.isnan(v).any():
            if k == "funfam":
                continue                                        # partial coverage: skip hierarchical clustering
            v = np.nan_to_num(v, nan=np.nanmax(v))
        v = (v + v.T) / 2; np.fill_diagonal(v, 0)
        Z = linkage(squareform(v, checks=False), method="average")
        for nc in (10, 20, 30):
            c = fcluster(Z, nc, "maxclust")
            if nc == 20: clusters[f"cl20_{k}"] = c
            ok = ffl != ""
            ct = pd.crosstab(c, labels)
            maj = [(c == cc).sum() for cc in ct.index if ct.loc[cc].get(True, 0) / ct.loc[cc].sum() >= 0.5 and ct.loc[cc].get(True, 0) > 0]
            cl_rows.append(dict(space=k, n_clusters=nc,
                                ARI_vs_funfam=round(adjusted_rand_score(ffl[ok], c[ok]), 3),
                                NMI_vs_funfam=round(normalized_mutual_info_score(ffl[ok], c[ok]), 3),
                                abhd_genes_in_majority_abhd_clusters=int(sum(ct.loc[cc].get(True, 0) for cc in ct.index
                                                                             if ct.loc[cc].get(True, 0) / ct.loc[cc].sum() >= 0.5)),
                                clusters_with_abhd=int((ct.get(True, pd.Series(dtype=int)) > 0).sum())))
    ca = pd.DataFrame(cl_rows); ca.to_csv(C.DERIVED / "global_clusters.csv", index=False)
    clusters.to_csv(C.DERIVED / "global_cluster_assignments.csv", index=False)
    print(ca.to_string(index=False))

    # ---- concordance between spaces ----
    names = list(mats); iu = np.triu_indices(len(accs), 1)
    conc = pd.DataFrame(index=names, columns=names, dtype=float)
    for a in names:
        for b in names:
            x, y = mats[a].values[iu], mats[b].values[iu]; ok = ~np.isnan(x) & ~np.isnan(y)
            conc.loc[a, b] = spearmanr(x[ok], y[ok])[0] if ok.sum() > 10 else np.nan
    conc.round(3).to_csv(C.DERIVED / "global_concordance.csv"); print(conc.round(2).to_string())

    # ---- discordant pairs: structurally close but sequence-distant (and vice versa) ----
    tm, sq, ff = mats["tmalign"].values[iu], mats["seqid"].values[iu], mats["funfam"].values[iu]
    ia, ib = iu
    df = pd.DataFrame(dict(a=[gene[accs[i]] for i in ia], b=[gene[accs[j]] for j in ib], tm_dist=tm, seq_dist=sq, funfam_dist=ff,
                           a_abhd=labels[ia], b_abhd=labels[ib]))
    df["tm_rank"] = df.tm_dist.rank(pct=True); df["seq_rank"] = df.seq_dist.rank(pct=True)
    df["kind"] = np.where((df.tm_rank <= 0.05) & (df.seq_rank >= 0.5), "structure_close_sequence_far",
                          np.where((df.seq_rank <= 0.05) & (df.tm_rank >= 0.5), "sequence_close_structure_far", ""))
    df[df.kind != ""].sort_values("tm_dist").to_csv(C.DERIVED / "discordant_pairs.csv", index=False)
    print(df.kind.value_counts().to_dict())

    # ---- figure: per-space kNN purity vs random expectation ----
    fig = go.Figure()
    fig.add_bar(x=lt.space, y=lt.knn_abhd_purity, name="observed (ABHD genes' 3 nearest neighbours that are ABHD)", marker_color="crimson")
    fig.add_bar(x=lt.space, y=lt.knn_expected, name="expected if the name carried no information", marker_color="lightgrey")
    fig.update_layout(barmode="group", template="simple_white", width=800, height=430,
                      title="Stage 6: does the ABHD name predict nearest neighbours?", yaxis_title="fraction")
    fig.write_html(C.FIGURES / "stage6_label_purity.html")


if __name__ == "__main__":
    main()
