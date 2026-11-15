#!/usr/bin/env python3
"""
Structure-derived, size-controlled cores instead of CATH ranges (docs/redesign_plan.md, P4).

Core(K) = the K confident (pLDDT >= 70) residues of the FULL AlphaFold model whose C-alpha are nearest the elbow C-alpha (3D); original order kept.
TM-align (two normalisations) and local sequence distance on cores for K in {150, 180, 200}; the ABHD nearest-neighbour tests are repeated and
compared with the CATH-domain results by a paired gene bootstrap with an equivalence margin of 0.10.
NOTE: cores are centred on the elbow, so they are not independent of the elbow-based measures.
Outputs: data/derived/core_residues.csv, dist_tmalign_core<K>.csv, dist_tmshort_core<K>.csv, dist_seqid_core<K>.csv, core_vs_cath_label_tests.csv
"""
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from scipy.stats import spearmanr
from tmtools import tm_align

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.elbow import THREE_TO_ONE

KS = (150, 180, 200)
PRIMARY = 200
MARGIN = 0.10
rng = np.random.default_rng(21)
_AL = None


def _init():
    global _AL
    _AL = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-11, extend_gap_score=-1)


def core(acc, elbow, K):
    s = PDBParser(QUIET=True).get_structure(acc, str(C.STRUCTURES / f"{acc}.pdb"))
    ids, ca, seq = [], [], []
    for r in s.get_residues():
        if "CA" in r and r["CA"].get_bfactor() >= 70:
            ids.append(r.get_id()[1]); ca.append(r["CA"].coord); seq.append(THREE_TO_ONE.get(r.get_resname(), "X"))
    ids = np.array(ids); ca = np.array(ca)
    e = ca[np.where(ids == elbow)[0][0]]
    order = np.argsort(np.linalg.norm(ca - e, axis=1))[:K]; order.sort()
    return ids[order], ca[order], "".join(seq[i] for i in order)


def _tm(a):
    (i, ca1, s1), (j, ca2, s2) = a
    r = tm_align(ca1, ca2, s1, s2)
    return i, j, (r.tm_norm_chain1 + r.tm_norm_chain2) / 2, max(r.tm_norm_chain1, r.tm_norm_chain2)


def _sw(a):
    i, j, s1, s2 = a
    return i, j, _AL.score(s1, s2)


def matrices(data, accs):
    n = len(accs); T = np.eye(n); Ts = np.eye(n)
    jobs = [((i, data[accs[i]][1], data[accs[i]][2]), (j, data[accs[j]][1], data[accs[j]][2])) for i in range(n) for j in range(i + 1, n)]
    with Pool() as p:
        for i, j, t, ts in p.imap_unordered(_tm, jobs, chunksize=50):
            T[i, j] = T[j, i] = t; Ts[i, j] = Ts[j, i] = ts
    seqs = [data[a][2].replace("X", "A") for a in accs]
    with Pool(initializer=_init) as p:
        res = p.map(_sw, [(i, j, seqs[i], seqs[j]) for i in range(n) for j in range(i, n)], chunksize=100)
    S = np.zeros((n, n))
    for i, j, sc in res: S[i, j] = S[j, i] = sc
    d = np.sqrt(np.outer(np.diag(S), np.diag(S))); D = 1 - S / d; np.fill_diagonal(D, 0)
    f = lambda M: pd.DataFrame(M, index=accs, columns=accs)
    return f(1 - T), f(1 - Ts), f(D)


def nn_stats(M, ab, genes, R, quantile_cut=None):
    v = M.values.copy().astype(float); np.fill_diagonal(v, np.inf)
    if quantile_cut is not None: v[v < quantile_cut] = np.inf
    out = []
    for i in np.where(ab)[0]:
        order = np.argsort(v[i]); order = order[np.isfinite(v[i, order])]
        out.append((float(ab[order[0]]), float(ab[order[:R]].mean())))
    return np.array(out)


def main():
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv").set_index("accession")
    uni = pd.read_csv(C.DERIVED / "universe.csv"); qc = pd.read_csv(C.DERIVED / "structure_qc.csv").set_index("accession")
    uni = uni[uni.accession.isin(cs.index[cs.elbow_resnum.notna()])].reset_index(drop=True)
    accs = list(uni.accession); ab = uni.abhd_named.to_numpy(bool); genes = uni.gene.to_numpy()
    cath = {k: pd.read_csv(C.DERIVED / f"dist_{k}.csv", index_col=0).reindex(index=accs, columns=accs) for k in ("tmalign", "tmshort", "seqid")}
    # core residues and descriptors
    rows, results = [], []
    for K in KS:
        data = {a: core(a, int(cs.loc[a, "elbow_resnum"]), K) for a in accs}
        for a in accs:
            ids = data[a][0]; lo, hi = int(qc.loc[a, "dom_start"]), int(qc.loc[a, "dom_end"])
            rows.append(dict(K=K, gene=uni.set_index("accession").loc[a, "gene"], accession=a, n_core=len(ids), frac_outside_cath=round(float(((ids < lo) | (ids > hi)).mean()), 3),
                             n_segments=int(1 + (np.diff(ids) > 1).sum()), core_lt_K=len(ids) < K))
        T, Ts, Sq = matrices(data, accs)
        for name, M in (("tmalign", T), ("tmshort", Ts), ("seqid", Sq)):
            M.to_csv(C.DERIVED / f"dist_{name}_core{K}.csv")
            iu = np.triu_indices(len(accs), 1)
            rho = spearmanr(M.values[iu], cath[name].values[iu])[0]
            R = int(ab.sum()) - 1
            s_core = nn_stats(M, ab, genes, R)
            s_cath = nn_stats(cath[name], ab, genes, R)
            # quantile-matched removal of the closest 3.4% pairs
            q = 0.034
            qc_core = nn_stats(M, ab, genes, R, np.quantile(M.values[iu], q)); qc_cath = nn_stats(cath[name], ab, genes, R, np.quantile(cath[name].values[iu], q))
            for variant, a1, a0 in (("all", s_core, s_cath), ("closest_3.4pct_removed", qc_core, qc_cath)):
                dif = a1 - a0                                   # per ABHD gene: core minus CATH
                bs = np.array([dif[rng.integers(0, len(dif), len(dif))].mean(0) for _ in range(2000)])
                lo, hi = np.percentile(bs, [2.5, 97.5], axis=0)
                results.append(dict(K=K, space=name, variant=variant, spearman_core_vs_cath=round(rho, 3), nn1_core=round(a1[:, 0].mean(), 3), nn1_cath=round(a0[:, 0].mean(), 3),
                                    nn1_diff=round(dif[:, 0].mean(), 3), nn1_ci_lo=round(lo[0], 3), nn1_ci_hi=round(hi[0], 3), nn1_equivalent_within_0p10=bool(lo[0] > -MARGIN and hi[0] < MARGIN),
                                    rprec_core=round(a1[:, 1].mean(), 3), rprec_cath=round(a0[:, 1].mean(), 3), rprec_diff=round(dif[:, 1].mean(), 3),
                                    rprec_ci_lo=round(lo[1], 3), rprec_ci_hi=round(hi[1], 3), rprec_equivalent_within_0p10=bool(lo[1] > -MARGIN and hi[1] < MARGIN)))
        print(f"K={K} done", flush=True)
    R_ = pd.DataFrame(rows); R_.to_csv(C.DERIVED / "core_residues.csv", index=False)
    T_ = pd.DataFrame(results); T_.to_csv(C.DERIVED / "core_vs_cath_label_tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(R_.groupby("K").agg(median_core=("n_core", "median"), lt_K=("core_lt_K", "sum"), outside=("frac_outside_cath", "median"), segs=("n_segments", "median")).to_string())
    print(T_.to_string(index=False))


if __name__ == "__main__":
    main()
