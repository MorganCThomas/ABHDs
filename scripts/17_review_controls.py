#!/usr/bin/env python3
"""
Controls added after the independent review (reports/09_review_revisions.md).

 A. Elbow calibration: held-out split. Phi/psi windows are re-derived (observed range +- margin) from a random half of the
    31 gold nucleophiles and scored on the other half; and elbow specificity (S/C/T/D residues passing per gene).
 B. ABHD-name tests under two label definitions (strict gene symbol; broad UniProt name 'abhydrolase domain-containing'),
    with gene-bootstrap CIs, exclusion of near neighbours, a length-stratified permutation, and a cluster-level null
    (random gene sets of the same size) across linkages.
Outputs: data/derived/review_elbow_heldout.csv, review_label_robustness.csv, review_cluster_null.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import squareform

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.catalytic import apply_elbow

rng = np.random.default_rng(1)


def heldout_elbow():
    et = pd.read_csv(C.DERIVED / "elbow_table.csv")
    g = pd.read_csv(C.DERIVED / "benchmark_stage3_genes.csv")
    gold = []
    for r in g.itertuples():
        row = et[(et.accession == r.accession) & (et.resnum == r.gold_nuc)]
        if len(row):
            gold.append(row.iloc[0])
    gold = pd.DataFrame(gold)
    out = []
    for rep in range(200):
        idx = rng.permutation(len(gold)); tr, te = gold.iloc[idx[:len(gold) // 2]], gold.iloc[idx[len(gold) // 2:]]
        m = 10
        phi = (tr.phi.min() - m, tr.phi.max() + m); psi = (tr.psi.min() - m, tr.psi.max() + m)
        hit = apply_elbow(te, phi=phi, psi=psi)
        out.append(dict(rep=rep, recall_heldout=hit.mean(), phi_lo=phi[0], phi_hi=phi[1], psi_lo=psi[0], psi_hi=psi[1]))
    ho = pd.DataFrame(out)
    # specificity: residues passing the default elbow definition, per gene and by residue type
    ok = et[apply_elbow(et)]
    ok_nuc = ok[ok.aa.isin(list("SCTD"))]
    per_gene = ok_nuc.groupby("accession").size()
    summ = dict(n_gold=len(gold), gold_phi_range=f"{gold.phi.min():.0f}..{gold.phi.max():.0f}",
                gold_psi_range=f"{gold.psi.min():.0f}..{gold.psi.max():.0f}",
                heldout_recall_mean=round(ho.recall_heldout.mean(), 3), heldout_recall_min=round(ho.recall_heldout.min(), 3),
                n_genes_with_SCTD_elbow=int(per_gene.size), n_genes_exactly_one=int((per_gene == 1).sum()),
                passing_residues_total=len(ok), passing_residues_SCTD=len(ok_nuc),
                passing_by_aa=";".join(f"{a}:{n}" for a, n in ok.aa.value_counts().items()))
    pd.DataFrame([summ]).to_csv(C.DERIVED / "review_elbow_heldout.csv", index=False)
    print(pd.DataFrame([summ]).T.to_string())


def stem(g):
    import re
    return re.sub(r"[A-Z]$", "", g) if re.match(r"ABHD\d+[A-Z]$", g) else g


def label_robustness():
    tm = pd.read_csv(C.DERIVED / "dist_tmalign.csv", index_col=0).values
    iu_ = np.triu_indices(len(tm), 1); QCUT = float((tm[iu_] < 0.2).mean())     # fraction of pairs the TM < 0.2 cut removes
    uni = pd.read_csv(C.DERIVED / "universe.csv").set_index("accession")
    uni["dom_len"] = pd.read_csv(C.DERIVED / "structure_qc.csv").set_index("accession").dom_len
    rows = []
    for defn in ("abhd_named", "abhd_family"):
        for sp in ("tmalign", "tmshort", "seqid", "funfam", "anchor_elbow"):
            M = pd.read_csv(C.DERIVED / f"dist_{sp}.csv", index_col=0)
            M = M.loc[M.notna().any(axis=1), M.notna().any(axis=0)]
            accs = list(M.index); ab = uni.loc[accs, defn].to_numpy(bool); genes = uni.loc[accs, "gene"].to_numpy()
            L = uni.loc[accs, "dom_len"].to_numpy(float)
            v = M.values.copy(); np.fill_diagonal(v, np.inf); v[np.isnan(v)] = np.inf
            R = int(ab.sum()) - 1; base = R / (len(accs) - 1)
            vals_ = M.values[np.triu_indices(len(M), 1)]; vals_ = vals_[np.isfinite(vals_)]
            qthr = float(np.quantile(vals_, QCUT))            # quantile-matched: drop the same share of closest pairs in every space
            ia = np.where(ab)[0]

            def per_gene(i, excl_near=None, lab=ab):
                order = np.argsort(v[i]); order = order[np.isfinite(v[i, order])]
                order = np.array([j for j in order if not (lab[j] and stem(genes[j]) == stem(genes[i]))]) if excl_near == "paralog" else order
                if excl_near is not None and excl_near != "paralog":
                    order = order[v[i, order] >= excl_near]
                return (float(lab[order[0]]) if len(order) else np.nan,
                        float(lab[order[:R]].mean()) if len(order) else np.nan)
            for variant, ex in (("all", None), ("skip_paralog_stem", "paralog"), ("skip_closest_3.4pct_pairs", qthr)):
                vals = np.array([per_gene(i, ex) for i in ia]); vals = vals[~np.isnan(vals[:, 0])]
                boots = [vals[rng.integers(0, len(vals), len(vals))].mean(0) for _ in range(2000)]
                lo, hi = np.percentile(boots, [2.5, 97.5], axis=0)
                rows.append(dict(definition=defn, space=sp, variant=variant, n_abhd=int(ab.sum()), nn1=round(vals[:, 0].mean(), 3),
                                 nn1_lo=round(lo[0], 3), nn1_hi=round(hi[0], 3), r_precision=round(vals[:, 1].mean(), 3),
                                 rp_lo=round(lo[1], 3), rp_hi=round(hi[1], 3), chance=round(base, 3)))
            # length-stratified permutation of R-precision (labels shuffled within length quartiles)
            q = pd.qcut(L, 4, labels=False, duplicates="drop")
            obs = np.mean([per_gene(i)[1] for i in ia])
            null = []
            for _ in range(1000):
                lab = ab.copy()
                for k in np.unique(q):
                    m = np.where(q == k)[0]; lab[m] = rng.permutation(ab[m])
                Rn = int(lab.sum()) - 1
                r = []
                for i in np.where(lab)[0]:
                    order = np.argsort(v[i]); order = order[np.isfinite(v[i, order])][:Rn]; r.append(lab[order].mean())
                null.append(np.mean(r))
            rows.append(dict(definition=defn, space=sp, variant="length_stratified_perm", n_abhd=int(ab.sum()), r_precision=round(obs, 3),
                             chance=round(float(np.mean(null)), 3), rp_lo=round(float(np.percentile(null, 2.5)), 3),
                             rp_hi=round(float(np.percentile(null, 97.5)), 3),
                             nn1=round(float((np.sum(np.array(null) >= obs) + 1) / 1001), 4)))   # nn1 column holds the perm p here
    df = pd.DataFrame(rows); df.to_csv(C.DERIVED / "review_label_robustness.csv", index=False)
    print(df.to_string(index=False))


def cluster_null():
    uni = pd.read_csv(C.DERIVED / "universe.csv").set_index("accession")
    rows = []
    for sp in ("tmalign", "tmshort"):
        M = pd.read_csv(C.DERIVED / f"dist_{sp}.csv", index_col=0); accs = list(M.index)
        v = (M.values + M.values.T) / 2; np.fill_diagonal(v, 0)
        for defn in ("abhd_named", "abhd_family"):
            ab = uni.loc[accs, defn].to_numpy(bool)
            for method in ("average", "complete", "single", "ward"):
                if method == "ward":
                    continue
                Z = linkage(squareform(v, checks=False), method=method)
                for nc in (10, 20, 30):
                    c = fcluster(Z, nc, "maxclust")
                    obs = len(set(c[ab]))
                    null = [len(set(c[rng.permutation(len(c))[:ab.sum()]])) for _ in range(2000)]
                    rows.append(dict(space=sp, definition=defn, linkage=method, n_clusters=nc, abhd_clusters_observed=obs,
                                     random_mean=round(float(np.mean(null)), 2), p_le_observed=round(float((np.sum(np.array(null) <= obs) + 1) / 2001), 4)))
    df = pd.DataFrame(rows); df.to_csv(C.DERIVED / "review_cluster_null.csv", index=False)
    print(df.to_string(index=False))


if __name__ == "__main__":
    heldout_elbow(); label_robustness(); cluster_null()
