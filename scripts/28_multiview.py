#!/usr/bin/env python3
"""
Multi-view agreement (docs/redesign_plan.md, P5, reduced): five views of relatedness, none treated as truth.
Views: structure (TM-align, CATH domain), structure (core K=200), local sequence, FunFam HMM, EC (hierarchical distance, complete EC only).
For every pair of views: Spearman correlation of pairwise distances over pairs defined in both, for ALL pairs and CROSS-FAMILY pairs (different paralog stem,
no shared FunFam); stem-cluster bootstrap CI (resample stems, de-duplicate genes, 500 draws). No Mantel p-values (pairs share paralog structure).
Per-gene discordance: EC3 of the nearest cross-family neighbour in each space.
Outputs: data/derived/multiview_association.csv, multiview_dependence_audit.csv, multiview_gene_discordance.csv
"""
import re
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

rng = np.random.default_rng(31)


def ec_list(x):
    out = []
    if isinstance(x, str):
        for e in x.split(";"):
            p = e.strip().split(".")
            if len(p) >= 3 and all(q not in ("-", "") for q in p[:3]): out.append(p[:4])
    return out


def ec_dist(a, b):
    if not a or not b: return np.nan
    best = 0
    for x in a:
        for y in b:
            d = 0
            for u, v in zip(x, y):
                if u != v or u.startswith("n"): break
                d += 1
            best = max(best, d)
    return 1 - best / 4


def family(gene):
    m = re.match(r"ABHD\d+", gene)
    return m.group(0) if m else re.match(r"[A-Z]+", gene).group(0)


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv"); accs = list(uni.accession); n = len(accs)
    ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
    mats = {k: pd.read_csv(C.DERIVED / f, index_col=0).reindex(index=accs, columns=accs).values for k, f in
            (("structure", "dist_tmalign.csv"), ("structure_core200", "dist_tmalign_core200.csv"), ("sequence", "dist_seqid.csv"), ("funfam", "dist_funfam.csv"))}
    ecs = [ec_list(ann.loc[a, "EC number"]) if a in ann.index else [] for a in accs]
    E = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(n):
            E[i, j] = 0.0 if i == j else ec_dist(ecs[i], ecs[j])
    mats["ec"] = E
    for k in ("funfam",):
        v = mats[k]; mats[k] = np.fmin(v, v.T) if np.isnan(v).any() else (v + v.T) / 2
    genes = uni.gene.to_numpy(); stem = np.array([family(g) for g in genes])
    ffs = [set(str(f).split(";")) - {"nan"} for f in uni.funfam_ids]
    cross = np.zeros((n, n), bool)
    for i in range(n):
        for j in range(n):
            cross[i, j] = stem[i] != stem[j] and not (ffs[i] & ffs[j])
    views = list(mats)
    iu = np.triu_indices(n, 1)

    def corr(a, b, mask_genes, only_cross):
        A, B = mats[a], mats[b]; ii, jj = iu
        sel = mask_genes[ii] & mask_genes[jj]
        if only_cross: sel &= cross[ii, jj]
        x, y = A[ii, jj][sel], B[ii, jj][sel]; ok = np.isfinite(x) & np.isfinite(y)
        if ok.sum() < 30: return np.nan, int(ok.sum())
        return spearmanr(x[ok], y[ok])[0], int(ok.sum())

    stems = np.unique(stem); idx_by_stem = {s: np.where(stem == s)[0] for s in stems}
    allg = np.ones(n, bool); rows = []
    boots = {}
    for _ in range(500):
        pick = rng.choice(stems, len(stems)); m = np.zeros(n, bool)
        for s in pick: m[idx_by_stem[s]] = True
        boots.setdefault("masks", []).append(m)
    for a, b in combinations(views, 2):
        for scope in ("all_pairs", "cross_family_pairs"):
            r0, npairs = corr(a, b, allg, scope == "cross_family_pairs")
            bs = [corr(a, b, m, scope == "cross_family_pairs")[0] for m in boots["masks"]]
            bs = [x for x in bs if np.isfinite(x)]
            rows.append(dict(view_a=a, view_b=b, scope=scope, n_pairs=npairs, spearman=round(r0, 3) if np.isfinite(r0) else np.nan,
                             ci_lo=round(np.percentile(bs, 2.5), 3) if bs else np.nan, ci_hi=round(np.percentile(bs, 97.5), 3) if bs else np.nan))
    A = pd.DataFrame(rows); A.to_csv(C.DERIVED / "multiview_association.csv", index=False)
    audit = pd.DataFrame([
        ("sequence", "funfam", "dependent by construction: FunFam is built from sequence similarity"),
        ("structure", "structure_core200", "dependent: same models; core is centred on the elbow"),
        ("ec", "funfam", "partly dependent: EC assigned by sequence similarity; FunFam functional families use EC-annotated members"),
        ("ec", "sequence", "partly dependent: EC propagated by sequence similarity"),
        ("structure", "sequence", "independent measurement, shared evolutionary signal"),
    ], columns=["view_a", "view_b", "note"]); audit.to_csv(C.DERIVED / "multiview_dependence_audit.csv", index=False)
    pd.set_option("display.width", 200)
    print(A.pivot_table(index=["view_a", "view_b"], columns="scope", values=["spearman", "n_pairs"]).round(3).to_string())
    print(A[A.scope == "cross_family_pairs"][["view_a", "view_b", "spearman", "ci_lo", "ci_hi", "n_pairs"]].to_string(index=False))
    # per-gene discordance of EC3 of nearest cross-family neighbour
    ec3 = [".".join(e[:3]) if e else None for e in [x[0] if x else None for x in ecs]]
    ec3s = [set(".".join(e[:3]) for e in x) for x in ecs]
    out = []
    for i in range(n):
        if not ec3s[i]: continue
        rec = dict(gene=genes[i], abhd_named=bool(uni.abhd_named[i]), ec3=";".join(sorted(ec3s[i])))
        for v in ("structure", "structure_core200", "sequence", "funfam"):
            d = mats[v][i].copy(); ok = np.where(cross[i] & np.isfinite(d) & np.array([bool(s) for s in ec3s]))[0]
            if len(ok):
                j = ok[np.argmin(d[ok])]; rec[f"nn_{v}"] = genes[j]; rec[f"agree_{v}"] = bool(ec3s[i] & ec3s[j])
        out.append(rec)
    G = pd.DataFrame(out); G["struct_only"] = G.agree_structure & ~G.agree_sequence; G["seq_only"] = ~G.agree_structure & G.agree_sequence
    G.to_csv(C.DERIVED / "multiview_gene_discordance.csv", index=False)
    print("cross-family nearest-neighbour EC3 agreement:", {v: round(float(G[f"agree_{v}"].mean()), 3) for v in ("structure", "structure_core200", "sequence", "funfam")}, "n =", len(G))
    print("structure agrees, sequence not:", G[G.struct_only].gene.tolist()); print("sequence agrees, structure not:", G[G.seq_only].gene.tolist())
    print("all-pair structure-sequence Spearman (check vs global_concordance 0.52):", A[(A.view_a == "structure") & (A.view_b == "sequence") & (A.scope == "all_pairs")].spearman.iat[0])


if __name__ == "__main__":
    main()
