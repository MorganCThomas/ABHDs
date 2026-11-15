#!/usr/bin/env python3
"""
Function transfer as a decay curve over the maximum allowed sequence identity (docs/redesign_plan.md, P2).

Identity: global Needleman-Wunsch (BLOSUM62, -11/-1) on trimmed domains; identical columns / aligned columns; coverage = aligned / shorter length.
Spaces (prespecified): tmalign, seqid, funfam, anchor_elbow. Modes: all neighbours, family-excluded (shared FunFam or paralog stem).
Statistics: difference to the majority-class baseline on the same queries, stem-cluster bootstrap CI (resample stems, 2000 draws).
Outputs: data/derived/identity_pairs.csv, function_decay.csv, function_decay_strata.csv, function_by_nn_identity.csv; figures/function_decay.html
"""
import re
import sys
from collections import Counter
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from Bio.Align import PairwiseAligner, substitution_matrices

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.structops import load

SPACES = ["tmalign", "seqid", "funfam", "anchor_elbow"]
CUTS = [100, 40, 35, 30, 25, 20]
MIN_POOL = 5
rng = np.random.default_rng(11)
_AL = None


def _init():
    global _AL
    _AL = PairwiseAligner(mode="global", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-11, extend_gap_score=-1)


def _ident(args):
    a, b, sa, sb = args
    aln = _AL.align(sa, sb)[0]
    same = tot = 0
    for (s0, e0), (s1, e1) in zip(*aln.aligned):
        tot += e0 - s0
        same += sum(x == y for x, y in zip(sa[s0:e0], sb[s1:e1]))
    return a, b, same / tot if tot else 0.0, tot / min(len(sa), len(sb))


def identity_matrix(uni):
    seqs = {a: load(a)[2].replace("X", "A") for a in uni.accession}
    accs = list(seqs)
    jobs = [(a, b, seqs[a], seqs[b]) for i, a in enumerate(accs) for b in accs[i + 1:]]
    with Pool(initializer=_init) as p:
        res = p.map(_ident, jobs, chunksize=40)
    I = pd.DataFrame(np.eye(len(accs)), index=accs, columns=accs); Cv = I.copy()
    rows = []
    for a, b, i, c in res:
        I.loc[a, b] = I.loc[b, a] = i; Cv.loc[a, b] = Cv.loc[b, a] = c; rows.append((a, b, i, c))
    pd.DataFrame(rows, columns=["a", "b", "identity", "coverage"]).to_csv(C.DERIVED / "identity_pairs.csv", index=False)
    return I


def ec_set(x, level=3):
    if not isinstance(x, str) or not x.strip(): return set()
    out = set()
    for e in x.split(";"):
        parts = e.strip().split(".")[:level]
        if e.strip() and len(parts) == level and "-" not in parts: out.add(".".join(parts))
    return out


def family(gene):
    m = re.match(r"ABHD\d+", gene)
    return m.group(0) if m else re.match(r"[A-Z]+", gene).group(0)


def wilson(k, n, z=1.96):
    if n == 0: return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return (c - h) / d, (c + h) / d


def run(I, lab, fam, ffs, maj, space_mats, q_ok, n_ok, cut, mode, space):
    """return per-query dataframe rows: acc, stem, correct, base_maj, base_rand for queries meeting the pool size"""
    M = space_mats[space]; out = []
    for a in lab:
        if a not in M.index or a not in q_ok: continue
        d = M.loc[a].drop(a)
        keep = [x for x in d.index if x in lab and x in n_ok and np.isfinite(d[x]) and I.loc[a, x] <= cut / 100.0 + 1e-9
                and (mode == "all" or (fam[x] != fam[a] and not (ffs[x] & ffs[a])))]
        if len(keep) < MIN_POOL: continue
        dd = d[keep].sort_values()
        nn = dd.index[0]
        out.append(dict(acc=a, stem=fam[a], correct=bool(lab[nn] & lab[a]), base_maj=bool(maj in lab[a]),
                        base_rand=float(np.mean([bool(lab[x] & lab[a]) for x in keep])), nn=nn, nn_ident=I.loc[a, nn]))
    return pd.DataFrame(out)


def boot_diff(df, nb=2000):
    stems = df.stem.unique(); groups = {s: g for s, g in df.groupby("stem")}
    d = df.correct.mean() - df.base_maj.mean()
    bs = []
    for _ in range(nb):
        pick = rng.choice(stems, len(stems))
        g = pd.concat([groups[s] for s in pick])
        bs.append(g.correct.mean() - g.base_maj.mean())
    lo, hi = np.percentile(bs, [2.5, 97.5])
    return d, lo, hi


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv").set_index("accession")
    prov = pd.read_csv(C.DERIVED / "label_provenance.csv").set_index("accession")
    ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
    I = identity_matrix(uni.reset_index())
    off = I.values[np.triu_indices(len(I), 1)]
    print("identity: median %.2f, p95 %.2f, share > 30%% %.3f, > 40%% %.3f" % (np.median(off), np.percentile(off, 95), (off > .3).mean(), (off > .4).mean()))
    lab = {a: ec_set(ann.loc[a, "EC number"]) for a in uni.index if a in ann.index}; lab = {a: s for a, s in lab.items() if s}
    fam = {a: family(uni.loc[a, "gene"]) for a in uni.index}
    ffs = {a: set(str(f).split(";")) - {"nan"} for a, f in uni.funfam_ids.items()}
    maj = Counter(c for s in lab.values() for c in s).most_common(1)[0][0]
    mats = {s: pd.read_csv(C.DERIVED / f"dist_{s}.csv", index_col=0) for s in SPACES}
    allq = set(lab); exp_any = set(prov.index[prov.evidence_any_experimental & prov.index.isin(lab)]); exp_all = set(prov.index[prov.evidence_all_experimental & prov.index.isin(lab)])

    rows, rows_s, rows_b = [], [], []
    for space in SPACES:
        for mode in ("all", "leave_family_out"):
            for cut in CUTS:
                df = run(I, lab, fam, ffs, maj, mats, allq, allq, cut, mode, space)
                if len(df) < 10: continue
                d, lo, hi = boot_diff(df)
                k = int(df.correct.sum()); wl, wh = wilson(k, len(df))
                rows.append(dict(space=space, mode=mode, max_identity=cut, n_queries=len(df), mean_pool_nn_identity=round(df.nn_ident.mean(), 3), accuracy=round(df.correct.mean(), 3),
                                 wilson_lo=round(wl, 3), wilson_hi=round(wh, 3), majority_baseline=round(df.base_maj.mean(), 3), random_baseline=round(df.base_rand.mean(), 3),
                                 diff_vs_majority=round(d, 3), diff_ci_lo=round(lo, 3), diff_ci_hi=round(hi, 3)))
    F = pd.DataFrame(rows); F.to_csv(C.DERIVED / "function_decay.csv", index=False)
    # strata (descriptive), cutoffs 100 and 30, both modes
    for name, qset, nset in (("all_complete_ec3", allq, allq), ("query_any_experimental", exp_any, allq), ("query_and_neighbour_all_experimental", exp_all, exp_all)):
        for space in SPACES:
            for mode in ("all", "leave_family_out"):
                for cut in (100, 30):
                    df = run(I, lab, fam, ffs, maj, mats, qset, nset, cut, mode, space)
                    if len(df) < 10: continue
                    wl, wh = wilson(int(df.correct.sum()), len(df))
                    rows_s.append(dict(stratum=name, space=space, mode=mode, max_identity=cut, n_queries=len(df), accuracy=round(df.correct.mean(), 3),
                                       wilson_lo=round(wl, 3), wilson_hi=round(wh, 3), majority_baseline=round(df.base_maj.mean(), 3)))
    pd.DataFrame(rows_s).to_csv(C.DERIVED / "function_decay_strata.csv", index=False)
    # accuracy by identity of the actual nearest labelled neighbour
    bins = [(0, .25), (.25, .30), (.30, .40), (.40, 1.01)]
    for space in SPACES:
        df = run(I, lab, fam, ffs, maj, mats, allq, allq, 100, "all", space)
        for lo, hi in bins:
            s = df[(df.nn_ident >= lo) & (df.nn_ident < hi)]
            if len(s): rows_b.append(dict(space=space, nn_identity_bin=f"{lo:.2f}-{min(hi,1):.2f}", n=len(s), accuracy=round(s.correct.mean(), 3)))
    pd.DataFrame(rows_b).to_csv(C.DERIVED / "function_by_nn_identity.csv", index=False)
    pd.set_option("display.width", 220)
    print(F[F["mode"] == "all"].query("max_identity in [100, 30, 25, 20]").to_string(index=False))
    print(F[F["mode"] == "leave_family_out"].to_string(index=False))
    print(pd.DataFrame(rows_b).to_string(index=False))
    # acceptance checks against script 14 output
    L = pd.read_csv(C.DERIVED / "function_loo.csv"); L = L[L.ec_level == 3]
    for space in ("tmalign", "seqid", "funfam", "anchor_elbow"):
        for mode in ("all", "leave_family_out"):
            ref = L[(L.space == space) & (L["mode"] == mode)].acc_1nn.iat[0]
            got = F[(F.space == space) & (F["mode"] == mode) & (F.max_identity == 100)].accuracy.iat[0]
            print(f"check {space:13s} {mode:17s} script14 {ref:.3f}  here {got:.3f}")
    fig = go.Figure()
    for space in SPACES:
        s = F[(F.space == space) & (F["mode"] == "all")]
        fig.add_scatter(x=s.max_identity, y=s.diff_vs_majority, mode="lines+markers", name=space, error_y=dict(type="data", symmetric=False, array=s.diff_ci_hi - s.diff_vs_majority, arrayminus=s.diff_vs_majority - s.diff_ci_lo))
    fig.add_hline(y=0, line_dash="dot"); fig.update_xaxes(autorange="reversed", title="maximum neighbour identity (%)"); fig.update_yaxes(title="accuracy minus majority-class baseline")
    fig.update_layout(template="simple_white", title="Function transfer vs allowed neighbour identity (stem-cluster bootstrap 95% CI)", width=900, height=500)
    fig.write_html(C.FIGURES / "function_decay.html")


if __name__ == "__main__":
    main()
