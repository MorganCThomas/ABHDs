#!/usr/bin/env python3
"""
Study stage 7 - can annotation be transferred between neighbours? Leave-one-out test, per space.

Label = EC class from UniProt (set per gene; levels 1-3 and 4). Prediction for gene i = label of its nearest annotated
other gene(s) in a distance space (k=1, and k=3 majority). Correct if any predicted label is in gene i's label set.
Incomplete ECs (3.-.-) are not classes. Modes: all neighbours, and leave-family-out (neighbours sharing a FunFam or paralog
stem are excluded). Baselines: majority class and random labelled neighbour. Wilson 95% CI on 1-NN accuracy.
Also writes hypotheses for genes without an EC annotation, from the best-performing structural space.
Outputs: data/derived/function_loo.csv, function_hypotheses.csv; figures/stage7_loo.html
"""
import re
import sys
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

SPACES = ["tmalign", "tmshort", "seqid", "protbert", "funfam", "soap", "triad", "anchor_elbow", "anchor_cat"]


def ec_set(x, level):
    """EC classes at `level`; incomplete ECs ('3.-.-', '3.1.-') are NOT classes and are dropped at levels they do not define."""
    if not isinstance(x, str) or not x.strip():
        return set()
    out = set()
    for e in x.split(";"):
        parts = e.strip().split(".")[:level]
        if e.strip() and len(parts) == level and "-" not in parts:
            out.add(".".join(parts))
    return out


def family(gene):
    """Paralog group used for leave-family-out: ABHD stem (ABHD17A -> ABHD17) or the letter prefix (CES1 -> CES)."""
    m = re.match(r"ABHD\d+", gene)
    return m.group(0) if m else re.match(r"[A-Z]+", gene).group(0)


def wilson(k, n, z=1.96):
    if n == 0:
        return (np.nan, np.nan)
    p = k / n; d = 1 + z * z / n; c = p + z * z / (2 * n); h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n))
    return ((c - h) / d, (c + h) / d)


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv").set_index("accession")
    ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
    uni["ec"] = ann["EC number"].reindex(uni.index)
    n_ec = uni.ec.notna().sum()
    print(f"genes with an EC number: {n_ec}/{len(uni)}")
    rows = []
    fam = {a: family(g) for a, g in uni.gene.items()}
    ffs = {a: set(str(f).split(";")) - {"nan"} for a, f in uni.funfam_ids.items()}
    for level in (1, 2, 3):
        lab = {a: ec_set(e, level) for a, e in uni.ec.items()}
        lab = {a: s_ for a, s_ in lab.items() if s_}
        allc = Counter(c for s_ in lab.values() for c in s_)
        maj = allc.most_common(1)[0][0]
        base = np.mean([maj in s_ for s_ in lab.values()])
        # random-neighbour baseline: probability that a random other labelled gene shares a label
        rnd = np.mean([np.mean([bool(lab[a] & lab[b]) for b in lab if b != a]) for a in lab])
        for sp in SPACES:
            M = pd.read_csv(C.DERIVED / f"dist_{sp}.csv", index_col=0)
            for mode in ("all", "leave_family_out"):
                ok1, ok3, n = [], [], 0
                for a in lab:
                    if a not in M.index:
                        continue
                    d = M.loc[a].drop(a)
                    keep = [x for x in d.index if x in lab and np.isfinite(d[x])
                            and (mode == "all" or (fam[x] != fam[a] and not (ffs[x] & ffs[a])))]
                    d = d[keep].sort_values()
                    if len(d) < 3:
                        continue
                    n += 1
                    ok1.append(bool(lab[d.index[0]] & lab[a]))
                    votes = Counter(c for x in d.index[:3] for c in lab[x])
                    ok3.append(votes.most_common(1)[0][0] in lab[a])
                lo, hi = wilson(sum(ok1), n)
                rows.append(dict(ec_level=level, space=sp, mode=mode, n_genes=n, acc_1nn=round(np.mean(ok1), 3),
                                 ci95_lo=round(lo, 3), ci95_hi=round(hi, 3), acc_3nn_vote=round(np.mean(ok3), 3),
                                 majority_baseline=round(base, 3), random_neighbour_baseline=round(rnd, 3), majority_label=maj))
    df = pd.DataFrame(rows); df.to_csv(C.DERIVED / "function_loo.csv", index=False)
    print(df.to_string(index=False))

    # hypotheses for genes without EC: majority EC3 among the 3 nearest annotated genes in the TM-align and SOAP-triad consensus
    lab3 = {a: ec_set(e, 3) for a, e in uni.ec.items()}; lab3 = {a: s_ for a, s_ in lab3.items() if s_}
    Mt = pd.read_csv(C.DERIVED / "dist_tmalign.csv", index_col=0)
    out = []
    for a in uni.index[~uni.index.isin(lab3)]:
        d = Mt.loc[a].drop(a); d = d[[x for x in d.index if x in lab3]].sort_values()[:3]
        votes = Counter(c for x in d.index for c in lab3[x])
        pred, cnt = votes.most_common(1)[0]
        out.append(dict(gene=uni.loc[a, "gene"], accession=a, abhd_named=uni.loc[a, "abhd_named"], predicted_ec3=pred,
                        support=f"{cnt}/3", neighbours=", ".join(f"{uni.loc[x, 'gene']}({d[x]:.2f})" for x in d.index)))
    h = pd.DataFrame(out)
    cat = pd.read_csv(C.DERIVED / "catalytic_sites.csv").set_index("accession").category
    h["catalytic_category"] = h.accession.map(cat)
    # an EC cannot be transferred to a protein whose elbow carries no nucleophile (pseudoenzyme-like): flag, do not delete
    h["transfer_applicable"] = ~h.catalytic_category.fillna("").str.startswith("no nucleophile")
    h.to_csv(C.DERIVED / "function_hypotheses.csv", index=False)
    print(f"{len(out)} genes without EC; hypotheses written")

    fig = go.Figure()
    d3 = df[(df.ec_level == 3) & (df["mode"] == "leave_family_out")]
    fig.add_bar(x=d3.space, y=d3.acc_1nn, name="nearest annotated neighbour", marker_color="steelblue")
    fig.add_bar(x=d3.space, y=d3.acc_3nn_vote, name="vote of 3 nearest", marker_color="darkorange")
    fig.add_hline(y=d3.majority_baseline.iloc[0], line_dash="dot", annotation_text="majority-class baseline")
    fig.update_layout(barmode="group", template="simple_white", width=800, height=430, yaxis_title="leave-one-out accuracy",
                      title="Stage 7: EC level-3 transfer, paralog families and shared FunFams excluded")
    fig.write_html(C.FIGURES / "stage7_loo.html")


if __name__ == "__main__":
    main()
