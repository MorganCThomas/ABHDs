#!/usr/bin/env python3
"""
Study stage 6b - stricter tests of "the ABHD name marks a clade", designed to separate paralog-group coherence from
superfamily-wide coherence.

For each space and each ABHD-named gene:
  nn1_abhd      is its single nearest other gene ABHD-named?
  r_precision   among its R nearest genes (R = number of other ABHD genes), the fraction that are ABHD-named
  nn_nonabhd_first  is the nearest non-ABHD gene closer than the nearest ABHD gene?
Also: FunFam purity (of the FunFams that contain an ABHD gene, how many contain only ABHD genes) and, per ABHD gene,
whether a paralog (same name stem, e.g. ABHD17A/B/C) is its nearest neighbour.
Outputs: data/derived/abhd_label_summary.csv, abhd_label_by_gene.csv, funfam_purity.csv; figures/stage6b_*.html
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

SPACES = ["tmalign", "tmshort", "seqid", "protbert", "funfam", "soap", "triad", "anchor_elbow", "anchor_cat"]


def stem(g):                      # ABHD17A -> ABHD17, ABHD14B -> ABHD14
    return re.sub(r"[A-Z]$", "", g) if re.match(r"ABHD\d+[A-Z]$", g) else g


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv").set_index("accession")
    qc = pd.read_csv(C.DERIVED / "structure_qc.csv").set_index("accession")
    uni["dom_len"] = qc.dom_len
    print("domain length: ABHD median %.0f, other median %.0f" % (uni[uni.abhd_named].dom_len.median(), uni[~uni.abhd_named].dom_len.median()))
    rows, summ = [], []
    for sp in SPACES:
        M = pd.read_csv(C.DERIVED / f"dist_{sp}.csv", index_col=0)
        M = M.loc[~M.isna().all(1), ~M.isna().all(0)]
        accs = list(M.index); ab = uni.loc[accs, "abhd_named"].to_numpy(bool); genes = uni.loc[accs, "gene"].to_numpy()
        v = M.values.copy(); np.fill_diagonal(v, np.inf); v[np.isnan(v)] = np.inf
        R = int(ab.sum()) - 1
        res = []
        for i in np.where(ab)[0]:
            order = np.argsort(v[i]); order = order[np.isfinite(v[i, order])]
            if not len(order): continue
            nn = order[0]
            first_ab = next((j for j in order if ab[j]), None); first_non = next((j for j in order if not ab[j]), None)
            top = order[:R]
            res.append(dict(space=sp, gene=genes[i], nn=genes[nn], nn1_abhd=bool(ab[nn]),
                            nn_is_paralog=bool(ab[nn] and stem(genes[nn]) == stem(genes[i]) and genes[nn] != genes[i]),
                            r_precision=float(ab[top].mean()),
                            nonabhd_first=bool(first_non is not None and (first_ab is None or v[i, first_non] < v[i, first_ab]))))
        r = pd.DataFrame(res); rows.append(r)
        # same, but only neighbours whose domain length is within 20% of the query's (controls for size effects in TM-score)
        L = uni.loc[accs, "dom_len"].to_numpy(float); lm = []
        for i in np.where(ab)[0]:
            order = [j for j in np.argsort(v[i]) if np.isfinite(v[i, j]) and abs(L[j] / L[i] - 1) <= 0.2 and not (ab[j] and stem(genes[j]) == stem(genes[i]))]
            if order: lm.append(bool(ab[order[0]]))
        # same, ignoring paralogs (neighbours with the same stem are skipped)
        excl = []
        for i in np.where(ab)[0]:
            order = [j for j in np.argsort(v[i]) if np.isfinite(v[i, j]) and not (ab[j] and stem(genes[j]) == stem(genes[i]))]
            if order: excl.append(bool(ab[order[0]]))
        summ.append(dict(space=sp, n_abhd=len(r), nn1_abhd=round(r.nn1_abhd.mean(), 3),
                         nn1_abhd_excluding_paralogs=round(float(np.mean(excl)), 3),
                         nn1_abhd_excl_paralogs_length_matched=round(float(np.mean(lm)), 3), n_length_matched=len(lm),
                         nn_is_paralog=round(r.nn_is_paralog.mean(), 3),
                         r_precision=round(r.r_precision.mean(), 3),
                         expected=round(R / (len(accs) - 1), 3)))
    pd.concat(rows).to_csv(C.DERIVED / "abhd_label_by_gene.csv", index=False)
    s = pd.DataFrame(summ); s.to_csv(C.DERIVED / "abhd_label_summary.csv", index=False); print(s.to_string(index=False))

    # FunFam purity
    fl = uni.dropna(subset=["funfam_ids"]).copy()
    fl["ff"] = fl.funfam_ids.str.split(";")
    ex = fl.explode("ff")
    g = ex.groupby("ff").agg(n=("gene", "size"), n_abhd=("abhd_named", "sum"), genes=("gene", lambda x: ",".join(sorted(x))))
    g = g[g.n_abhd > 0]; g["pure_abhd"] = g.n == g.n_abhd
    g.to_csv(C.DERIVED / "funfam_purity.csv")
    print(f"FunFams with an ABHD gene: {len(g)}; ABHD-only: {int(g.pure_abhd.sum())}; mixed: {int((~g.pure_abhd).sum())}")
    print(g.to_string())

    fig = go.Figure()
    fig.add_bar(x=s.space, y=s.nn1_abhd, name="nearest neighbour is ABHD", marker_color="crimson")
    fig.add_bar(x=s.space, y=s.nn1_abhd_excluding_paralogs, name="same, skipping same-stem paralogs", marker_color="salmon")
    fig.add_bar(x=s.space, y=s.r_precision, name="R-precision (compactness of the ABHD set)", marker_color="steelblue")
    fig.add_bar(x=s.space, y=s.expected, name="expected by chance", marker_color="lightgrey")
    fig.update_layout(barmode="group", template="simple_white", width=900, height=450,
                      title="Stage 6b: paralog coherence vs superfamily-wide coherence of the ABHD name", yaxis_title="fraction")
    fig.write_html(C.FIGURES / "stage6b_abhd_label.html")


if __name__ == "__main__":
    main()
