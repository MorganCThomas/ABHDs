#!/usr/bin/env python3
"""
Motif audit (reports/14_motif_audit.md): do the historical motifs HX4D ("acyltransferase") and GXSXG (nucleophile elbow) mean what they were taken to mean?

  1. HX4D at the catalytic His: how many triads start with H-x-x-x-x-D, and is that more often than at other histidines of the same proteins
     (null: one random non-catalytic His per protein, 20,000 draws), and than D at the neighbouring offsets (His+4, His+6)?
  2. Is the Asp of HX4D ever the triad acid? (In chain order the acid precedes the His, so it should not be.)
  3. Does HX4D go with acyltransferase annotation (EC 2.3.*) in this set?
  4. The 20 ABHD-named genes: elbow context, GXSXG / small-X-Nu-X-small, HX4D, and orthologue conservation of the His and the HX4D Asp where orthologue sets exist.
Outputs: data/derived/motif_audit_triads.csv, motif_audit_offsets.csv, motif_audit_summary.csv, motif_audit_orthologues.csv
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.atlas import Atlas
from abhd.orthologs import ORTH, conservation
from abhd.structops import load_pdb

SMALL = set("GASCT")
rng = np.random.default_rng(0)


def main():
    A = Atlas(); cs = A.cs
    rows, offs_t, offs_bg, bg_all, bg_acidic, bg_small = [], [], [], [], [], []
    seqs = {}
    for acc in cs.index:
        ids, _, s = load_pdb(C.STRUCTURES / f"{acc}.pdb"); seqs[acc] = (np.array(ids), s)
    for acc, c in cs.iterrows():
        if pd.isna(c["his"]):
            continue
        ids, s = seqs[acc]; g = A.gene(acc)
        k = int(np.where(ids == int(c["his"]))[0][0])
        has = lambda i, aa: 0 <= i < len(s) and s[i] == aa
        hx4d = has(k + 5, "D"); hx4e = has(k + 5, "E")
        acid_k = int(np.where(ids == int(c["acid"]))[0][0]) if not pd.isna(c["acid"]) else None
        n = None
        if not pd.isna(c["elbow_resnum"]):
            n = int(np.where(ids == int(c["elbow_resnum"]))[0][0])
        ctx = s[n - 2:n + 3] if n is not None and n >= 2 else ""
        gxsxg = bool(re.fullmatch(r"G.[SC].G", ctx))
        smallxnuxsmall = len(ctx) == 5 and ctx[0] in SMALL and ctx[4] in SMALL
        bg_small += [s[i - 2] in SMALL and s[i + 2] in SMALL for i in range(2, len(s) - 2) if s[i] == "S"]
        others = [i for i, a in enumerate(s) if a == "H" and i != k]
        rows.append(dict(gene=g, accession=acc, abhd_named=bool(A.uni.set_index("accession").loc[acc, "abhd_named"]), triad=c["category"], his=f"H{int(c['his'])}",
                         window=s[k:k + 6], hx4d=hx4d, hx4e=hx4e, funfam=c["funfam"], d5_is_acid=(acid_k == k + 5), n_other_his=len(others),
                         bg_hx4d_frac=np.mean([has(i + 5, "D") for i in others]) if others else np.nan,
                         elbow=ctx, gxsxg=gxsxg, small_x_nu_x_small=smallxnuxsmall, ec=A.uni.set_index("accession").loc[acc, "ec"]))
        for off in range(-10, 11):
            offs_t.append(dict(offset=off, is_D=has(k + off, "D")))
        for i in others:
            bg_all.append(has(i + 5, "D")); bg_acidic.append(0 <= i + 5 < len(s) and s[i + 5] in "DE")
            for off in range(-10, 11):
                offs_bg.append(dict(offset=off, is_D=has(i + off, "D")))
    T = pd.DataFrame(rows)
    T.to_csv(C.DERIVED / "motif_audit_triads.csv", index=False)
    n, h = len(T), int(T.hx4d.sum())
    # null: one random non-catalytic His per protein
    # vectorised: per-protein probability that a random non-catalytic His has D at +5; the sum of independent Bernoullis is Poisson-binomial
    p = T.bg_hx4d_frac.fillna(0).values
    sims = (rng.random((20000, len(p))) < p).sum(1)
    p_perm = float((sims >= h).mean() + 1 / 20001)
    offs = pd.concat([pd.DataFrame(offs_t).groupby("offset").is_D.mean().rename("triad_his"), pd.DataFrame(offs_bg).groupby("offset").is_D.mean().rename("other_his")], axis=1)
    bg_acid = float(np.mean([has_ for has_ in bg_acidic]))
    offs.to_csv(C.DERIVED / "motif_audit_offsets.csv")
    ec23 = T.ec.fillna("").str.contains(r"\b2\.3\.")
    a, b = int((T.hx4d & ec23).sum()), int((~T.hx4d & ec23).sum())
    c_, d_ = int((T.hx4d & ~ec23).sum()), int((~T.hx4d & ~ec23).sum())
    fe = fisher_exact([[a, b], [c_, d_]])
    ham = T.abhd_named
    summ = {
        "triads": n, "hx4d_at_triad_his": h, "hx4d_fraction": round(h / n, 3),
        "hx4d_at_other_his_fraction": round(float(np.mean(bg_all)), 3), "null_expected_triads": round(float(p.sum()), 1), "perm_p_one_sided": round(p_perm, 4),
        "d_at_his+4": round(float(offs.loc[4, "triad_his"]), 3), "d_at_his+5": round(float(offs.loc[5, "triad_his"]), 3), "d_at_his+6": round(float(offs.loc[6, "triad_his"]), 3),
        "hx4e_at_triad_his": int(T.hx4e.sum()), "hx4_acidic_at_triad_his": int((T.hx4d | T.hx4e).sum()),
        "funfams_with_triad": int(T.funfam.nunique()), "funfams_with_hx4d": int(T[T.hx4d].funfam.nunique()),
        "acidic_at_other_his_fraction": round(bg_acid, 3), "small_x_nu_x_small_at_random_ser": round(float(np.mean(bg_small)), 3),
        "hx4d_D_is_triad_acid": int(T.d5_is_acid.sum()),
        "ec23_proteins": int(ec23.sum()), "ec23_with_hx4d": a, "non_ec23_with_hx4d": c_, "non_ec23": c_ + d_, "fisher_ec23_p": round(float(fe[1]), 3),
        "hx4d_in_abhd_named": f"{int((T.hx4d & ham).sum())}/{int(ham.sum())}", "hx4d_in_non_abhd": f"{int((T.hx4d & ~ham).sum())}/{int((~ham).sum())}",
        "gxsxg_all": f"{int(T.gxsxg.sum())}/{n}", "small_x_nu_x_small_all": f"{int(T.small_x_nu_x_small.sum())}/{n}",
        "gxsxg_abhd_named": f"{int((T.gxsxg & ham).sum())}/{int(ham.sum())}",
        "small_x_nu_x_small_abhd_named": f"{int((T.small_x_nu_x_small & ham).sum())}/{int(ham.sum())}",
    }
    pd.Series(summ).rename("value").to_csv(C.DERIVED / "motif_audit_summary.csv")
    print(pd.Series(summ).to_string())
    print(offs.round(3).T.to_string())
    print(T[T.abhd_named][["gene", "triad", "his", "window", "hx4d", "elbow", "gxsxg", "small_x_nu_x_small"]].to_string(index=False))
    # ABHD5 / ABHD15 have no nucleophile, so no triad call; their His is the residue TM-align places on the partner's catalytic His (followup A2)
    A2 = pd.read_csv(C.DERIVED / "followup_a2_catalytic_positions.csv"); extra = {}
    for r in A2[A2.role == "His"].itertuples():
        g = r.pair.split(" vs ")[0]
        if g in ("ABHD5", "ABHD15") and g not in extra:
            extra[g] = int(r.partner_residue[1:])
    for g, hnum in extra.items():
        ids, s = seqs[A.acc(g)]; k = int(np.where(ids == hnum)[0][0])
        print(f"{g}: His{hnum} window {s[k:k+6]} (HX4D: {s[k+5] == 'D'}, HX4[DE]: {s[k+5] in 'DE'})")
    # orthologue conservation for genes with cached orthologue sets
    out = []
    for f in sorted(ORTH.glob("*_v2.tsv")):
        g = f.stem[:-3]
        try:
            r = conservation(A, g)
        except Exception as e:
            print("skip", g, e); continue
        if "error" in r:
            continue
        pp = r["per_position"].set_index("resnum")
        if g in set(T.gene):
            t = T[T.gene == g].iloc[0]; hnum = int(t.his[1:]); hx = bool(t.hx4d)
        elif g in extra:
            hnum = extra[g]; hx = pp.loc[hnum + 5, "aa"] == "D"
        else:
            continue
        out.append(dict(gene=g, n_orthologues=r["n_orthologues"], his=f"H{hnum}", his_conservation=pp.loc[hnum, "conservation"], d_plus5=f"{pp.loc[hnum + 5, 'aa']}{hnum + 5}",
                        d_plus5_conservation=pp.loc[hnum + 5, "conservation"], d_plus5_percentile=pp.loc[hnum + 5, "percentile"], human_hx4d=hx))
    if out:
        O = pd.DataFrame(out); O.to_csv(C.DERIVED / "motif_audit_orthologues.csv", index=False); print(O.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
