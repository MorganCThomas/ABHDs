"""Pairwise comparison of two proteins from the study: similarity in every view, aligned catalytic residues and pocket-lining residues."""
import numpy as np
import pandas as pd
from scipy.stats import percentileofscore
from tmtools import tm_align

from .atlas import VIEWS, Atlas
from .structops import align, load

CLASS = {**{x: "hydrophobic" for x in "AVILMC"}, **{x: "aromatic" for x in "FWY"}, **{x: "polar" for x in "STNQH"}, **{x: "charged" for x in "DEKR"}, "G": "small", "P": "small"}


def compare(A: Atlas, ga, gb):
    a, b = A.acc(ga), A.acc(gb)
    ga, gb = A.gene(a), A.gene(b)
    out = dict(a=ga, b=gb)
    (ia, ca, sa), (ib, cb, sb), X, pairs, r = align(a, b)
    out["structure"] = dict(tm_a=float(r.tm_norm_chain1), tm_b=float(r.tm_norm_chain2), rmsd=float(r.rmsd), n_aligned=len(pairs), len_a=len(ia), len_b=len(ib))
    out["identity"] = A.identity(a, b)
    out["same_family"] = A.same_family(a, b); out["shared_funfams"] = sorted(A.ffs[a] & A.ffs[b])
    # distances and where the pair sits among all pairs / each other's neighbour lists
    rows = []
    for v in A.available_views():
        M = A.dist(v)
        d = M.loc[a, b]
        if np.isnan(d):
            continue
        allp = M.values[np.triu_indices(len(M), 1)]; allp = allp[np.isfinite(allp)]
        ra, rb = A.rank_in(ga, v).get(b, np.nan), A.rank_in(gb, v).get(a, np.nan)
        rows.append(dict(view=v, distance=float(d), percentile_of_all_pairs=float(percentileofscore(allp, d)), rank_of_b_for_a=ra, rank_of_a_for_b=rb))
    out["views"] = pd.DataFrame(rows)
    # aligned catalytic residues
    i2j = {i: j for i, j in pairs}; j2i = {j: i for i, j in pairs}
    num_a = {int(n): k for k, n in enumerate(ia)}; num_b = {int(n): k for k, n in enumerate(ib)}
    disp = lambda i, j: float(np.linalg.norm(X[i] - cb[j]))
    cat = []
    for role, key in (("nucleophile position", "elbow_resnum"), ("His", "his"), ("acid", "acid")):
        va = A.cs.loc[a, key] if a in A.cs.index else np.nan
        vb = A.cs.loc[b, key] if b in A.cs.index else np.nan
        row = dict(role=role, a_residue="—", b_aligned_residue="—", b_own_residue="—", ca_displacement="—")
        if not np.isnan(va) and int(va) in num_a:
            ka = num_a[int(va)]; row["a_residue"] = f"{sa[ka]}{int(va)}"
            if ka in i2j:
                row["b_aligned_residue"] = f"{sb[i2j[ka]]}{int(ib[i2j[ka]])}"; row["ca_displacement"] = f"{disp(ka, i2j[ka]):.1f} Å"
            else:
                row["b_aligned_residue"] = "unaligned"
        if not np.isnan(vb) and int(vb) in num_b:
            row["b_own_residue"] = f"{sb[num_b[int(vb)]]}{int(vb)}"
        cat.append(row)
    out["catalytic"] = pd.DataFrame(cat)
    out["category"] = (str(A.cs.loc[a, "category"]), str(A.cs.loc[b, "category"]))
    d2 = A.table("followup_b1_donor2_all.csv")
    donors = {}
    if d2 is not None:
        for g in (ga, gb):
            r_ = d2[d2.gene == g]
            donors[g] = f"{r_.donor2_aa.iat[0]}{int(r_.donor2_resnum.iat[0])}" if len(r_) and isinstance(r_.donor2_aa.iat[0], str) else "—"
    out["second_donor"] = donors
    # pocket-lining residues of A mapped onto B
    pl = A.table("pocket_lining.csv")
    if pl is not None:
        la = pl[pl.gene == ga]; lb = set(pl[pl.gene == gb].resnum)
        rows, same, cls = [], 0, 0
        for rr in la.itertuples():
            k = num_a.get(int(rr.resnum))
            if k is None or k not in i2j:
                continue
            j = i2j[k]; bres = sb[j]; same += rr.aa == bres; cls += CLASS.get(rr.aa) == CLASS.get(bres)
            rows.append(dict(a=f"{rr.aa}{rr.resnum}", b_aligned=f"{bres}{int(ib[j])}", identical=("yes" if rr.aa == bres else ""), b_also_lines_pocket=("yes" if int(ib[j]) in lb else "")))
        out["pocket"] = pd.DataFrame(rows)
        out["pocket_summary"] = dict(n=len(rows), identical=same, same_class=cls)
    cs = A.table("cavity_summary.csv")
    if cs is not None:
        out["cavity"] = {g: (float(cs[cs.gene == g]["median_3.2"].iat[0]) if (cs.gene == g).any() else None) for g in (ga, gb)}
    out["ec3"] = {ga: sorted(A.ec3[a]), gb: sorted(A.ec3[b])}
    return out


def _md(df):
    """Markdown table without the optional tabulate dependency."""
    cols = list(df.columns)
    rows = [[("" if (isinstance(v, float) and np.isnan(v)) else (f"{v:.3f}" if isinstance(v, float) else str(v))) for v in r] for r in df.itertuples(index=False)]
    return "\n".join(["| " + " | ".join(cols) + " |", "|" + "---|" * len(cols)] + ["| " + " | ".join(r) + " |" for r in rows])


def to_markdown(c):
    a, b = c["a"], c["b"]; s = c["structure"]
    L = [f"# {a} vs {b}", "",
         f"- Structure: TM-score {s['tm_a']:.2f} (normalised by {a}) / {s['tm_b']:.2f} (by {b}); RMSD {s['rmsd']:.1f} Å over {s['n_aligned']} aligned residues ({s['len_a']} and {s['len_b']} in the domains)",
         f"- Sequence identity (global alignment of the trimmed domains): {c['identity']:.0%}",
         f"- Same paralog group or FunFam: {'yes' if c['same_family'] else 'no'}" + (f" (shared FunFams {', '.join(c['shared_funfams'])})" if c["shared_funfams"] else ""),
         f"- Catalytic categories: {a} {c['category'][0]}; {b} {c['category'][1]}",
         f"- EC3: {a} {', '.join(c['ec3'][a]) or 'none'}; {b} {', '.join(c['ec3'][b]) or 'none'}",
         f"- Second oxyanion donor: {a} {c['second_donor'].get(a, '—')}; {b} {c['second_donor'].get(b, '—')}"]
    if c.get("cavity"):
        L.append(f"- Active-site cavity (median grid units, 3.2 Å, like-for-like): {a} {c['cavity'][a] or '—'}; {b} {c['cavity'][b] or '—'} (model-dependent; see report 12)")
    L += ["", "## Distance in each view", "", _md(c["views"]), "",
          "## Catalytic residues (aligned by TM-align)", "", _md(c["catalytic"])]
    if "pocket" in c and len(c["pocket"]):
        ps = c["pocket_summary"]
        L += ["", f"## Pocket-lining residues of {a} mapped onto {b}", "",
              f"{ps['identical']} of {ps['n']} identical, {ps['same_class']} of the same physicochemical class (random core positions of this protein pair would be identical in roughly 15–25%).", "",
              _md(c["pocket"])]
    return "\n".join(L)
