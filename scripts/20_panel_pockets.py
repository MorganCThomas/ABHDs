#!/usr/bin/env python3
"""
Follow-up B (docs/followup_design.md): monoacylglycerol-hydrolase panel - oxyanion donors and aligned pocket-lining positions.

B1 oxyanion donors: Nu+1 backbone N (local) and the upstream non-local backbone N (second donor), residue identity, motif context,
   whether the second donor's atoms line the active-site cavity.
B3 pocket positions: each gene's cavity-lining residues mapped by TM-align onto every other panel gene; pairwise pocket identity
   vs global aligned identity (pocket-specific similarity = difference).
Outputs: data/derived/followup_b1_oxyanion_panel.csv, b1_donor2_all.csv, b3_pocket_pair_identity.csv, b3_mgll_positions.csv;
figures/followup_b_pockets.html
"""
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

from abhd import structops as a19

PANEL = ["MGLL", "ABHD6", "ABHD10", "ABHD12", "ABHD16A", "ABHD2", "ABHD1", "ABHD3", "ABHD4", "ABHD13", "ABHD14B", "ABHD17A", "LYPLA1"]
CLASS = {**{a: "hydrophobic" for a in "AVILMC"}, **{a: "aromatic" for a in "FWY"}, **{a: "polar" for a in "STNQH"},
         **{a: "charged" for a in "DEKR"}, "G": "small", "P": "small"}


rng = np.random.default_rng(3)


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv"); gacc = dict(zip(uni.gene, uni.accession))
    pd_ = pd.read_csv(C.DERIVED / "pocket_descriptors.csv").set_index("gene")
    ox = pd.read_csv(C.DERIVED / "oxyanion_donors.csv")
    lin = pd.read_csv(C.DERIVED / "pocket_lining.csv")
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv").set_index("gene")

    # ---- B1: second donor = nearest upstream (offset<0) non-local backbone N within 5.5 A of the nucleophile atom ----
    from Bio.PDB import PDBParser
    rows = []
    for g in uni.gene:
        if g not in cs.index or pd.isna(cs.loc[g, "elbow_resnum"]) or str(cs.loc[g, "category"]).startswith("no nucleophile"):
            continue
        acc = gacc[g]; n = int(cs.loc[g, "elbow_resnum"])
        s = PDBParser(QUIET=True).get_structure(acc, str(C.DERIVED / "domains" / f"{acc}.pdb")); res = {x.get_id()[1]: x for x in s.get_residues()}
        nuc = [a for a in res[n] if a.get_name() in ("OG", "SG", "OG1", "OD1")][0].coord
        best = None
        for i, r in res.items():
            if i < n - 4 and "N" in r:
                d = float(np.linalg.norm(r["N"].coord - nuc))
                if d <= 5.5 and (best is None or d < best[0]):
                    best = (d, i)
        d1 = float(np.linalg.norm(res[n + 1]["N"].coord - nuc)) if n + 1 in res else np.nan
        row = dict(gene=g, abhd_named=bool(cs.loc[g, "abhd_named"]), nuc=n, donor1_nu_plus1_dist=round(d1, 2))
        if best:
            i = best[1]
            from abhd.elbow import THREE_TO_ONE as T
            ctx = "".join(T.get(res[j].get_resname(), "-") if j in res else "-" for j in range(i - 3, i + 3))
            row.update(donor2_resnum=i, donor2_aa=T.get(res[i].get_resname(), "X"), donor2_dist=round(best[0], 2), donor2_context=ctx,
                       donor2_lines_cavity=bool(((lin.gene == g) & (lin.resnum == i)).any()))
        rows.append(row)
    B1 = pd.DataFrame(rows); B1.to_csv(C.DERIVED / "followup_b1_donor2_all.csv", index=False)
    pan = B1[B1.gene.isin(PANEL)].set_index("gene").reindex(PANEL).reset_index()
    pan = pan.merge(pd_[["vol_3.0", "vol_3.2", "vol_3.6", "touches_box_3.2", "n_lining", "frac_hydrophobic", "n_aromatic", "n_charged", "vol_resid_z"]].reset_index(), on="gene")
    pan.to_csv(C.DERIVED / "followup_b1_oxyanion_panel.csv", index=False)
    pd.set_option("display.width", 250)
    print(pan[["gene", "donor1_nu_plus1_dist", "donor2_aa", "donor2_resnum", "donor2_dist", "donor2_context", "donor2_lines_cavity",
               "vol_3.6", "vol_3.2", "touches_box_3.2", "n_lining", "frac_hydrophobic", "vol_resid_z"]].to_string(index=False))
    print("donor2 residue counts, all nucleophile genes:", B1.donor2_aa.value_counts().to_dict())
    print("donor2 by ABHD-named:", B1.groupby("abhd_named").donor2_aa.value_counts().unstack(0).fillna(0).astype(int).to_dict())
    print("share of genes whose second donor lines the cavity:", round(B1.donor2_lines_cavity.mean(), 2))

    # ---- B3: pairwise pocket identity vs global aligned identity ----
    L = {g: lin[lin.gene == g].set_index("resnum").aa.to_dict() for g in PANEL}
    seqs = {g: a19.load(gacc[g]) for g in PANEL}
    out, mg = [], []
    for g1, g2 in combinations(PANEL, 2):
        a, b = gacc[g1], gacc[g2]
        (ia, ca, sa), (ib, cb, sb), X, pairs, r = a19.align(a, b)
        i2j = {i: j for i, j in pairs}
        idn = np.mean([sa[i] == sb[j] for i, j in pairs])
        aa_of = {int(n): s for n, s in zip(ia, sa)}; bb_of = {int(n): s for n, s in zip(ib, sb)}
        num2idx = {int(n): k for k, n in enumerate(ia)}
        same, cls, m = 0, 0, 0
        for resnum in L[g1]:
            k = num2idx.get(resnum)
            if k is not None and k in i2j:
                m += 1; x, y = sa[k], sb[i2j[k]]
                same += x == y; cls += CLASS.get(x) == CLASS.get(y)
                if g1 == "MGLL":
                    mg.append(dict(mgll_pos=f"{x}{resnum}", partner=g2, partner_aa=y))
        # and the reverse direction
        i2j_r = {j: i for i, j in pairs}; num2idx_b = {int(n): k for k, n in enumerate(ib)}
        for resnum in L[g2]:
            k = num2idx_b.get(resnum)
            if k is not None and k in i2j_r:
                m += 1; x, y = sb[k], sa[i2j_r[k]]
                same += x == y; cls += CLASS.get(x) == CLASS.get(y)
        # null: same number of positions drawn at random from the structurally superposed core (CA displacement < 3 A)
        core = [(i, j) for i, j in pairs if np.linalg.norm(X[i] - cb[j]) < 3.0]
        ident = np.array([sa[i] == sb[j] for i, j in core]); k = min(max(m // 2, 1), len(ident))
        null = np.array([ident[rng.choice(len(ident), k, replace=False)].mean() for _ in range(500)])
        out.append(dict(core_identity=round(float(ident.mean()), 3), core_null_p95=round(float(np.percentile(null, 95)), 3), a=g1, b=g2, tm_rmsd=round(r.rmsd, 2), global_identity=round(float(idn), 3), n_pocket_pos=m,
                        pocket_identity=round(same / m, 3) if m else np.nan, pocket_class_identity=round(cls / m, 3) if m else np.nan))
    P = pd.DataFrame(out); P["pocket_minus_global"] = (P.pocket_identity - P.global_identity).round(3)
    P["pocket_minus_core"] = (P.pocket_identity - P.core_identity).round(3)
    P["pocket_above_core_null95"] = P.pocket_identity > P.core_null_p95
    P.to_csv(C.DERIVED / "followup_b3_pocket_pair_identity.csv", index=False)
    pd.DataFrame(mg).to_csv(C.DERIVED / "followup_b3_mgll_positions.csv", index=False)
    print(P.sort_values("pocket_identity", ascending=False).head(12).to_string(index=False))
    mgp = P[(P.a == "MGLL")].set_index("b")[["global_identity", "pocket_identity", "pocket_class_identity", "pocket_minus_global", "n_pocket_pos"]]
    print(mgp.to_string()); print("mean pocket_minus_global over all panel pairs:", round(P.pocket_minus_global.mean(), 3),
                              "| positive in", int((P.pocket_minus_global > 0).sum()), "of", len(P))
    print("pocket identity above the core-position null (95th pct) in", int(P.pocket_above_core_null95.sum()), "of", len(P),
          "pairs; mean pocket-minus-core:", round(P.pocket_minus_core.mean(), 3))
    from scipy.stats import fisher_exact
    d2 = B1.dropna(subset=["donor2_aa"])
    tab = pd.crosstab(d2.abhd_named, d2.donor2_aa == "N")
    print("Asn second donor, ABHD-named vs other:", tab.to_dict(), "Fisher p =", round(fisher_exact(tab.values)[1], 5))

    # ---- figure ----
    f = go.Figure()
    f.add_bar(x=pan.gene, y=pan["vol_3.6"], name="cavity volume, clearance 3.6 A", marker_color="steelblue")
    f.add_bar(x=pan.gene, y=pan["vol_3.2"], name="clearance 3.2 A", marker_color="lightsteelblue")
    f.update_layout(barmode="group", template="simple_white", width=900, height=420, yaxis_title="grid volume (A^3)",
                    title="Follow-up B: active-site cavity volume across the panel (apo AlphaFold models)")
    C.FIGURES.mkdir(exist_ok=True); f.write_html(C.FIGURES / "followup_b_pockets.html")


if __name__ == "__main__":
    main()
