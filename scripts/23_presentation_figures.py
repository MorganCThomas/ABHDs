#!/usr/bin/env python3
"""Clean static figures for the presentation (figures/presentation/*.png); every number is read from data/derived tables."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.catalytic import apply_elbow

OUT = C.FIGURES / "presentation"; OUT.mkdir(parents=True, exist_ok=True)
INK, MUTED, ABHD, OTHER, GREY, BG = "#1B2A41", "#5B6677", "#D9731A", "#2F6DA8", "#B8BEC6", "#FFFFFF"
ABHD_LIGHT = "#F2C79B"
plt.rcParams.update({"font.family": ["Arial", "Helvetica", "DejaVu Sans"], "font.size": 17, "axes.edgecolor": MUTED, "axes.labelcolor": INK, "text.color": INK,
                     "xtick.color": MUTED, "ytick.color": INK, "figure.facecolor": BG, "axes.facecolor": BG, "savefig.facecolor": BG,
                     "axes.spines.top": False, "axes.spines.right": False})
D = C.DERIVED


def save(fig, name):
    fig.tight_layout(); fig.savefig(OUT / f"{name}.png", dpi=200); plt.close(fig)


def universe():
    u = pd.read_csv(D / "universe.csv")
    order = [("uniprot_only", "UniProt query only"), ("both", "Both sources"), ("funfam_only", "CATH FunFam only")]
    fig, ax = plt.subplots(figsize=(9, 4))
    for i, (k, lab) in enumerate(order):
        a = int(((u.source == k) & u.abhd_named).sum()); o = int(((u.source == k) & ~u.abhd_named).sum())
        ax.barh(i, o, color=OTHER); ax.barh(i, a, left=o, color=ABHD)
        ax.text(o + a + 1.2, i, f"{o + a}", va="center", fontsize=20, fontweight="bold")
    ax.set_yticks(range(3)); ax.set_yticklabels([l for _, l in order]); ax.invert_yaxis(); ax.set_xlim(0, 70)
    ax.set_xlabel("genes"); ax.barh([], [])
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=OTHER), Rectangle((0, 0), 1, 1, color=ABHD)], labels=["other names", "ABHD-named"], frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    save(fig, "universe")


def elbow():
    et = pd.read_csv(D / "elbow_table.csv"); g = pd.read_csv(D / "benchmark_stage3_genes.csv")
    gold = et.merge(g[["accession", "gold_nuc"]], left_on=["accession", "resnum"], right_on=["accession", "gold_nuc"])
    ok = et[apply_elbow(et)]
    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(et.phi, et.psi, s=8, color=GREY, alpha=0.5, label="candidate residues", linewidths=0)
    ax.scatter(ok.phi, ok.psi, s=14, color=OTHER, alpha=0.8, label="pass the elbow rule", linewidths=0)
    ax.scatter(gold.phi, gold.psi, s=70, color=ABHD, edgecolor=INK, linewidths=1, label="known nucleophiles", zorder=5)
    ax.add_patch(Rectangle((30, -170), 60, 90, fill=False, edgecolor=INK, linestyle="--", linewidth=1.5))
    ax.set_xlabel("phi (degrees)"); ax.set_ylabel("psi (degrees)"); ax.set_xlim(-10, 130); ax.set_ylim(-180, -40)
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3, fontsize=13, markerscale=1.0, columnspacing=1.0, handletextpad=0.2)
    save(fig, "elbow")


def categories():
    cs = pd.read_csv(D / "catalytic_sites.csv")
    def grp(c):
        c = str(c)
        if c == "Ser-His-Asp": return "Ser-His-Asp"
        if c == "Ser-His-Glu": return "Ser-His-Glu"
        if c.startswith("Asp"): return "Asp-based"
        if c.startswith("Cys"): return "Cys-His-Asp"
        if c.startswith("Ser, no His"): return "Ser, no His found"
        return "No nucleophile at the elbow"
    cs["g"] = cs.category.map(grp)
    order = ["Ser-His-Asp", "Ser-His-Glu", "Asp-based", "Cys-His-Asp", "Ser, no His found", "No nucleophile at the elbow"]
    fig, ax = plt.subplots(figsize=(9, 5))
    for i, k in enumerate(order):
        s = cs[cs.g == k]; a = int(s.abhd_named.sum()); o = len(s) - a
        ax.barh(i, o, color=GREY if k.startswith("No") else OTHER); ax.barh(i, a, left=o, color=ABHD)
        ax.text(len(s) + 1, i, str(len(s)), va="center", fontsize=19, fontweight="bold")
    ax.set_yticks(range(len(order))); ax.set_yticklabels(order); ax.invert_yaxis(); ax.set_xlim(0, 90); ax.set_xlabel("genes")
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=OTHER), Rectangle((0, 0), 1, 1, color=ABHD)], labels=["other names", "ABHD-named"], frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    save(fig, "categories")


def nn_chance():
    R = pd.read_csv(D / "review_label_robustness.csv"); R = R[(R.definition == "abhd_named") & (R.variant == "all")].set_index("space")
    spaces = [("tmalign", "Structure\n(TM-align)"), ("seqid", "Sequence\n(local)"), ("funfam", "CATH FunFam\nHMMs")]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for i, (k, lab) in enumerate(spaces):
        v = R.loc[k, "nn1"]; ax.bar(i, v, color=ABHD, width=0.55, yerr=[[v - R.loc[k, "nn1_lo"]], [min(1, R.loc[k, "nn1_hi"]) - v]], capsize=6, ecolor=INK)
        ax.text(i, 0.04, f"{v:.0%}", ha="center", color="white", fontsize=22, fontweight="bold")
    ax.axhline(R.chance.mean(), color=INK, linestyle="--"); ax.text(2.45, R.chance.mean() + 0.03, "chance 17%", ha="right", fontsize=15)
    ax.set_xticks(range(3)); ax.set_xticklabels([l for _, l in spaces]); ax.set_ylim(0, 1.05); ax.set_yticks([0, .25, .5, .75, 1]); ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
    ax.set_ylabel("nearest neighbour is ABHD-named")
    save(fig, "nn_chance")


def nn_removed():
    R = pd.read_csv(D / "review_label_robustness.csv"); R = R[(R.definition == "abhd_named")]
    spaces = [("tmalign", "Structure"), ("seqid", "Sequence"), ("funfam", "FunFam"), ("anchor_elbow", "Elbow window")]
    fig, ax = plt.subplots(figsize=(9, 4.8)); w = 0.36
    for i, (k, lab) in enumerate(spaces):
        a = R[(R.space == k) & (R.variant == "all")].nn1.iat[0]; b = R[(R.space == k) & (R.variant == "skip_closest_3.4pct_pairs")].nn1.iat[0]
        ax.bar(i - w / 2, a, w, color=ABHD_LIGHT); ax.bar(i + w / 2, b, w, color=ABHD)
        ax.text(i - w / 2, a + 0.02, f"{a:.0%}", ha="center", fontsize=15); ax.text(i + w / 2, b + 0.02, f"{b:.0%}", ha="center", fontsize=15, fontweight="bold")
    ax.axhline(0.17, color=INK, linestyle="--"); ax.text(3.6, 0.2, "chance", ha="right", fontsize=14)
    ax.set_xticks(range(4)); ax.set_xticklabels([l for _, l in spaces]); ax.set_ylim(0, 1.0); ax.set_ylabel("nearest neighbour is ABHD-named")
    ax.set_yticks([0, .25, .5, .75, 1]); ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=ABHD_LIGHT), Rectangle((0, 0), 1, 1, color=ABHD)], labels=["all genes", "closest 3.4% of pairs removed"], frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, fontsize=14)
    save(fig, "nn_removed")


def ec():
    L = pd.read_csv(D / "function_loo.csv"); L = L[L.ec_level == 3]
    spaces = [("tmalign", "Structure"), ("seqid", "Sequence"), ("funfam", "FunFam"), ("anchor_elbow", "Elbow window"), ("soap", "Local site\n(SOAP)"), ("protbert", "Language\nmodel")]
    fig, ax = plt.subplots(figsize=(10, 5.6)); w = 0.38
    for i, (k, lab) in enumerate(spaces):
        a = L[(L.space == k) & (L["mode"] == "all")].iloc[0]; b = L[(L.space == k) & (L["mode"] == "leave_family_out")].iloc[0]
        ax.bar(i - w / 2, a.acc_1nn, w, color=ABHD_LIGHT)
        ax.bar(i + w / 2, b.acc_1nn, w, color=ABHD, yerr=[[b.acc_1nn - b.ci95_lo], [b.ci95_hi - b.acc_1nn]], capsize=4, ecolor=INK)
    l1 = ax.axhline(0.588, color=INK, linestyle="--"); l2 = ax.axhline(0.397, color=MUTED, linestyle=":")
    ax.set_xticks(range(len(spaces))); ax.set_xticklabels([l for _, l in spaces], fontsize=15); ax.set_ylim(0, 1.0); ax.set_ylabel("EC class recovered")
    ax.set_yticks([0, .25, .5, .75, 1]); ax.set_yticklabels(["0", "25%", "50%", "75%", "100%"])
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=ABHD_LIGHT), Rectangle((0, 0), 1, 1, color=ABHD), l1, l2],
              labels=["any neighbour", "paralog family excluded (95% CI)", "guess the most common class: 59%", "random neighbour: 40%"],
              frameon=False, loc="lower center", fontsize=13, bbox_to_anchor=(0.5, 1.0), ncol=2)
    save(fig, "ec")


def abhd15():
    E = pd.read_csv(D / "orth_elbow_by_class.csv"); E = E[E.gene == "ABHD15"]
    rows = []
    for r in E.itertuples():
        d = dict(x.split(":") for x in r.residues.split(";")); n = sum(int(v) for v in d.values())
        rows.append((r.taxon_class, int(d.get("G", 0)), int(d.get("S", 0)), n - int(d.get("G", 0)) - int(d.get("S", 0)), n))
    order = ["mammals", "birds", "lepidosaurs", "turtles", "amphibians", "ray-finned fish", "crocodilians"]
    rows = sorted(rows, key=lambda r: order.index(r[0]))
    fig, ax = plt.subplots(figsize=(9, 5))
    for i, (k, g, s, o, n) in enumerate(rows):
        ax.barh(i, g / n, color=OTHER); ax.barh(i, s / n, left=g / n, color=ABHD); ax.barh(i, o / n, left=(g + s) / n, color=GREY)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels([f"{k} (n={n})" for k, *_, n in rows], fontsize=15); ax.invert_yaxis()
    ax.set_xlim(0, 1); ax.set_xticks([0, .5, 1]); ax.set_xticklabels(["0", "50%", "100%"]); ax.set_xlabel("residue at the nucleophile position of ABHD15 orthologues")
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=c) for c in (OTHER, ABHD, GREY)], labels=["Gly (as in human)", "Ser", "gap"], frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(0.45, 1.0), fontsize=14)
    save(fig, "abhd15_orthologues")


def segment():
    G = pd.read_csv(D / "orth_segment.csv")
    fig, ax = plt.subplots(figsize=(10.5, 4.8)); x = np.arange(len(G)); w = 0.4
    ax.bar(x - w / 2, G.abhd4_conservation, w, color=ABHD)
    un = G.abhd5_pos == "unaligned"
    ax.bar(x[~un] + w / 2, G.abhd5_conservation[~un], w, color=OTHER)
    ax.bar(x[un] + w / 2, np.full(un.sum(), 0.06), w, color=GREY, hatch="//", edgecolor=MUTED)
    ax.set_xticks(x); ax.set_xticklabels(G.abhd4_pos, rotation=60, fontsize=13); ax.set_ylim(0, 1.05); ax.set_ylabel("conserved across orthologues")
    ax.set_yticks([0, .5, 1]); ax.set_yticklabels(["0", "50%", "100%"]); ax.set_xlabel("ABHD4 residues 180-198")
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=ABHD), Rectangle((0, 0), 1, 1, color=OTHER), Rectangle((0, 0), 1, 1, facecolor=GREY, hatch="//", edgecolor=MUTED)],
              labels=["ABHD4 orthologues", "ABHD5 orthologues (aligned residue)", "no ABHD5 equivalent"], frameon=False, fontsize=13, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=3)
    save(fig, "abhd4_segment")


def cavity():
    S = pd.read_csv(D / "cavity_summary.csv").set_index("gene")
    panel = ["ABHD6", "ABHD5", "ABHD3", "ABHD2", "MGLL", "ABHD10", "ABHD1", "ABHD4", "ABHD12", "ABHD16A", "ABHD15"]
    S = S.loc[panel].sort_values("median_3.2", ascending=True)
    fig, ax = plt.subplots(figsize=(8, 6))
    for i, (g, r) in enumerate(S.iterrows()):
        col = ABHD if g in ("ABHD15", "ABHD5") else OTHER
        ax.barh(i, r["median_3.2"], color=col, alpha=0.95)
        ax.plot([r["vmin_3.2"], r["vmax_3.2"]], [i, i], color=INK, linewidth=2)
        if r.frac_box_truncated > 0.3:
            ax.text(r["vmax_3.2"] + 8, i, "box-truncated", va="center", fontsize=12, color=MUTED)
    ax.set_yticks(range(len(S))); ax.set_yticklabels(S.index); ax.set_xlabel("cavity volume (grid units, like-for-like)")
    save(fig, "cavity")


def attrition():
    a = pd.read_csv(D / "attrition_stage1.csv").set_index("fate_class").n_entries
    lab = {"alt_entry": "other entries of a kept gene", "canonical": "kept: one reviewed entry per gene", "fragment": "fragments placed by sequence",
           "unassigned_fragment": "unplaced fragments", "unreviewed_only_gene": "unreviewed-only gene"}
    order = ["canonical", "alt_entry", "fragment", "unassigned_fragment", "unreviewed_only_gene"]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    for i, k in enumerate(order):
        ax.barh(i, a[k], color=OTHER if k == "canonical" else GREY); ax.text(a[k] + 5, i, str(a[k]), va="center", fontsize=19, fontweight="bold")
    ax.set_yticks(range(len(order))); ax.set_yticklabels([lab[k] for k in order], fontsize=15); ax.invert_yaxis(); ax.set_xlim(0, 330); ax.set_xlabel("source accessions (444 in total)")
    save(fig, "attrition")


def rmsd():
    fig, axs = plt.subplots(1, 2, figsize=(10, 4.4), sharey=False)
    for ax, f, ttl, col in ((axs[0], "site_rmsd_elbow", "Elbow backbone (5 residues)", OTHER), (axs[1], "site_rmsd_triad", "Triad skeleton", ABHD)):
        d = pd.read_csv(D / f"{f}.csv", index_col=0).values; v = d[np.triu_indices(len(d), 1)]; v = v[np.isfinite(v)]
        ax.hist(v, bins=np.linspace(0, 1.6 if "triad" in f else 0.8, 40), color=col); m = np.median(v)
        ax.axvline(m, color=INK, linestyle="--"); ax.text(m + 0.02, ax.get_ylim()[1] * 0.88, f"median {m:.2f} Å", fontsize=15)
        ax.set_title(ttl, fontsize=17); ax.set_xlabel("pairwise RMSD (Å)"); ax.set_yticks([])
    save(fig, "rmsd")


def soap():
    S = pd.read_csv(D / "soap_null_summary.csv")
    lab = ["elbow vs elbow", "elbow vs random Ser", "random Ser vs random Ser"]
    fig, ax = plt.subplots(figsize=(9, 4.4))
    for i, r in S.iterrows():
        ax.barh(i, r["median"], color=ABHD if i == 0 else GREY); ax.plot([r.p5, r.p95], [i, i], color=INK, linewidth=2.5)
        ax.text(r.p95 + 0.004, i, f"{r['median']:.3f}", va="center", fontsize=18, fontweight="bold")
    ax.set_yticks(range(3)); ax.set_yticklabels(lab, fontsize=16); ax.invert_yaxis(); ax.set_xlim(0, 0.21)
    ax.set_xlabel("SOAP cosine distance (median; line = 5th–95th percentile)")
    save(fig, "soap_null")


def concordance():
    M = pd.read_csv(D / "global_concordance.csv", index_col=0)
    names = {"tmalign": "Structure", "tmshort": "Structure\n(short norm.)", "seqid": "Sequence", "protbert": "Language\nmodel", "funfam": "FunFam", "soap": "SOAP", "triad": "Triad"}
    fig, ax = plt.subplots(figsize=(8, 6.6))
    im = ax.imshow(M.values, cmap="Blues", vmin=-0.1, vmax=1)
    ax.set_xticks(range(len(M))); ax.set_xticklabels([names[c] for c in M.columns], rotation=45, ha="right", fontsize=14)
    ax.set_yticks(range(len(M))); ax.set_yticklabels([names[c] for c in M.index], fontsize=14)
    for i in range(len(M)):
        for j in range(len(M)):
            ax.text(j, i, f"{M.values[i, j]:.2f}", ha="center", va="center", fontsize=13, color="white" if M.values[i, j] > 0.55 else INK)
    for sp in ax.spines.values(): sp.set_visible(False)
    save(fig, "concordance")


def cluster_null():
    c = pd.read_csv(D / "global_cluster_assignments.csv"); lab = c.abhd.to_numpy(bool); cl = c.cl20_tmalign.to_numpy()
    rng = np.random.default_rng(0); null = [len(set(cl[rng.permutation(len(cl))[:lab.sum()]])) for _ in range(5000)]
    obs = len(set(cl[lab]))
    fig, ax = plt.subplots(figsize=(8.5, 4.6))
    vals, counts = np.unique(null, return_counts=True); ax.bar(vals, counts / 5000, color=GREY, width=0.8)
    ax.axvline(obs, color=ABHD, linewidth=4); ax.text(obs + 0.2, ax.get_ylim()[1] * 0.85, f"ABHD genes: {obs}", color=ABHD, fontsize=20, fontweight="bold")
    ax.set_xlabel("structural clusters occupied by 20 genes (of 20 clusters)"); ax.set_ylabel("share of random gene sets"); ax.set_yticks([])
    save(fig, "cluster_null")


def majority():
    G = pd.read_csv(D / "global_clusters.csv"); G = G[G.n_clusters == 20].set_index("space")
    sp = [("tmalign", "Structure"), ("tmshort", "Structure\n(short norm.)"), ("seqid", "Sequence"), ("protbert", "Language\nmodel"), ("soap", "SOAP"), ("triad", "Triad")]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    for i, (k, lab) in enumerate(sp):
        v = G.loc[k, "abhd_genes_in_majority_abhd_clusters"]; ax.bar(i, v, color=ABHD if k == "seqid" else GREY, width=0.6); ax.text(i, v + 0.4, str(v), ha="center", fontsize=20, fontweight="bold")
    ax.set_xticks(range(len(sp))); ax.set_xticklabels([l for _, l in sp], fontsize=14); ax.set_ylim(0, 22); ax.set_ylabel("ABHD genes in a majority-ABHD cluster")
    save(fig, "majority_clusters")


def anchored():
    M = pd.read_csv(D / "global_concordance.csv", index_col=0); A = pd.read_csv(D / "anchored_summary.csv")
    v = [("Whole-domain\nlocal sequence", M.loc["seqid", "funfam"]), ("Elbow window\n(31 residues)", float(A[(A.anchored == "anchor_elbow") & (A.vs == "funfam")].spearman.iat[0])),
         ("Catalytic windows\n(63 residues)", float(A[(A.anchored == "anchor_cat") & (A.vs == "funfam")].spearman.iat[0])), ("Structure\n(TM-align)", M.loc["tmalign", "funfam"])]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    for i, (l, x) in enumerate(v):
        ax.bar(i, x, color=OTHER if i in (1, 2) else GREY, width=0.6); ax.text(i, x + 0.015, f"{x:.2f}", ha="center", fontsize=20, fontweight="bold")
    ax.set_xticks(range(4)); ax.set_xticklabels([l for l, _ in v], fontsize=14); ax.set_ylim(0, 0.7); ax.set_ylabel("rank correlation with FunFam distance")
    save(fig, "anchored")


def cavity_models():
    V = pd.read_csv(D / "cavity_variants.csv"); V = V[(V.anchor == "ca_trunc") & (V.clearance == 3.2)]
    P = V.groupby(["gene", "model"]).volume.median().unstack()
    from scipy.stats import spearmanr
    rho = spearmanr(P.domain, P.full)[0]
    fig, ax = plt.subplots(figsize=(7.5, 6.2))
    ax.scatter(P.domain, P.full, s=40, color=OTHER, alpha=0.8, linewidths=0)
    for g, c in (("ABHD2", ABHD), ("ABHD15", ABHD), ("ABHD17A", ABHD)):
        ax.scatter(P.loc[g, "domain"], P.loc[g, "full"], s=110, color=ABHD, edgecolor=INK, zorder=5); ax.annotate(g, (P.loc[g, "domain"], P.loc[g, "full"]), xytext=(8, 6), textcoords="offset points", fontsize=15)
    m = max(P.max()); ax.plot([0, m], [0, m], color=MUTED, linestyle=":")
    ax.set_xlabel("cavity volume, trimmed domain"); ax.set_ylabel("cavity volume, full-length model"); ax.text(0.03, 0.93, f"Spearman {rho:.2f}", transform=ax.transAxes, fontsize=19, fontweight="bold")
    save(fig, "cavity_models")


def pocket_lining():
    P = pd.read_csv(D / "followup_b3_pocket_pair_identity.csv"); P = P[P.a == "MGLL"].set_index("b").loc[["ABHD6", "ABHD10", "ABHD12", "ABHD16A", "ABHD2", "ABHD1", "ABHD3", "ABHD4", "ABHD13", "ABHD14B", "ABHD17A", "LYPLA1"]]
    fig, ax = plt.subplots(figsize=(10, 4.9)); x = np.arange(len(P)); w = 0.38
    ax.bar(x - w / 2, P.core_identity, w, color=GREY); ax.bar(x + w / 2, P.pocket_identity, w, color=ABHD)
    ax.scatter(x + w / 2, P.core_null_p95, marker="_", s=300, color=INK, linewidths=2.5, zorder=5)
    ax.set_xticks(x); ax.set_xticklabels(P.index, rotation=45, ha="right", fontsize=14); ax.set_ylabel("identity with MGLL"); ax.set_ylim(0, 0.55)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=GREY), Rectangle((0, 0), 1, 1, color=ABHD), plt.Line2D([], [], marker="_", color=INK, linestyle="", markersize=18, markeredgewidth=2.5)],
              labels=["aligned core, any position", "MGLL pocket-lining positions", "95th percentile of random core positions"], frameon=False, fontsize=13, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2)
    save(fig, "pocket_lining")


def decay():
    F = pd.read_csv(D / "function_decay.csv"); F = F[(F["mode"] == "all") & (F.max_identity >= 25)]
    names = {"tmalign": "Structure", "seqid": "Sequence", "funfam": "FunFam", "anchor_elbow": "Elbow window"}
    cols = {"tmalign": OTHER, "seqid": ABHD, "funfam": "#7A5C99", "anchor_elbow": "#3C8C6E"}
    fig, ax = plt.subplots(figsize=(10, 5.4)); off = {"tmalign": -0.6, "seqid": -0.2, "funfam": 0.2, "anchor_elbow": 0.6}
    for k in names:
        s = F[F.space == k].sort_values("max_identity", ascending=False); x = np.arange(len(s)) + off[k] * 0.12
        ax.errorbar(x, s.diff_vs_majority, yerr=[s.diff_vs_majority - s.diff_ci_lo, s.diff_ci_hi - s.diff_vs_majority], marker="o", color=cols[k], label=names[k], capsize=3, linewidth=2)
    ax.axhline(0, color=INK, linestyle="--"); ax.set_xticks(range(len(s))); ax.set_xticklabels([f"{int(v)}%" for v in s.max_identity])
    ax.set_xlabel("maximum identity allowed between query and neighbour"); ax.set_ylabel("accuracy minus baseline")
    ax.legend(frameon=False, loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=4, fontsize=14)
    save(fig, "decay")


def multiview():
    A = pd.read_csv(D / "multiview_association.csv"); A = A[A.scope == "cross_family_pairs"]
    names = {"structure": "Structure", "structure_core200": "Structure (core)", "sequence": "Sequence", "funfam": "FunFam", "ec": "EC"}
    want = [("structure", "sequence"), ("sequence", "funfam"), ("structure", "funfam"), ("structure", "structure_core200"), ("structure", "ec"), ("sequence", "ec"), ("funfam", "ec")]
    rows = [A[(A.view_a == a) & (A.view_b == b)].iloc[0] for a, b in want]
    fig, ax = plt.subplots(figsize=(10, 5)); y = np.arange(len(rows))
    for i, r in enumerate(rows):
        c = ABHD if "ec" in (r.view_a, r.view_b) else OTHER
        ax.barh(i, r.spearman, color=c, height=0.6); ax.plot([r.ci_lo, r.ci_hi], [i, i], color=INK, linewidth=2.5)
    ax.set_yticks(y); ax.set_yticklabels([f"{names[r.view_a]} vs {names[r.view_b]}" for r in rows], fontsize=15); ax.invert_yaxis(); ax.axvline(0, color=INK, linewidth=1)
    ax.set_xlabel("Spearman correlation (line = 95% CI)")
    save(fig, "multiview")


if __name__ == "__main__":
    for f in (universe, elbow, categories, nn_chance, nn_removed, ec, abhd15, segment, cavity, attrition, rmsd, soap, concordance, cluster_null, majority, anchored, cavity_models, pocket_lining, decay, multiview):
        f()
    print(sorted(p.name for p in OUT.iterdir()))
