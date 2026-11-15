#!/usr/bin/env python3
"""
Study stages 3-4 - catalytic-site calls and categories for every universe gene.

Per gene: nucleophile elbow (any residue type, so a replaced nucleophile is still seen), His/acid triad with the
chain-order constraint nucleophile < acid < His, oxyanion-hole candidate. Scored against UniProt "Active site"
annotations (Nucleophile / Charge relay system) and against the legacy regex. Threshold sensitivity is computed from the
raw table without re-reading structures.
Outputs (data/derived/): elbow_table.csv, catalytic_sites.csv, benchmark_stage3.csv, sensitivity_stage3.csv,
attrition_stage3.csv; figures/stage3_*.html, figures/stage4_categories.html
"""
import itertools
import re
import sys
from pathlib import Path

import pandas as pd
import plotly.graph_objects as go

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.catalytic import apply_elbow, his_near_ca, scan_elbow_table, triad_and_oxyanion

NUC_AA = set("SCTD")
LEGACY = re.compile(C.LEGACY_REGEX)
ONE2THREE = {"S": "Ser", "C": "Cys", "T": "Thr", "D": "Asp", "E": "Glu", "H": "His", "N": "Asn", "G": "Gly", "A": "Ala"}


def gold_sites(text):
    """UniProt ACT_SITE string -> (nucleophile positions, other catalytic positions)."""
    nuc, other = [], []
    for pos, note in re.findall(r'ACT_SITE (\d+); /note="([^"]+)"', str(text)):
        (nuc if note.startswith("Nucleophile") else other).append(int(pos))
    return nuc, other


def build_table(qc):
    rows = []
    for r in qc.itertuples():
        if r.flag == "no_model":
            continue
        rows += scan_elbow_table(C.STRUCTURES / f"{r.accession}.pdb", lo=r.dom_start, hi=r.dom_end)
    return pd.DataFrame(rows)


def main():
    qc = pd.read_csv(C.DERIVED / "structure_qc.csv")
    uni = pd.read_csv(C.DERIVED / "universe.csv").set_index("accession")
    ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
    tab = build_table(qc)
    tab.to_csv(C.DERIVED / "elbow_table.csv", index=False)
    tab["elbow"] = apply_elbow(tab)
    el = tab[tab.elbow].copy()
    el["nuc_type"] = el.aa.isin(NUC_AA)

    rows = []
    for q in qc.itertuples():
        acc = q.accession
        nuc_gold, oth_gold = gold_sites(ann.loc[acc, "Active site"]) if acc in ann.index else ([], [])
        sub = el[el.accession == acc] if q.flag != "no_model" else el.iloc[0:0]
        # annotate every elbow candidate with a triad so the primary call can prefer a complete site
        info = []
        for e in sub.itertuples():
            pdb = C.STRUCTURES / f"{acc}.pdb"
            t = triad_and_oxyanion(pdb, e.resnum) if e.nuc_type else {"his_near_ca": his_near_ca(pdb, e.resnum)}
            info.append((e, t))
        def rank(x):
            e, t = x
            return (not e.nuc_type, not (t.get("his") and t.get("order_ok")), e.aa == "G", not t.get("his_near_ca", False), -e.plddt)
        info.sort(key=rank)
        row = dict(gene=q.gene, accession=acc, abhd_named=q.abhd_named, source=q.source, funfam=q.funfam,
                   n_elbow=len(sub), gold_nuc=";".join(map(str, nuc_gold)), gold_other=";".join(map(str, oth_gold)))
        if info:
            e, t = info[0]
            row.update(elbow_resnum=e.resnum, elbow_aa=e.aa, elbow_context=e.context, phi=round(e.phi, 1),
                       psi=round(e.psi, 1), plddt=round(e.plddt, 1), **{k: v for k, v in t.items() if k != 'his_near_ca'})
        # legacy regex: any 5-mer in the domain model that matches (as in the original notebook)
        rows.append(row)
    cs = pd.DataFrame(rows)

    # regex baseline: does the regex match at the gold nucleophile, and how many S/C/T residues does it flag per gene
    reg = {}
    for q in qc.itertuples():
        if q.flag == "no_model":
            continue
        from Bio.PDB import PDBParser, PPBuilder
        from abhd.elbow import THREE_TO_ONE
        s = PDBParser(QUIET=True).get_structure(q.accession, str(C.STRUCTURES / f"{q.accession}.pdb"))
        hits = []
        for pp in PPBuilder().build_peptides(s):
            res = list(pp)
            for i in range(2, len(res) - 2):
                if res[i].get_resname() in ("SER", "CYS", "THR") and q.dom_start <= res[i].get_id()[1] <= q.dom_end:
                    ctx = "".join(THREE_TO_ONE.get(x.get_resname(), "X") for x in res[i - 2:i + 3])
                    if LEGACY.match(ctx):
                        hits.append(res[i].get_id()[1])
        reg[q.accession] = hits
    cs["regex_hits"] = cs.accession.map(lambda a: len(reg.get(a, [])))
    cs["regex_at_gold"] = cs.apply(lambda r: bool(r.gold_nuc) and int(r.gold_nuc.split(";")[0]) in reg.get(r.accession, []), axis=1)

    # ---- categories (stage 4) ----
    def category(r):
        if pd.isna(r.elbow_resnum):
            return "no elbow found"
        aa = r.elbow_aa
        if aa not in NUC_AA:
            return f"no nucleophile: {ONE2THREE.get(aa, aa)} at elbow"
        n = ONE2THREE[aa]
        if pd.isna(r.his):
            return f"{n}, no His in range"
        if pd.isna(r.acid):
            return f"{n}-His, no acid"
        a = "Asp" if r.acid_type == "ASP" else "Glu"
        return f"{n}-His-{a}" + ("" if r.order_ok else " (non-canonical order)")
    cs["category"] = cs.apply(category, axis=1)
    cs.to_csv(C.DERIVED / "catalytic_sites.csv", index=False)

    # ---- benchmark ----
    g = cs[cs.gold_nuc.ne("") & cs.gold_nuc.notna()].copy()
    g["gold_first"] = g.gold_nuc.str.split(";").str[0].astype(int)
    g["elbow_hit"] = g.elbow_resnum == g.gold_first
    # an elbow call anywhere equal to gold (not only the primary)
    anyhit = {a: set(el[el.accession == a].resnum) for a in g.accession}
    g["any_elbow_hit"] = g.apply(lambda r: r.gold_first in anyhit[r.accession], axis=1)
    bench = dict(
        n_genes=len(cs), n_with_gold_nucleophile=len(g),
        elbow_primary_match=int(g.elbow_hit.sum()), elbow_any_match=int(g.any_elbow_hit.sum()),
        regex_match=int(g.regex_at_gold.sum()),
        elbow_calls_per_gene_median=float(cs.n_elbow.median()),
        regex_hits_per_gene_median=float(cs.regex_hits.median()),
        n_no_elbow=int((cs.n_elbow == 0).sum()))
    # triad: His/acid vs annotated charge-relay residues
    def triad_match(r):
        if pd.isna(r.his) or not r.gold_other:
            return None
        gold = set(map(int, r.gold_other.split(";")))
        return int(r.his) in gold and int(r.acid) in gold if not pd.isna(r.acid) else int(r.his) in gold
    g["triad_matches_gold"] = g.apply(triad_match, axis=1)
    g["his_matches_gold"] = g.apply(lambda r: None if pd.isna(r.his) or not r.gold_other else int(r.his) in set(map(int, r.gold_other.split(";"))), axis=1)
    bench_his = int(g.his_matches_gold.fillna(False).sum()); bench_his_n = int(g.his_matches_gold.notna().sum())
    bench["triad_with_call_and_gold"] = int(g.triad_matches_gold.notna().sum())
    bench["triad_matches_gold"] = int(g.triad_matches_gold.fillna(False).sum())
    bench["his_matches_gold"] = bench_his
    pd.Series(bench).rename("value").to_csv(C.DERIVED / "benchmark_stage3.csv")
    g[["gene", "accession", "gold_nuc", "gold_other", "elbow_resnum", "elbow_aa", "elbow_hit", "any_elbow_hit",
       "regex_at_gold", "his", "acid", "triad_matches_gold", "his_matches_gold", "category"]].to_csv(
        C.DERIVED / "benchmark_stage3_genes.csv", index=False)

    # ---- sensitivity: vary thresholds on the raw table ----
    gold_pos = {r.accession: int(r.gold_nuc.split(";")[0]) for r in cs.itertuples() if isinstance(r.gold_nuc, str) and r.gold_nuc}
    sens = []
    grid = itertools.product([(30, 90), (40, 75), (20, 100)], [(-170, -80), (-150, -95)], [60, 70, 80], [True, False])
    for phi, psi, pl, ss in grid:
        m = tab[apply_elbow(tab, phi=phi, psi=psi, plddt_min=pl, need_ss=ss)]
        calls = m.groupby("accession").resnum.apply(set)
        hit = sum(gold_pos[a] in calls.get(a, set()) for a in gold_pos)
        genes_with = sum(1 for a in calls.index)
        nuc_calls = m[m.aa.isin(NUC_AA)].groupby("accession").size()
        sens.append(dict(phi=str(phi), psi=str(psi), plddt_min=pl, need_ss=ss, genes_with_call=genes_with,
                         gold_recovered=hit, gold_total=len(gold_pos),
                         median_calls_per_gene=float(nuc_calls.median()) if len(nuc_calls) else 0,
                         genes_with_multiple_nuc_calls=int((nuc_calls > 1).sum())))
    pd.DataFrame(sens).to_csv(C.DERIVED / "sensitivity_stage3.csv", index=False)

    att = cs.category.value_counts().rename_axis("category").reset_index(name="n_genes")
    att.to_csv(C.DERIVED / "attrition_stage3.csv", index=False)

    # ---- figures ----
    C.FIGURES.mkdir(exist_ok=True)
    gp = set(zip(g.accession, g.gold_first))
    allp = tab[tab.aa.isin(NUC_AA)].copy()
    allp["is_gold"] = [(a, n) in gp for a, n in zip(allp.accession, allp.resnum)]
    fig = go.Figure()
    for flag, name, col, sz in ((False, "other S/T/C/D with phi>0", "lightgrey", 5), (True, "UniProt nucleophile", "crimson", 8)):
        d = allp[allp.is_gold == flag]
        fig.add_scatter(x=d.phi, y=d.psi, mode="markers", name=name, marker=dict(color=col, size=sz),
                        text=d.accession + " " + d.aa + d.resnum.astype(str) + " " + d.context,
                        hovertemplate="%{text}<extra></extra>")
    fig.add_shape(type="rect", x0=C.ELBOW_PHI[0], x1=C.ELBOW_PHI[1], y0=C.ELBOW_PSI[0], y1=C.ELBOW_PSI[1],
                  line=dict(dash="dot"))
    fig.update_layout(template="simple_white", width=800, height=550, title="Stage 3: Ramachandran position of annotated nucleophiles",
                      xaxis_title="phi", yaxis_title="psi")
    fig.write_html(C.FIGURES / "stage3_elbow_phi_psi.html")
    cat = pd.crosstab(cs.category, cs.abhd_named).reindex(cs.category.value_counts().index)
    fig = go.Figure()
    for named, name, col in ((True, "ABHD-named", "crimson"), (False, "other name", "lightgrey")):
        fig.add_bar(y=cat.index, x=cat.get(named, 0), orientation="h", name=name, marker_color=col)
    fig.update_layout(barmode="stack", template="simple_white", width=850, height=450, yaxis=dict(autorange="reversed"),
                      title="Stage 4: catalytic-site category by gene name", xaxis_title="genes")
    fig.write_html(C.FIGURES / "stage4_categories.html")
    print(pd.Series(bench).to_string()); print(att.to_string(index=False))


if __name__ == "__main__":
    main()
