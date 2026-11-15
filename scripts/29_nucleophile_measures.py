#!/usr/bin/env python3
"""
Nucleophile calls as three measures on the 31 gold genes (docs/redesign_plan.md, P6, reduced; no latent-class model):
  elbow     the elbow call is exactly the annotated nucleophile
  motif     the legacy sequence motif matches at the annotated nucleophile
  transfer  leave-one-out structural transfer: the residue aligned (TM-align, CATH domains) to the annotated nucleophile of the structurally nearest OTHER
            gold gene (nearest by core TM distance) is the annotated nucleophile
Stratified by the evidence of the UniProt nucleophile annotation (experimental = ECO:0000269; curator 0000305; similarity 0000250; rule 0000255).
Also: motif hit counts per no-nucleophile gene. Outputs: data/derived/nucleophile_measures.csv, nucleophile_measures_summary.csv
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.structops import align


def evidence_class(eco):
    if "0000269" in eco: return "experimental"
    if "0000305" in eco: return "curator"
    if "0000250" in eco: return "similarity"
    if "0000255" in eco: return "rule"
    return "none"


def main():
    b = pd.read_csv(C.DERIVED / "benchmark_stage3_genes.csv").set_index("accession")
    ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv").set_index("accession")
    D = pd.read_csv(C.DERIVED / "dist_tmalign_core200.csv", index_col=0)
    gold = list(b.index)
    rows = []
    for a in gold:
        site = str(ann.loc[a, "Active site"])
        ev = re.search(r'ACT_SITE (\d+); /note="Nucleophile[^"]*"; /evidence="([^"]+)"', site)
        eco = set(re.findall(r"ECO:(\d+)", ev.group(2))) if ev else set()
        d = D.loc[a, [g for g in gold if g != a]].astype(float).dropna().sort_values()
        h = d.index[0]
        (ia, _, sa), (ib, _, sb), X, pairs, r = align(a, h)
        j = int(np.where(ib == int(b.loc[h, "gold_nuc"]))[0][0]); i2j = {j2: i2 for i2, j2 in pairs}
        pos = int(ia[i2j[j]]) if j in i2j else None
        gn = int(b.loc[a, "gold_nuc"])
        rows.append(dict(gene=b.loc[a, "gene"], accession=a, gold_nuc=gn, evidence=evidence_class(eco), nearest_gold_gene=b.loc[h, "gene"], transfer_pos=pos,
                         elbow=bool(b.loc[a, "elbow_hit"]), motif=bool(b.loc[a, "regex_at_gold"]), transfer=bool(pos == gn), transfer_within1=bool(pos is not None and abs(pos - gn) <= 1)))
    T = pd.DataFrame(rows); T.to_csv(C.DERIVED / "nucleophile_measures.csv", index=False)
    S = []
    for name, sub in (("all", T), *[(k, g) for k, g in T.groupby("evidence")]):
        S.append(dict(stratum=name, n=len(sub), elbow=int(sub.elbow.sum()), motif=int(sub.motif.sum()), transfer_exact=int(sub.transfer.sum()), transfer_within1=int(sub.transfer_within1.sum())))
    S = pd.DataFrame(S); S.to_csv(C.DERIVED / "nucleophile_measures_summary.csv", index=False)
    pd.set_option("display.width", 200); print(S.to_string(index=False))
    print("elbow 31/31 and motif 26/31 check:", int(T.elbow.sum()), int(T.motif.sum()))
    print(T[~T.transfer][["gene", "gold_nuc", "evidence", "nearest_gold_gene", "transfer_pos"]].to_string(index=False))
    nn = cs[cs.category.astype(str).str.startswith("no nucleophile")]
    print("motif hits in the 14 no-nucleophile genes (per gene): median", nn.regex_hits.median() if "regex_hits" in nn else "n/a", "| all genes mean", round(cs.regex_hits.mean(), 2))


if __name__ == "__main__":
    main()
