#!/usr/bin/env python3
"""
Label provenance (docs/redesign_plan.md, P1): how well-supported are the EC / Rhea / GO labels used as 'truth'?

Per reaction block of UniProt 'Catalytic activity': EC, Rhea IDs and the ECO evidence set. Reaction class: experimental (ECO:0000269 present),
else curator (0000305), similarity (0000250), rule (0000255), none. Per gene: reaction counts, complete EC levels, EC3 set, evidence class.
Outputs: data/derived/label_provenance.csv, label_provenance_summary.csv
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

CLASS_ORDER = ["experimental", "curator", "similarity", "rule", "none"]


def classify(eco):
    if "0000269" in eco: return "experimental"
    if "0000305" in eco: return "curator"
    if "0000250" in eco: return "similarity"
    if "0000255" in eco: return "rule"
    return "none"


def parse_blocks(text):
    out = []
    if not isinstance(text, str):
        return out
    for blk in text.split("CATALYTIC ACTIVITY:")[1:]:
        ec = re.search(r"EC=([0-9n.\-]+)", blk)
        ev = re.search(r"Evidence=\{([^}]*)\}", blk)
        eco = set(re.findall(r"ECO:(\d+)", ev.group(1))) if ev else set()
        out.append(dict(ec=ec.group(1) if ec else "", rhea=set(re.findall(r"RHEA:(\d+)", blk)), cls=classify(eco), eco=eco))
    return out


def ec_levels(ec):
    parts = ec.split(".")
    n = 0
    for p in parts:
        if p in ("-", "") or "n" in p: break
        n += 1
    return n


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv")
    ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
    rows = []
    for r in uni.itertuples():
        a = ann.loc[r.accession] if r.accession in ann.index else None
        ecs = [e.strip() for e in str(a["EC number"]).split(";") if e.strip() and e.strip() != "nan"] if a is not None else []
        complete = [e for e in ecs if ec_levels(e) >= 3]
        ec3 = sorted({".".join(e.split(".")[:3]) for e in complete})
        blocks = parse_blocks(a["Catalytic activity"]) if a is not None else []
        # reaction evidence for the EC3 classes the gene is labelled with
        rel = [b for b in blocks if b["ec"] and ".".join(b["ec"].split(".")[:3]) in ec3]
        classes = [b["cls"] for b in rel] or [b["cls"] for b in blocks]
        any_exp = any(c == "experimental" for c in classes)
        all_exp = bool(classes) and all(c == "experimental" for c in classes)
        if not ec3: cls = "no_complete_ec3"
        elif not blocks: cls = "no_record"
        elif all_exp: cls = "experimental_only"
        elif any_exp: cls = "mixed"
        else: cls = "non_experimental"
        go = str(a["Gene Ontology (molecular function)"]).split(";") if a is not None and isinstance(a["Gene Ontology (molecular function)"], str) else []
        pdb = str(a["PDB"]).split(";") if a is not None and isinstance(a["PDB"], str) else []
        rows.append(dict(gene=r.gene, accession=r.accession, abhd_named=r.abhd_named, ec_all=";".join(ecs), ec3=";".join(ec3), n_complete_ec3=len(ec3),
                         n_reactions=len(blocks), n_rhea=len(set().union(*[b["rhea"] for b in blocks])) if blocks else 0,
                         reaction_classes=";".join(sorted(set(b["cls"] for b in blocks))), evidence_any_experimental=any_exp, evidence_all_experimental=all_exp,
                         evidence_class=cls, n_go_mf=len([g for g in go if g.strip()]), n_pdb=len([p for p in pdb if p.strip()])))
    P = pd.DataFrame(rows); P.to_csv(C.DERIVED / "label_provenance.csv", index=False)
    e3 = P[P.n_complete_ec3 > 0]
    summ = {
        "genes": len(P), "genes_with_catalytic_blocks": int((P.n_reactions > 0).sum()), "genes_with_rhea": int((P.n_rhea > 0).sum()),
        "genes_complete_ec3": len(e3), "ec3_any_experimental": int(e3.evidence_any_experimental.sum()), "ec3_all_experimental": int(e3.evidence_all_experimental.sum()),
        "ec3_no_record": int((e3.evidence_class == "no_record").sum()), "genes_with_go_mf": int((P.n_go_mf > 0).sum()), "genes_with_pdb": int((P.n_pdb > 0).sum()),
        "abhd_complete_ec3": int(e3.abhd_named.sum()), "abhd_ec3_any_experimental": int(e3[e3.abhd_named].evidence_any_experimental.sum()),
    }
    pd.DataFrame([summ]).to_csv(C.DERIVED / "label_provenance_summary.csv", index=False)
    print(pd.Series(summ).to_string()); print(e3.evidence_class.value_counts().to_string())
    # acceptance: same EC3 set as function_loo (n=80)
    print("complete-EC3 genes:", len(e3), "(function_loo n=80)")


if __name__ == "__main__":
    main()
