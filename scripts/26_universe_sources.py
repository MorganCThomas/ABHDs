#!/usr/bin/env python3
"""
Third universe source (docs/redesign_plan.md, P3): reviewed human UniProt entries annotated with InterPro IPR029058 (alpha/beta hydrolase fold)
or any of several Pfam alpha/beta-hydrolase families. Reports list overlap and candidate genes not in the 110; no population-size estimate.
Cached in data/raw/uniprot/third_source_<name>.tsv with the fetch date in third_source_fetched.txt.
Outputs: data/derived/universe_sources.csv (overlap patterns), universe_candidates.csv
"""
import sys
from datetime import date
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
PFAM = ["PF00561", "PF12697", "PF12146", "PF00135", "PF00326", "PF01738", "PF07859", "PF02230", "PF03096", "PF06441", "PF02450", "PF00975"]
FIELDS = "accession,gene_primary,protein_name,reviewed,length,xref_pfam,xref_interpro"
QUERIES = {"interpro_IPR029058": "(xref:interpro-IPR029058) AND (organism_id:9606) AND (reviewed:true)",
           "pfam_aby": "(" + " OR ".join(f"xref:pfam-{p}" for p in PFAM) + ") AND (organism_id:9606) AND (reviewed:true)"}


def fetch(name, q):
    p = C.RAW / "uniprot" / f"third_source_{name}.tsv"
    if not p.exists():
        parts, r = [], requests.get(SEARCH, params={"query": q, "format": "tsv", "size": 500, "fields": FIELDS}, timeout=90)
        while True:
            r.raise_for_status(); parts.append(r.text if not parts else r.text.split("\n", 1)[1])
            nxt = r.links.get("next", {}).get("url")
            if not nxt: break
            r = requests.get(nxt, timeout=90)
        p.write_text("".join(parts))
        (C.RAW / "uniprot" / "third_source_fetched.txt").write_text(f"fetched {date.today().isoformat()} from the UniProt REST API (current release)\n")
    return pd.read_csv(p, sep="\t")


def main():
    uni = pd.read_csv(C.DERIVED / "universe.csv")
    ent = pd.read_csv(C.DERIVED / "universe_entries.csv").set_index("accession")
    acc2gene = {a: g for a, g in ent.assigned_gene.items() if isinstance(g, str) and g}
    third = {}
    for name, q in QUERIES.items():
        d = fetch(name, q); third[name] = d
        print(f"{name}: {len(d)} reviewed human entries")
    allthird = pd.concat(third.values()).drop_duplicates("Entry")
    allthird["gene"] = allthird.Entry.map(acc2gene).fillna(allthird["Gene Names (primary)"])
    genes3 = set(allthird.gene.dropna())
    u = uni.set_index("gene")
    l1 = set(u.index[(u.in_uniprot_2025 | u.in_uniprot_2026)]); l2 = set(u.index[u.in_funfam]); l3 = genes3
    universe = set(u.index)
    pats = {}
    for g in l1 | l2 | l3:
        key = ("UniProt" if g in l1 else "-", "FunFam" if g in l2 else "-", "InterPro/Pfam" if g in l3 else "-")
        pats.setdefault(key, []).append(g)
    S = pd.DataFrame([dict(uniprot=k[0], funfam=k[1], third=k[2], n_genes=len(v), in_110=sum(g in universe for g in v)) for k, v in sorted(pats.items())])
    S.to_csv(C.DERIVED / "universe_sources.csv", index=False)
    print(S.to_string(index=False)); print("union", len(l1 | l2 | l3), "| first two lists reproduce:", len(l1 | l2), "(110)", len(l1 & l2), "both")
    cand = allthird[~allthird.gene.isin(universe)].copy()
    cand = cand[["Entry", "gene", "Protein names", "Length", "Pfam", "InterPro"]].sort_values("gene")
    cand.to_csv(C.DERIVED / "universe_candidates.csv", index=False)
    print(f"{len(cand)} candidate genes in the third source but not in the 110"); print(cand[["Entry", "gene", "Protein names"]].head(60).to_string(index=False))
    print("110 genes absent from the third source:", len(universe - l3), "of 110; ABHD-named absent:", sum(g.startswith("ABHD") for g in universe - l3))


if __name__ == "__main__":
    main()
