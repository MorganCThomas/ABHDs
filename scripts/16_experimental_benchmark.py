#!/usr/bin/env python3
"""
Study stage 3 addendum - do experimental structures give the same elbow as the AlphaFold models?

For every universe gene with a PDB cross-reference (UniProt), download up to MAX_PDB entries (PDB format), apply the same
elbow definition (no pLDDT filter) to all chains, and ask whether an elbow with the same 5-residue context and the same
residue type as the AlphaFold call is found.
Output: data/derived/experimental_benchmark.csv; experimental PDBs are cached in data/raw/pdb/.
"""
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.catalytic import apply_elbow, scan_elbow_table

MAX_PDB = 3
PDB_DIR = C.RAW / "pdb"


def fetch(pid):
    f = PDB_DIR / f"{pid}.pdb"
    if f.exists():
        return f
    r = requests.get(f"https://files.rcsb.org/download/{pid}.pdb", timeout=60)
    if r.status_code != 200:
        return None
    f.write_text(r.text)
    return f


def main():
    PDB_DIR.mkdir(parents=True, exist_ok=True)
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv").set_index("accession")
    ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
    rows = []
    for acc, r in cs.iterrows():
        ids = [x for x in str(ann.loc[acc, "PDB"]).split(";") if x and x != "nan"][:MAX_PDB]
        if not ids or pd.isna(r.elbow_context):
            continue
        tested = hit = 0
        for pid in ids:
            f = fetch(pid)
            if f is None:
                continue
            t = pd.DataFrame(scan_elbow_table(f, accession=pid))
            if t.empty:
                tested += 1; continue
            t = t[apply_elbow(t, plddt_min=0)]
            tested += 1
            hit += int(((t.context == r.elbow_context) & (t.aa == r.elbow_aa)).any())
        rows.append(dict(gene=r.gene, accession=acc, af_context=r.elbow_context, af_aa=r.elbow_aa, n_pdb_tested=tested, n_pdb_same_elbow=hit))
    df = pd.DataFrame(rows); df.to_csv(C.DERIVED / "experimental_benchmark.csv", index=False)
    ok = df[df.n_pdb_tested > 0]
    print(f"{len(ok)} genes with experimental structures; {int((ok.n_pdb_same_elbow > 0).sum())} reproduce the AlphaFold elbow "
          f"({int(ok.n_pdb_same_elbow.sum())}/{int(ok.n_pdb_tested.sum())} entries)")
    print(df[df.n_pdb_same_elbow == 0].to_string())


if __name__ == "__main__":
    main()
