#!/usr/bin/env python3
"""
Study stages 2/3/7 input - UniProt annotations for the canonical accession of every universe gene:
active-site features (gold standard for catalytic-site calls), signal/transmembrane features, EC, catalytic activity
(with Rhea ids), family line and PDB cross-references.
Output: data/raw/uniprot/universe_annotations.tsv (frozen on first fetch; delete to refresh).
"""
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

OUT = C.RAW / "uniprot" / "universe_annotations.tsv"
FIELDS = "accession,gene_primary,length,ec,cc_catalytic_activity,cc_similarity,ft_act_site,ft_signal,ft_transmem,xref_pdb,cc_function,go_f"


def main():
    if OUT.exists():
        print(f"{OUT} exists; delete to refresh")
        return
    accs = pd.read_csv(C.DERIVED / "universe.csv").accession.tolist()
    chunks = []
    for i in range(0, len(accs), 30):               # UniProt rejects very long OR queries
        r = requests.get("https://rest.uniprot.org/uniprotkb/search", timeout=120, params={
            "query": "accession:(" + " OR ".join(accs[i:i + 30]) + ")", "fields": FIELDS, "format": "tsv", "size": 500})
        r.raise_for_status()
        lines = r.text.rstrip("\n").splitlines()
        chunks += lines if not chunks else lines[1:]
    OUT.write_text("\n".join(chunks) + "\n")
    df = pd.read_csv(OUT, sep="\t")
    print(f"{len(df)} annotated entries; {df['Active site'].notna().sum()} with active-site features -> {OUT}")


if __name__ == "__main__":
    main()
