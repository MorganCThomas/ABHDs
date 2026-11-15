"""Map every UniProt accession in the FunFam metadata to entry name + gene (UniProt REST).

Existing entries in data/derived/ff_uniprot_map.json are reused; only missing accessions hit the network.
"""
import json
import sys
from pathlib import Path

import pandas as pd
import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C


def lookup(acc):
    r = requests.get(f"https://rest.uniprot.org/uniprotkb/{acc}",
                     params={"fields": "accession,id,gene_primary", "format": "tsv"}, timeout=30)
    r.raise_for_status()
    _, name, gene = r.text.strip().splitlines()[1].split("\t")
    return {"name": name, "gene": gene}


def main():
    m = json.load(open(C.FF_UNIPROT_MAP)) if C.FF_UNIPROT_MAP.exists() else {}
    accs = {a for s in pd.read_csv(C.FUNFAM_META).human_uniprots for a in s.split(";")}
    for a in sorted(accs - set(m)):
        m[a] = lookup(a)
    json.dump(m, open(C.FF_UNIPROT_MAP, "w"))
    print(f"{len(accs)} FunFam accessions mapped ({len(accs - set(m))} missing)")


if __name__ == "__main__":
    main()
