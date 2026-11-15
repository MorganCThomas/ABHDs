#!/usr/bin/env python3
"""
Fetch AlphaFold predicted structures for entries in data/raw/uniprot/abhd_entries.json.
Structures are saved as PDB files in data/raw/structures/ (existing files are skipped).
Accessions with no AlphaFold entry are reported at the end.
`--universe` fetches the canonical accession of every gene in data/derived/universe.csv instead.
"""

import json
import os
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction"
ENTRIES_FILE = C.ENTRIES_JSON
OUT_DIR = str(C.STRUCTURES)


def fetch_structure(accession: str, out_dir: str) -> bool:
    """Download PDB for accession. Returns True if successful, False if not found."""
    out_path = os.path.join(out_dir, f"{accession}.pdb")
    if os.path.exists(out_path):
        return True

    resp = requests.get(f"{ALPHAFOLD_API}/{accession}", timeout=30)
    if resp.status_code == 404:
        return False
    resp.raise_for_status()

    data = resp.json()
    if not data:
        return False

    pdb_url = data[0].get("pdbUrl")
    if not pdb_url:
        return False

    pdb_resp = requests.get(pdb_url, timeout=60)
    pdb_resp.raise_for_status()

    with open(out_path, "w") as f:
        f.write(pdb_resp.text)
    return True


def main():
    universe = "--universe" in sys.argv    # study stage 2: canonical accessions of the merged universe (stage 1)
    if universe:
        import csv
        with open(C.DERIVED / "universe.csv") as f:
            entries = [{"accession": r["accession"], "entry_name": r["entry_name"]} for r in csv.DictReader(f)]
    else:
        with open(ENTRIES_FILE) as f:
            entries = json.load(f)

    os.makedirs(OUT_DIR, exist_ok=True)

    missing = []
    for i, entry in enumerate(entries, 1):
        accession = entry["accession"]
        entry_name = entry.get("entry_name", accession)
        print(f"[{i}/{len(entries)}] {entry_name} ({accession})", end=" ... ", flush=True)
        found = fetch_structure(accession, OUT_DIR)
        if found:
            print("ok")
        else:
            print("not found")
            missing.append(entry_name)

    log = [f"Done. {len(entries) - len(missing)}/{len(entries)} structures saved to {OUT_DIR}/"]
    if missing:
        log.append(f"No AlphaFold entry for: {', '.join(missing)}")
    open(os.path.join(OUT_DIR, "fetch_log_universe.txt" if universe else "fetch_log.txt"), "w").write("\n".join(log) + "\n")
    print(f"\nDone. {len(entries) - len(missing)}/{len(entries)} structures saved to {OUT_DIR}/")
    if missing:
        print(f"No AlphaFold entry for: {', '.join(missing)}")


if __name__ == "__main__":
    main()
