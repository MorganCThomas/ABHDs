#!/usr/bin/env python3
"""
Fetch UniProt entries + sequences for the frozen set of human AB hydrolase superfamily accessions
(data/raw/uniprot/accessions.txt). `--from-query` re-runs the family query instead (see main()).
Outputs:
  data/raw/uniprot/abhd_entries.json    — list of entry dicts
  data/raw/uniprot/abhd_sequences.fasta — FASTA file with all sequences
"""

import json
import sys
from pathlib import Path

import requests

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

SEARCH_URL = "https://rest.uniprot.org/uniprotkb/search"
QUERY = C.UNIPROT_QUERY            # only used with --from-query
ACCESSIONS = C.ENTRIES_JSON.parent / "accessions.txt"   # the frozen study set
PAGE_SIZE = 500


def fetch_entries(QUERY) -> list[dict]:
    entries = []
    params = {
        "query": QUERY,
        "fields": "accession,id,protein_name,gene_primary,cc_catalytic_activity,cc_function,ec,reviewed",
        "format": "tsv",
        "size": PAGE_SIZE,
    }
    url = SEARCH_URL
    while url:
        resp = requests.get(url, params=params if url == SEARCH_URL else None, timeout=30)
        resp.raise_for_status()
        lines = resp.text.rstrip("\n").splitlines()
        for line in lines[1:]:  # skip header
            parts = line.split("\t")
            parts += [""] * (8 - len(parts))   # trailing empty TSV fields
            entries.append({
                "accession": parts[0].strip(),
                "entry_name": parts[1].strip(),
                "name": parts[2].strip(),
                "gene_name": parts[3].strip(),
                "catalytic_activity": parts[4].strip() if len(parts) > 4 else "",
                "function": parts[5].strip() if len(parts) > 5 else "",
                "ec": parts[6].strip() if len(parts) > 6 else "",
                "reviewed": parts[7].strip() if len(parts) > 7 else "",
            })
        link = resp.headers.get("Link", "")
        url = None
        if link:
            for part in link.split(","):
                part = part.strip()
                if 'rel="next"' in part:
                    url = part.split(";")[0].strip().strip("<>")
        params = None
    return entries


def fetch_fasta(QUERY) -> str:
    chunks = []
    params = {
        "query": QUERY,
        "format": "fasta",
        "size": PAGE_SIZE,
    }
    url = SEARCH_URL
    while url:
        resp = requests.get(url, params=params if url == SEARCH_URL else None, timeout=30)
        resp.raise_for_status()
        chunks.append(resp.text if resp.text.endswith("\n") else resp.text + "\n")
        link = resp.headers.get("Link", "")
        url = None
        if link:
            for part in link.split(","):
                part = part.strip()
                if 'rel="next"' in part:
                    url = part.split(";")[0].strip().strip("<>")
        params = None
    return "".join(chunks)


def fill_deleted(accs, fasta):
    """Entries deleted from UniProtKB are absent from FASTA search results; take their sequence from UniParc."""
    have = {line.split("|")[1] for line in fasta.splitlines() if line.startswith(">")}
    for a in accs:
        if a in have:
            continue
        info = requests.get(f"https://rest.uniprot.org/uniprotkb/{a}", timeout=30).json()
        upi = info.get("extraAttributes", {}).get("uniParcId")
        if not upi:
            raise SystemExit(f"No sequence found for {a}")
        seq = requests.get(f"https://rest.uniprot.org/uniparc/{upi}.fasta", timeout=30).text.splitlines()[1:]
        fasta += f">tr|{a}|{a}_HUMAN deleted from UniProtKB, sequence from UniParc {upi}\n" + "\n".join(seq) + "\n"
        print(f"  {a}: deleted from UniProtKB, sequence taken from UniParc {upi}")
    return fasta


def accession_query(accs):
    return "(" + " OR ".join(f"accession:{a}" for a in accs) + ")"


def main():
    """Default: fetch the frozen study set in accessions.txt (reproducible).
    --from-query: re-run the UniProt family query and write a NEW accessions.txt (the result drifts between
    UniProt releases: it returned 210 entries in 2025 and 378 in Nov 2026)."""
    from_query = "--from-query" in sys.argv
    C.ENTRIES_JSON.parent.mkdir(parents=True, exist_ok=True)
    if from_query:
        entries, fasta = fetch_entries(QUERY), fetch_fasta(QUERY)
        ACCESSIONS.write_text("\n".join(e["accession"] for e in entries) + "\n")
    else:
        accs = ACCESSIONS.read_text().split()
        entries, fasta = [], ""
        for i in range(0, len(accs), 50):
            q = accession_query(accs[i:i + 50])
            entries += fetch_entries(q)
            fasta += fetch_fasta(q)
        fasta = fill_deleted(accs, fasta)
        order = {a: i for i, a in enumerate(accs)}
        entries.sort(key=lambda e: order[e["accession"]])
    print(f"{len(entries)} entries")
    with open(C.ENTRIES_JSON, "w") as f:
        json.dump(entries, f, indent=2)
    sanitized = "\n".join(line.replace("/", "_") if line.startswith(">") else line for line in fasta.splitlines())
    with open(C.SEQUENCES_FASTA, "w") as f:
        f.write(sanitized + "\n")
    print(f"Saved {fasta.count('>')} sequences to {C.SEQUENCES_FASTA}")


if __name__ == "__main__":
    main()
