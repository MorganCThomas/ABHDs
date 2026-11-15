#!/usr/bin/env python3
"""
Query the RCSB PDB for protein sequence similarity.

Searches against experimentally solved structures (X-ray or CryoEM) and
returns the top N hits with identity, resolution, organism, and deposit date.

Usage:
    python pdb_sequence_search.py SEQUENCE [-n N]
    python pdb_sequence_search.py -f sequence.fasta [-n N]
"""

import argparse
import concurrent.futures
import sys

import requests

SEARCH_URL = "https://search.rcsb.org/rcsbsearch/v2/query"
DATA_URL = "https://data.rcsb.org/rest/v1/core"


def parse_fasta(path: str) -> str:
    """Return the first sequence from a FASTA file (no headers)."""
    with open(path) as fh:
        lines = [l.strip() for l in fh if not l.startswith(">") and l.strip()]
    return "".join(lines)


def search_by_sequence(sequence: str, top_n: int) -> list[dict]:
    """POST a sequence similarity query to the RCSB search API."""
    query = {
        "query": {
            "type": "group",
            "logical_operator": "and",
            "nodes": [
                {
                    "type": "terminal",
                    "service": "sequence",
                    "parameters": {
                        "evalue_cutoff": 10000,
                        "identity_cutoff": 0,
                        "sequence_type": "protein",
                        "value": sequence,
                    },
                },
                {
                    "type": "terminal",
                    "service": "text",
                    "parameters": {
                        "attribute": "exptl.method",
                        "operator": "in",
                        "value": ["X-RAY DIFFRACTION", "ELECTRON MICROSCOPY"],
                        "negation": False,
                    },
                },
            ],
        },
        "request_options": {
            "scoring_strategy": "sequence",
            "results_verbosity": "verbose",
            "paginate": {"start": 0, "rows": top_n},
        },
        "return_type": "polymer_entity",
    }

    resp = requests.post(SEARCH_URL, json=query, timeout=60)
    resp.raise_for_status()
    if resp.status_code == 204 or not resp.text.strip():
        return []
    return resp.json().get("result_set", [])


def fetch_hit_details(pdb_id: str, entity_id: str) -> dict:
    """Fetch entry- and entity-level metadata from the RCSB Data API."""
    entry_resp = requests.get(f"{DATA_URL}/entry/{pdb_id}", timeout=30)
    entity_resp = requests.get(
        f"{DATA_URL}/polymer_entity/{pdb_id}/{entity_id}", timeout=30
    )
    entry_resp.raise_for_status()
    entity_resp.raise_for_status()

    entry = entry_resp.json()
    entity = entity_resp.json()

    method = entry.get("exptl", [{}])[0].get("method", "N/A")

    resolution_list = entry.get("rcsb_entry_info", {}).get("resolution_combined", [])
    resolution = f"{resolution_list[0]:.2f}" if resolution_list else "N/A"

    raw_date = entry.get("rcsb_accession_info", {}).get("deposit_date", "")
    deposit_date = raw_date[:10] if raw_date else "N/A"

    source_organisms = entity.get("rcsb_entity_source_organism", [])
    organism = source_organisms[0].get("scientific_name", "N/A") if source_organisms else "N/A"

    chains = entity.get("entity_poly", {}).get("pdbx_strand_id", "N/A")

    return {
        "method": method,
        "resolution": resolution,
        "deposit_date": deposit_date,
        "organism": organism,
        "chains": chains,
    }


def _fetch_one(args: tuple) -> dict | None:
    pdb_id, entity_id, identifier, score = args
    try:
        details = fetch_hit_details(pdb_id, entity_id)
    except Exception as exc:
        print(f"  Warning: could not fetch details for {identifier}: {exc}", file=sys.stderr)
        details = {
            "method": "N/A",
            "resolution": "N/A",
            "deposit_date": "N/A",
            "organism": "N/A",
            "chains": "N/A",
        }
    return {"pdb_id": pdb_id, "entity_id": entity_id, "identity": score, **details}


def print_table(rows: list[dict]) -> None:
    headers = [
        "PDB ID", "Entity", "Chains", "Identity (%)",
        "Method", "Resolution (Å)", "Organism", "Deposited",
    ]

    table = []
    for r in rows:
        table.append([
            r["pdb_id"],
            r["entity_id"],
            r["chains"],
            f"{r['identity'] * 100:.1f}",
            r["method"],
            r["resolution"],
            r["organism"],
            r["deposit_date"],
        ])

    col_widths = [
        max(len(headers[i]), max(len(row[i]) for row in table))
        for i in range(len(headers))
    ]

    def fmt_row(row):
        return "  ".join(cell.ljust(col_widths[i]) for i, cell in enumerate(row))

    separator = "  ".join("-" * w for w in col_widths)
    print(fmt_row(headers))
    print(separator)
    for row in table:
        print(fmt_row(row))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Search the PDB for structurally solved proteins similar to a query sequence."
    )
    seq_group = parser.add_mutually_exclusive_group(required=True)
    seq_group.add_argument(
        "sequence",
        nargs="?",
        help="Protein sequence in single-letter amino acid code.",
    )
    seq_group.add_argument(
        "-f", "--fasta",
        metavar="FILE",
        help="Path to a FASTA file (first sequence used).",
    )
    parser.add_argument(
        "-n", "--top-n",
        type=int,
        default=10,
        metavar="N",
        help="Number of top hits to return (default: 10).",
    )
    args = parser.parse_args()

    sequence = parse_fasta(args.fasta) if args.fasta else args.sequence.strip().upper()

    if not sequence:
        parser.error("No sequence provided.")

    print(f"Querying PDB for top {args.top_n} sequence similarity hits (X-ray / CryoEM)...\n")

    hits = search_by_sequence(sequence, args.top_n)

    if not hits:
        print("No hits found.")
        return

    fetch_args = []
    for hit in hits:
        identifier = hit["identifier"]   # e.g. "4HHB_1"
        # Extract sequence identity from the verbose match context
        identity = None
        for svc in hit.get("services", []):
            if svc.get("service_type") == "sequence":
                for node in svc.get("nodes", []):
                    ctx = node.get("match_context", [])
                    if ctx:
                        identity = ctx[0].get("sequence_identity")
                        break
        if identity is None:
            identity = 0.0
        parts = identifier.split("_")
        pdb_id, entity_id = parts[0], parts[1]
        fetch_args.append((pdb_id, entity_id, identifier, identity))

    print(f"Fetching details for {len(fetch_args)} hits...\n")

    with concurrent.futures.ThreadPoolExecutor(max_workers=16) as pool:
        rows = list(pool.map(_fetch_one, fetch_args))

    rows.sort(key=lambda r: r["identity"], reverse=True)
    print_table(rows)


if __name__ == "__main__":
    main()
