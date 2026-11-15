#!/usr/bin/env python3
"""
Protein report: what the ABHD study knows about a protein and how it relates to the other 109.

  uv run python scripts/protein_report.py ABHD2                      # HTML report + 3D view in out/protein_reports/
  uv run python scripts/protein_report.py ABHD2 MGLL --orthologues   # also fetch/clean orthologues (needs network the first time)
  uv run python scripts/protein_report.py --batch genes.txt          # one gene per line; adds out/protein_reports/summary.csv
  uv run python scripts/protein_report.py --all --no-structure       # every gene in the universe (fast)

Needs the derived tables of the study (run `make study` or have data/derived present). Outputs go to --out (default out/protein_reports, git-ignored).
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd.atlas import Atlas
from abhd.protein_report import build_report


def summary_row(A, res):
    i = res["info"]; T = res["neighbours"]; g = i["gene"]
    row = dict(gene=g, accession=i["accession"], abhd_named=i["abhd_named"], category=i.get("category"), ec3=";".join(i["ec3"]), evidence_class=i["evidence_class"])
    for v in ("structure", "sequence", "funfam"):
        if v in A.available_views():
            n = A.neighbours(g, v, k=3)
            row[f"nearest_{v}"] = ", ".join(f"{r.gene} ({r.identity:.0%})" for r in n.itertuples())
    p = A.predict_ec(g, "structure")
    row.update(predicted_ec3=("gated: no nucleophile" if p.get("gated") else p.get("ec3")), prediction_neighbour_identity=p.get("identity"), prediction_identity_tier=p.get("tier"))
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("genes", nargs="*", help="gene symbols (case-insensitive) or UniProt accessions in the universe")
    ap.add_argument("--batch", help="text file with one gene per line")
    ap.add_argument("--all", action="store_true", help="every gene in the universe")
    ap.add_argument("--orthologues", action="store_true", help="add the orthologue conservation section (fetches UniProt the first time)")
    ap.add_argument("--no-structure", action="store_true", help="skip the superposed 3D view files")
    ap.add_argument("--plotly-cdn", action="store_true", help="load plotly.js from a CDN instead of a local plotly.min.js next to the reports")
    ap.add_argument("--k", type=int, default=6, help="neighbours per view in the table (default 6)")
    ap.add_argument("--out", default="out/protein_reports")
    a = ap.parse_args()
    A = Atlas()
    genes = list(a.genes)
    if a.batch:
        genes += [x.strip() for x in Path(a.batch).read_text().splitlines() if x.strip() and not x.startswith("#")]
    if a.all:
        genes = list(A.uni.gene)
    if not genes:
        ap.error("give at least one gene, --batch or --all")
    rows = []
    resolved = set()
    for g in genes:
        try:
            resolved.add(A.gene(A.acc(g)))
        except KeyError:
            pass
    known = resolved | {p.stem for p in Path(a.out).glob("*.html") if not p.stem.endswith("_view")}      # reports that exist or are being built: clickable in the plots
    for g in genes:
        try:
            res = build_report(A, g, a.out, orthologues=a.orthologues, structure=not a.no_structure, k=a.k, known_reports=known, plotly_cdn=a.plotly_cdn)
        except KeyError as e:
            print(f"skip: {e.args[0]}", file=sys.stderr); continue
        print(f"{res['gene']}: {res['path']}", flush=True)
        rows.append(summary_row(A, res))
    if len(rows) > 1:
        p = Path(a.out) / "summary.csv"; pd.DataFrame(rows).to_csv(p, index=False); print("summary:", p)


if __name__ == "__main__":
    main()
