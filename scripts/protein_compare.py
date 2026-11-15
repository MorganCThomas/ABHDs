#!/usr/bin/env python3
"""
Compare two proteins of the ABHD study: similarity in every view, aligned catalytic residues, second oxyanion donors, pocket-lining residues.

  uv run python scripts/protein_compare.py ABHD2 MGLL
  uv run python scripts/protein_compare.py ABHD2 ABHD6 --out out/compare_ABHD2_ABHD6.md
"""
import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd.atlas import Atlas
from abhd.compare import compare, to_markdown


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("a"); ap.add_argument("b"); ap.add_argument("--out", help="write the Markdown here as well as printing it")
    args = ap.parse_args()
    md = to_markdown(compare(Atlas(), args.a, args.b))
    print(md)
    if args.out:
        Path(args.out).parent.mkdir(parents=True, exist_ok=True); Path(args.out).write_text(md + "\n")


if __name__ == "__main__":
    main()
