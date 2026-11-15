"""Find the nucleophile elbow (and catalytic triad) in every AlphaFold model.

Output: data/derived/structure_summary.csv (length, mean pLDDT per model) and
        data/derived/nucleophile_elbow.csv, one row per candidate residue (elbow-like phi/psi or a
match of the legacy notebook regex). `elbow` is the final call; `regex` is the old motif test;
`has_triad` says whether a nucleophile-His-Asp/Glu triad fits the distance cutoffs in config.py.
"""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.elbow import scan_structure, structure_summary


def main():
    rows, summ = [], []
    files = sorted(C.STRUCTURES.glob("*.pdb"))
    for f in files:
        rows += scan_structure(f)
        summ.append(structure_summary(f))
    df = pd.DataFrame(rows)
    C.DERIVED.mkdir(parents=True, exist_ok=True)
    df.to_csv(C.ELBOW_CSV, index=False)
    pd.DataFrame(summ).to_csv(C.STRUCT_SUMMARY, index=False)
    print(f"{len(files)} structures; {df.accession.nunique()} with candidates; "
          f"{df[df.elbow].accession.nunique()} with an elbow -> {C.ELBOW_CSV}")


if __name__ == "__main__":
    main()
