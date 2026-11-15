#!/usr/bin/env python3
"""Smoke test for the protein tools (needs data/derived from the study). Run: uv run python tests/smoke.py"""
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd.atlas import Atlas
from abhd.compare import compare, to_markdown
from abhd.protein_report import build_report

A = Atlas()
i = A.info("abhd2")                                       # case-insensitive lookup
assert i["accession"] == "P08910" and i["category"] == "Ser-His-Asp" and int(i["elbow_resnum"]) == 207
n = A.neighbours("ABHD2", "structure", k=5)
assert len(n) == 5 and n.gene.iloc[0] == "ABHD1" and (n.distance.diff().dropna() >= 0).all()
assert A.predict_ec("ABHD15")["gated"] is True             # no-nucleophile proteins get no EC transfer
assert abs(A.identity(A.acc("ABHD2"), A.acc("ABHD3")) - A.identity(A.acc("ABHD3"), A.acc("ABHD2"))) < 1e-12
try:
    A.acc("NOTAGENE"); raise SystemExit("expected KeyError")
except KeyError:
    pass
c = compare(A, "ABHD2", "MGLL")
assert c["catalytic"].b_aligned_residue.tolist() == ["S122", "H269", "D239"]   # catalytic residues map onto each other
assert "ABHD2 vs MGLL" in to_markdown(c)
with tempfile.TemporaryDirectory() as d:
    r = build_report(A, "SPG21", d, structure=False)
    assert Path(r["path"]).stat().st_size > 50_000
print("smoke test passed")
