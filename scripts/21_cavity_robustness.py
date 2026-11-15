#!/usr/bin/env python3
"""
Cavity redo (reports/12): active-site cavity volumes under model, anchor, grid rotation, probe clearance and box-size variants.

Variants per gene: model {trimmed domain, full-length AlphaFold model with pLDDT < 50 atoms removed} x anchor {native (nucleophile
atom, else elbow CA), ca_trunc (elbow CA with the elbow residue's side chain removed for every gene: like-for-like)} x 3 grid
orientations (identity + 2 random rotations with sub-voxel shifts) x clearance {3.0, 3.2, 3.6 A}; grid half-width 18 A.
Outputs: data/derived/cavity_variants.csv, cavity_summary.csv, cavity_contrasts.csv
"""
import sys
from itertools import product
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.structops import cavity_variants, random_rotation

CLEAR = (3.0, 3.2, 3.6)
NUC_ATOM = {"SER": "OG", "CYS": "SG", "THR": "OG1", "ASP": "OD1"}
BACKBONE = {"N", "CA", "C", "O"}
CONTRASTS = [("ABHD15", "ABHD3"), ("ABHD15", "ABHD1"), ("ABHD15", "ABHD2"), ("ABHD5", "ABHD4"), ("ABHD2", "MGLL"), ("ABHD2", "ABHD6"),
             ("ABHD6", "MGLL"), ("ABHD12", "MGLL"), ("ABHD10", "MGLL"), ("ABHD16A", "MGLL"), ("ABHD3", "ABHD2")]


def job(args):
    gene, acc, elbow, native_ok, model = args
    path = (C.DERIVED / "domains" / f"{acc}.pdb") if model == "domain" else (C.STRUCTURES / f"{acc}.pdb")
    s = PDBParser(QUIET=True).get_structure(acc, str(path))
    res = {r.get_id()[1]: r for r in s.get_residues()}
    if elbow not in res:
        return []
    atoms = [a for a in s.get_atoms() if a.element != "H" and (model == "domain" or a.get_bfactor() >= 50)]
    rows = []
    for anchor_mode in ("native", "ca_trunc"):
        el = res[elbow]
        if anchor_mode == "native":
            nm = NUC_ATOM.get(el.get_resname())
            anchor = el[nm].coord if (native_ok and nm in el) else el["CA"].coord
            use = atoms
        else:
            anchor = el["CA"].coord
            use = [a for a in atoms if not (a.get_parent() is el and a.get_name() not in BACKBONE)]
        xyz = np.array([a.coord for a in use])
        rng = np.random.default_rng(abs(hash((acc, model, anchor_mode))) % (2 ** 32))
        for k in range(3):
            rot = None if k == 0 else random_rotation(rng)
            shift = None if k == 0 else rng.uniform(-0.5, 0.5, 3)
            out = cavity_variants(xyz, anchor, CLEAR, half=18, rot=rot, shift=shift)
            for c, (pts, touch) in out.items():
                rows.append(dict(gene=gene, accession=acc, model=model, anchor=anchor_mode, orientation=k, clearance=c,
                                 volume=len(pts), touches_box=touch))
    return rows


def main():
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv")
    cs = cs[cs.elbow_resnum.notna()]
    jobs = [(r.gene, r.accession, int(r.elbow_resnum), not str(r.category).startswith("no nucleophile"), m)
            for r in cs.itertuples() for m in ("domain", "full")]
    with Pool() as p:
        rows = [x for part in p.imap_unordered(job, jobs, chunksize=2) for x in part]
    V = pd.DataFrame(rows); V.to_csv(C.DERIVED / "cavity_variants.csv", index=False)
    print("variant rows:", len(V))

    L = V[V.anchor == "ca_trunc"]
    S = (L.groupby(["gene", "clearance"]).volume.agg(median="median", vmin="min", vmax="max").unstack("clearance"))
    S.columns = [f"{a}_{b}" for a, b in S.columns]
    S["frac_box_truncated"] = L.groupby("gene").touches_box.mean()
    S.reset_index().to_csv(C.DERIVED / "cavity_summary.csv", index=False)
    print("fraction of ca_trunc variants touching the box: %.2f" % L.touches_box.mean())

    # agreement between model variants and anchors (genome-wide)
    piv = V.pivot_table(index="gene", columns=["model", "anchor", "clearance"], values="volume", aggfunc="median")
    ref = piv[("domain", "ca_trunc", 3.2)]
    for col in piv.columns:
        print("Spearman vs domain/ca_trunc/3.2:", col, round(spearmanr(ref, piv[col], nan_policy="omit")[0], 2))

    # contrasts: in how many (model, orientation, clearance) combos is the first gene's cavity larger?
    rows = []
    for anchor_mode in ("ca_trunc", "native"):
        for a, b in CONTRASTS:
            sub = V[V.anchor == anchor_mode].pivot_table(index=["model", "orientation", "clearance"], columns="gene", values="volume")
            if a not in sub or b not in sub:
                continue
            d = sub[a] - sub[b]
            rows.append(dict(anchor=anchor_mode, first=a, second=b, n_combos=len(d), first_larger=int((d > 0).sum()),
                             second_larger=int((d < 0).sum()), median_diff=float(d.median()),
                             consistent=bool((d > 0).all() or (d < 0).all())))
    Cn = pd.DataFrame(rows); Cn.to_csv(C.DERIVED / "cavity_contrasts.csv", index=False)
    print(Cn.to_string(index=False))
    pd.set_option("display.width", 220)
    panel = ["MGLL", "ABHD6", "ABHD10", "ABHD12", "ABHD16A", "ABHD2", "ABHD1", "ABHD3", "ABHD4", "ABHD5", "ABHD15", "ABHD13", "ABHD14B", "ABHD17A", "LYPLA1"]
    print(S.reindex(panel)[["median_3.2", "vmin_3.2", "vmax_3.2", "median_3.6", "vmin_3.6", "vmax_3.6", "frac_box_truncated"]].to_string())


if __name__ == "__main__":
    main()
