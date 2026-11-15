#!/usr/bin/env python3
"""
Follow-up B (docs/followup_design.md): oxyanion-hole donors and grid-cavity pocket descriptors for all genes.

Cavity: 1 A grid around the anchor (nucleophile O/S, else elbow CA); points with clearance >= CLEAR from every heavy atom and
buried (>= BUR of 30 rays blocked within 12 A); the 26-connected component nearest the anchor is the active-site cavity.
Outputs: data/derived/pocket_descriptors.csv, pocket_lining.csv, oxyanion_donors.csv, pocket_vs_length.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBParser
from scipy.ndimage import label
from scipy.spatial import cKDTree
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.elbow import THREE_TO_ONE

HALF = 14
BUR = 21
HYD = set("AVILMFWYC")
NUC_ATOM = {"SER": "OG", "CYS": "SG", "THR": "OG1", "ASP": "OD1"}
DIRS = None


def fib(n=30):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def cavity(atoms_xyz, anchor, clear):
    tree = cKDTree(atoms_xyz)
    g = np.arange(-HALF, HALF + 1)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    pts = np.c_[X.ravel(), Y.ravel(), Z.ravel()] + np.round(anchor)
    d, _ = tree.query(pts)
    free = d >= clear
    idx = np.where(free)[0]
    P = pts[idx]
    blocked = np.zeros(len(P), int)
    for u in DIRS:
        ray = P[:, None, :] + u[None, None, :] * np.arange(1, 13)[None, :, None]
        dd, _ = tree.query(ray.reshape(-1, 3))
        blocked += (dd.reshape(len(P), 12) < 1.9).any(axis=1)
    keep = np.zeros(len(pts), bool); keep[idx[blocked >= BUR]] = True
    vol = keep.reshape(X.shape)
    lab, n = label(vol, structure=np.ones((3, 3, 3)))
    if n == 0:
        return np.zeros((0, 3)), False
    coords = pts.reshape(*X.shape, 3)
    best, bestn = 0, 0
    for k in range(1, n + 1):
        c = coords[lab == k]
        m = int((np.linalg.norm(c - anchor, axis=1) < 7).sum())
        if m > bestn:
            best, bestn = k, m
    if best == 0:
        return np.zeros((0, 3)), False
    sel = lab == best
    touches = bool(sel[0].any() or sel[-1].any() or sel[:, 0].any() or sel[:, -1].any() or sel[:, :, 0].any() or sel[:, :, -1].any())
    return coords[sel], touches


def main():
    global DIRS
    DIRS = fib()
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv")
    qc = pd.read_csv(C.DERIVED / "structure_qc.csv").set_index("accession")
    parser = PDBParser(QUIET=True)
    desc, lin, oxy = [], [], []
    for r in cs.itertuples():
        if pd.isna(r.elbow_resnum):
            continue
        path = C.DERIVED / "domains" / f"{r.accession}.pdb"
        s = parser.get_structure(r.accession, str(path))
        res = {x.get_id()[1]: x for x in s.get_residues()}
        n = int(r.elbow_resnum)
        if n not in res:
            continue
        atoms = [a for a in s.get_atoms() if a.element != "H"]
        xyz = np.array([a.coord for a in atoms]); owner = np.array([a.get_parent().get_id()[1] for a in atoms])
        el = res[n]
        nucname = NUC_ATOM.get(el.get_resname())
        has_nuc = r.category and not str(r.category).startswith("no nucleophile") and nucname in el
        anchor = el[nucname].coord if has_nuc else el["CA"].coord
        row = dict(gene=r.gene, accession=r.accession, abhd_named=r.abhd_named, category=r.category, anchor="nuc" if has_nuc else "CA",
                   dom_len=int(qc.loc[r.accession, "dom_len"]))
        for clear in (3.0, 3.2, 3.6):
            pts, touch = cavity(xyz, anchor, clear)
            tag = f"{clear:.1f}"
            row[f"vol_{tag}"] = len(pts); row[f"touches_box_{tag}"] = touch
            if clear == 3.2:
                if len(pts):
                    t = cKDTree(pts)
                    near = np.array([t.query(a.coord)[0] <= 1.0 + 3.4 for a in atoms])      # heavy atom within ~4.4 A of a cavity point
                    ids = sorted(set(owner[near]))
                    for i in ids:
                        rr = res[i]
                        lin.append(dict(gene=r.gene, accession=r.accession, resnum=i, aa=THREE_TO_ONE.get(rr.get_resname(), "X")))
                    seq = [THREE_TO_ONE.get(res[i].get_resname(), "X") for i in ids]
                    row.update(n_lining=len(ids), frac_hydrophobic=round(np.mean([a in HYD for a in seq]), 3),
                               n_aromatic=sum(a in "FWY" for a in seq), n_charged=sum(a in "DEKR" for a in seq),
                               n_polar=sum(a in "STNQH" for a in seq), max_extent=round(float(np.linalg.norm(pts - anchor, axis=1).max()), 1),
                               centroid_offset=round(float(np.linalg.norm(pts.mean(0) - anchor)), 1))
                    np.save(C.DERIVED / "domains" / f"{r.accession}.cavity.npy", pts)
        desc.append(row)
        # oxyanion donors: backbone N within 3.6 A of the nucleophile heavy atom (nucleophile genes only)
        if has_nuc:
            near = []
            for i, rr in res.items():
                if i == n or "N" not in rr:
                    continue
                d = float(np.linalg.norm(rr["N"].coord - anchor))
                if d <= 5.0:
                    ctx = "".join(THREE_TO_ONE.get(res[j].get_resname(), "-") if j in res else "-" for j in range(i - 3, i + 4))
                    near.append((d, i, THREE_TO_ONE.get(rr.get_resname(), "X"), i - n, ctx))
            near.sort()
            for d, i, a, off, ctx in near[:4]:
                oxy.append(dict(gene=r.gene, accession=r.accession, nuc=n, donor_resnum=i, donor_aa=a, offset_from_nuc=off,
                                dist=round(d, 2), within_3p6=d <= 3.6, context=ctx))
    D = pd.DataFrame(desc); D.to_csv(C.DERIVED / "pocket_descriptors.csv", index=False)
    pd.DataFrame(lin).to_csv(C.DERIVED / "pocket_lining.csv", index=False)
    O = pd.DataFrame(oxy); O.to_csv(C.DERIVED / "oxyanion_donors.csv", index=False)
    rho = spearmanr(D["vol_3.2"], D.dom_len)
    print(f"genes: {len(D)}; median volume {D['vol_3.2'].median():.0f}; vol vs domain length rho={rho[0]:.2f} p={rho[1]:.1e}")
    sl = np.polyfit(np.log(D.dom_len), np.log1p(D["vol_3.2"]), 1)
    D["vol_resid_z"] = ((np.log1p(D["vol_3.2"]) - np.polyval(sl, np.log(D.dom_len))))
    D["vol_resid_z"] = (D.vol_resid_z - D.vol_resid_z.mean()) / D.vol_resid_z.std()
    D.to_csv(C.DERIVED / "pocket_descriptors.csv", index=False)
    print("volume sensitivity (Spearman 3.0 vs 3.6):", round(spearmanr(D["vol_3.0"], D["vol_3.6"])[0], 2),
          "| 3.0 vs 3.2:", round(spearmanr(D["vol_3.0"], D["vol_3.2"])[0], 2))
    print("boxes touched at 3.2:", int(D["touches_box_3.2"].sum()))
    D[["gene", "vol_3.2"]].assign(vol_dom=D.dom_len).to_csv(C.DERIVED / "pocket_vs_length.csv", index=False)
    n_in = O[O.within_3p6].groupby("gene").size()
    print("genes with >=1 N within 3.6 A:", len(n_in), "of", O.gene.nunique(), "; with >=2:", int((n_in >= 2).sum()))


if __name__ == "__main__":
    main()
