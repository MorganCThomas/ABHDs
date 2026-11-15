#!/usr/bin/env python3
"""
Study stage 5 - local comparison of the catalytic site across the universe.

(a) Superposition RMSD (Kabsch) of the elbow backbone (N, CA, C of residues i-2..i+2; all genes with an elbow call) and
    of a triad "skeleton" (CA + functional-group atom of nucleophile, His, acid; genes with a complete triad).
(b) SOAP descriptors (r_cut/n_max/l_max from config) on the domain model, centred on the CA of the elbow residue, so genes
    with and without a nucleophile are comparable. Null baseline: SOAP at CA of random non-elbow serines.
(c) Crude pocket-environment proxy: composition of residues within 10 A of the nucleophile.
Outputs (data/derived/): site_rmsd_elbow.csv, site_rmsd_triad.csv, soap_site.npz, soap_site_cosine.csv,
soap_null_summary.csv, site_env.csv; figures/stage5_*.html
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from ase import Atoms
from Bio.PDB import PDBParser
from dscribe.descriptors import SOAP

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

DOMAINS = C.DERIVED / "domains"
HYD, CHG, POL = set("AVLIMFWC"), set("DEKR"), set("STNQHY")
FUNC = {"SER": ["OG"], "CYS": ["SG"], "THR": ["OG1"], "ASP": ["OD1", "OD2"], "GLU": ["OE1", "OE2"],
        "HIS": ["CG", "ND1", "CD2", "CE1", "NE2"]}
ONE = {"ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H",
       "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
       "TYR": "Y", "VAL": "V"}


def kabsch_rmsd(P, Q):
    P, Q = P - P.mean(0), Q - Q.mean(0)
    U, S, Vt = np.linalg.svd(P.T @ Q)
    d = np.sign(np.linalg.det(U @ Vt))
    D = np.diag([1, 1, d])
    R = U @ D @ Vt
    return float(np.sqrt(((P @ R - Q) ** 2).sum() / len(P)))


def centroid(res, names):
    xs = [res[n].coord for n in names if n in res]
    return np.mean(xs, axis=0) if xs else None


def rmsd_matrix(sets):
    keys = list(sets)
    M = pd.DataFrame(np.nan, index=keys, columns=keys)
    for i, a in enumerate(keys):
        for b in keys[i:]:
            v = kabsch_rmsd(sets[a], sets[b]); M.loc[a, b] = M.loc[b, a] = v
    return M


def to_ase(structure):
    atoms = [a for a in structure.get_atoms()]
    sym = [(a.element.strip() or a.get_name()[0]).capitalize() for a in atoms]
    return atoms, Atoms(symbols=sym, positions=np.array([a.coord for a in atoms]))


def main():
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv")
    cs = cs[cs.elbow_resnum.notna()].copy()
    parser = PDBParser(QUIET=True)
    elbow_pts, triad_pts, env, loaded = {}, {}, [], {}
    for r in cs.itertuples():
        s = parser.get_structure(r.accession, str(DOMAINS / f"{r.accession}.pdb"))
        res = {x.get_id()[1]: x for x in s.get_residues()}
        n = int(r.elbow_resnum)
        try:
            elbow_pts[r.accession] = np.array([res[k][a].coord for k in range(n - 2, n + 3) for a in ("N", "CA", "C")])
        except KeyError:
            continue
        loaded[r.accession] = (s, res, n)
        nuc = res[n]
        # triad skeleton
        if pd.notna(r.his) and pd.notna(r.acid) and nuc.get_resname() in FUNC and int(r.his) in res and int(r.acid) in res:
            h, ac = res[int(r.his)], res[int(r.acid)]
            pts = [nuc["CA"].coord, centroid(nuc, FUNC[nuc.get_resname()]), h["CA"].coord, centroid(h, FUNC["HIS"]),
                   ac["CA"].coord, centroid(ac, FUNC[ac.get_resname()])]
            if all(p is not None for p in pts):
                triad_pts[r.accession] = np.array(pts)
        # pocket proxy: residues within 10 A of the nucleophile functional atom (or CA)
        c = centroid(nuc, FUNC.get(nuc.get_resname(), [])) if nuc.get_resname() in FUNC else None
        c = nuc["CA"].coord if c is None else c
        near = [x for x in res.values() if any(np.linalg.norm(a.coord - c) <= 10.0 for a in x) and x is not nuc]
        aas = [ONE.get(x.get_resname(), "X") for x in near]
        env.append(dict(accession=r.accession, gene=r.gene, n_res_10A=len(near),
                        frac_hydrophobic=np.mean([a in HYD for a in aas]), frac_charged=np.mean([a in CHG for a in aas]),
                        frac_polar=np.mean([a in POL for a in aas]),
                        n_atoms_10A=sum(1 for x in res.values() for a in x if np.linalg.norm(a.coord - c) <= 10.0)))
    gene = dict(zip(cs.accession, cs.gene))
    for name, pts in (("elbow", elbow_pts), ("triad", triad_pts)):
        M = rmsd_matrix(pts); M.to_csv(C.DERIVED / f"site_rmsd_{name}.csv")
        print(name, M.shape, "median off-diagonal RMSD %.2f A" % np.nanmedian(M.values[np.triu_indices(len(M), 1)]))
    pd.DataFrame(env).to_csv(C.DERIVED / "site_env.csv", index=False)

    # ---- SOAP at CA of the elbow residue + null baseline ----
    species, asedata = set(), {}
    for acc, (s, res, n) in loaded.items():
        atoms, ase = to_ase(s)
        asedata[acc] = (atoms, ase); species.update(ase.get_chemical_symbols())
    soap = SOAP(species=sorted(species), periodic=False, r_cut=C.SOAP_RCUT, n_max=C.SOAP_NMAX, l_max=C.SOAP_LMAX)
    rng = np.random.default_rng(0)
    X, accs, N, Nlab = [], [], [], []
    for acc, (atoms, ase) in asedata.items():
        n = loaded[acc][2]
        idx = next(i for i, a in enumerate(atoms) if a.get_parent().get_id()[1] == n and a.get_name() == "CA")
        X.append(soap.create(ase, centers=[idx])[0]); accs.append(acc)
        sers = [i for i, a in enumerate(atoms) if a.get_name() == "CA" and a.get_parent().get_resname() == "SER"
                and abs(a.get_parent().get_id()[1] - n) > 5]
        for i in rng.choice(sers, size=min(3, len(sers)), replace=False):
            N.append(soap.create(ase, centers=[int(i)])[0]); Nlab.append(acc)
    X, N = np.array(X), np.array(N)
    np.savez(C.DERIVED / "soap_site.npz", accessions=np.array(accs), X=X)

    def cos(A, B):
        A = A / np.linalg.norm(A, axis=1, keepdims=True); B = B / np.linalg.norm(B, axis=1, keepdims=True)
        return 1 - A @ B.T
    D = cos(X, X); pd.DataFrame(D, index=accs, columns=accs).to_csv(C.DERIVED / "soap_site_cosine.csv")
    iu = np.triu_indices(len(X), 1)
    en = cos(X, N); mask = np.array([[a != b for b in Nlab] for a in accs])
    nn = cos(N, N); mk = np.array([[a != b for b in Nlab] for a in Nlab]); nu = np.triu(np.ones_like(nn, bool), 1) & mk
    summ = pd.DataFrame([
        dict(comparison="elbow site vs elbow site", n=len(iu[0]), median=np.median(D[iu]), p5=np.percentile(D[iu], 5), p95=np.percentile(D[iu], 95)),
        dict(comparison="elbow site vs random Ser CA", n=int(mask.sum()), median=np.median(en[mask]), p5=np.percentile(en[mask], 5), p95=np.percentile(en[mask], 95)),
        dict(comparison="random Ser CA vs random Ser CA", n=int(nu.sum()), median=np.median(nn[nu]), p5=np.percentile(nn[nu], 5), p95=np.percentile(nn[nu], 95))])
    summ.round(4).to_csv(C.DERIVED / "soap_null_summary.csv", index=False)
    print(summ.round(4).to_string(index=False))

    # figure: distance distributions
    fig = go.Figure()
    for lab, v in (("elbow vs elbow", D[iu]), ("elbow vs random Ser", en[mask]), ("random vs random", nn[nu])):
        fig.add_histogram(x=v, name=lab, opacity=0.6, histnorm="probability", nbinsx=60)
    fig.update_layout(barmode="overlay", template="simple_white", width=800, height=420,
                      title="Stage 5: SOAP cosine distance, elbow sites vs random serine sites", xaxis_title="cosine distance")
    fig.write_html(C.FIGURES / "stage5_soap_null.html")


if __name__ == "__main__":
    main()
