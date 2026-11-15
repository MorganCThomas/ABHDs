"""Shared structure helpers: domain loading, TM-align residue mapping, and a grid-cavity calculation."""
import numpy as np
from Bio.PDB import PDBParser
from scipy.ndimage import label
from scipy.spatial import cKDTree
from tmtools import tm_align

from . import config as C
from .elbow import THREE_TO_ONE


def load_pdb(path):
    """(residue numbers, CA coordinates, one-letter sequence) of a PDB file."""
    s = PDBParser(QUIET=True).get_structure("x", str(path))
    ids, ca, seq = [], [], []
    for r in s.get_residues():
        if "CA" in r:
            ids.append(r.get_id()[1]); ca.append(r["CA"].coord); seq.append(THREE_TO_ONE.get(r.get_resname(), "X"))
    return np.array(ids), np.array(ca), "".join(seq)


def load(acc):
    """Trimmed hydrolase domain (stage 2)."""
    return load_pdb(C.DERIVED / "domains" / f"{acc}.pdb")


def align(a, b):
    """TM-align of two trimmed domains. Returns ((ids,ca,seq)_a, (ids,ca,seq)_b, a_superposed_on_b, aligned index pairs, result)."""
    ia, ca, sa = load(a); ib, cb, sb = load(b)
    r = tm_align(ca, cb, sa, sb)
    X = ca @ r.u.T + r.t
    pairs, i, j = [], 0, 0
    for x, y in zip(r.seqxA, r.seqyA):
        if x != "-" and y != "-":
            pairs.append((i, j))
        i += x != "-"; j += y != "-"
    return (ia, ca, sa), (ib, cb, sb), X, pairs, r


def fib_dirs(n=30):
    i = np.arange(n) + 0.5
    phi = np.arccos(1 - 2 * i / n); th = np.pi * (1 + 5 ** 0.5) * i
    return np.c_[np.cos(th) * np.sin(phi), np.sin(th) * np.sin(phi), np.cos(phi)]


def random_rotation(rng):
    q = rng.normal(size=4); q /= np.linalg.norm(q)
    a, b, c, d = q
    return np.array([[a*a+b*b-c*c-d*d, 2*(b*c-a*d), 2*(b*d+a*c)],
                     [2*(b*c+a*d), a*a-b*b+c*c-d*d, 2*(c*d-a*b)],
                     [2*(b*d-a*c), 2*(c*d+a*b), a*a-b*b-c*c+d*d]])


def cavity_variants(atoms_xyz, anchor, clearances=(3.0, 3.2, 3.6), half=18, bur=21, rot=None, shift=None, n_dirs=30):
    """Active-site cavity at several probe clearances from one grid.

    A 1 A grid (optionally rotated by `rot` about the anchor and shifted by `shift`) spans +-half A around the anchor. Points
    with clearance >= the smallest clearance from every atom and buried (>= `bur` of `n_dirs` rays blocked within 12 A) are
    kept; for each clearance the 26-connected component closest to the anchor (most points within 7 A) is the cavity.
    Returns {clearance: (points in the original frame, touches_box)}.
    """
    R = np.eye(3) if rot is None else rot
    shift = np.zeros(3) if shift is None else shift
    X = (atoms_xyz - anchor) @ R.T                       # work in a frame centred on the anchor, rotated
    tree = cKDTree(X)
    g = np.arange(-half, half + 1)
    G = np.stack(np.meshgrid(g, g, g, indexing="ij"), -1).astype(float) + shift
    pts = G.reshape(-1, 3)
    d, _ = tree.query(pts)
    cmin = min(clearances)
    idx = np.where(d >= cmin)[0]
    P = pts[idx]; dirs = fib_dirs(n_dirs)
    blocked = np.zeros(len(P), int)
    steps = np.arange(1, 13)
    for u in dirs:
        ray = P[:, None, :] + u[None, None, :] * steps[None, :, None]
        dd, _ = tree.query(ray.reshape(-1, 3))
        blocked += (dd.reshape(len(P), 12) < 1.9).any(axis=1)
    buried = np.zeros(len(pts), bool); buried[idx[blocked >= bur]] = True
    out = {}
    for c in clearances:
        keep = buried & (d >= c)
        lab, n = label(keep.reshape(G.shape[:3]), structure=np.ones((3, 3, 3)))
        if n == 0:
            out[c] = (np.zeros((0, 3)), False); continue
        flat = lab.ravel()
        best, bestn = 0, 0
        near = np.linalg.norm(pts, axis=1) < 7
        counts = np.bincount(flat[near & (flat > 0)], minlength=n + 1)
        best = int(counts.argmax()); bestn = int(counts[best]) if best > 0 else 0
        if bestn == 0:
            out[c] = (np.zeros((0, 3)), False); continue
        sel = lab == best
        touches = bool(sel[0].any() or sel[-1].any() or sel[:, 0].any() or sel[:, -1].any() or sel[:, :, 0].any() or sel[:, :, -1].any())
        out[c] = (pts[sel.ravel()] @ R + anchor, touches)
    return out
