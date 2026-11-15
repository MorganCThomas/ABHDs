"""Catalytic-site calls used by study stages 3-4: elbow at any residue type, triad with chain-order constraint,
oxyanion-hole candidate. Thresholds that define the elbow live in config.py; `scan_elbow_table` keeps wide raw
features so thresholds can be varied afterwards (sensitivity analysis) without re-reading structures."""
import math

import numpy as np
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.NeighborSearch import NeighborSearch

from . import config as C
from .elbow import THREE_TO_ONE, _deg

NUC_ATOMS = {"SER": ["OG"], "CYS": ["SG"], "THR": ["OG1"], "ASP": ["OD1", "OD2"]}


def scan_elbow_table(path, accession=None, lo=None, hi=None):
    """One row per residue with phi in [0,120] and psi in [-180,-40] (wide net), any residue type, with the
    neighbouring-residue conformations needed to apply the elbow definition later."""
    accession = accession or path.stem
    structure = PDBParser(QUIET=True).get_structure(accession, str(path))
    rows = []
    for pp in PPBuilder().build_peptides(structure):
        res = list(pp)
        ang = [_deg(p) for p in pp.get_phi_psi_list()]
        for i in range(1, len(res) - 1):
            phi, psi = ang[i]
            if phi is None or psi is None or not (0 <= phi <= 120 and -180 <= psi <= -40):
                continue
            r = res[i]; num = r.get_id()[1]
            if (lo and num < lo) or (hi and num > hi):
                continue
            ctx = "".join(THREE_TO_ONE.get(x.get_resname(), "X") for x in res[max(0, i - 2):i + 3])
            rows.append(dict(accession=accession, resnum=num, resname=r.get_resname(),
                             aa=THREE_TO_ONE.get(r.get_resname(), "X"), context=ctx, phi=phi, psi=psi,
                             plddt=r["CA"].get_bfactor(),
                             prev_phi=ang[i - 1][0], prev_psi=ang[i - 1][1],
                             next_phi=ang[i + 1][0], next_psi=ang[i + 1][1]))
    return rows


def apply_elbow(df, phi=C.ELBOW_PHI, psi=C.ELBOW_PSI, plddt_min=C.PLDDT_MIN, need_ss=True):
    """Boolean Series: row satisfies the elbow definition under the given thresholds."""
    ok = df.phi.between(*phi) & df.psi.between(*psi) & (df.plddt >= plddt_min)
    if need_ss:
        beta = (df.prev_phi <= C.BETA_PHI_MAX) & ((df.prev_psi >= C.BETA_PSI_MIN) | (df.prev_psi <= C.BETA_PSI_WRAP))
        helix = df.next_phi.between(*C.HELIX_PHI) & df.next_psi.between(*C.HELIX_PSI)
        ok &= beta & helix
    return ok


def triad_and_oxyanion(path, nuc_resnum, accession=None):
    """For the residue `nuc_resnum`: best His/acid pair (distance cutoffs from config) preferring the canonical chain
    order nucleophile < acid < His, plus backbone-N atoms near the nucleophile (oxyanion-hole candidates)."""
    structure = PDBParser(QUIET=True).get_structure(accession or path.stem, str(path))
    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)
    nuc_res = next((r for r in structure.get_residues() if r.get_id()[1] == nuc_resnum), None)
    out = dict(his=None, acid=None, acid_type=None, d_nuc_his=None, d_his_acid=None, order_ok=None,
               n_order_ok_triads=0, oxy_dist=None, oxy_res=None)
    if nuc_res is None:
        return out
    names = NUC_ATOMS.get(nuc_res.get_resname(), [])
    cands = []
    for an in names:
        if an not in nuc_res:
            continue
        nuc = nuc_res[an]
        for a in ns.search(nuc.coord, C.D_NUC_HIS, level="A"):
            if a.get_name() not in ("NE2", "ND1"):
                continue
            his = a.get_parent()
            if his.get_resname() != "HIS" or his is nuc_res:
                continue
            nd = his["ND1" if a.get_name() == "NE2" else "NE2"]
            for o in ns.search(nd.coord, C.D_HIS_ACID, level="A"):
                if o.get_name() not in ("OD1", "OD2", "OE1", "OE2"):
                    continue
                acid = o.get_parent()
                if acid.get_resname() not in ("ASP", "GLU") or acid is his or acid is nuc_res:
                    continue
                d1 = float(np.linalg.norm(nuc.coord - a.coord)); d2 = float(np.linalg.norm(nd.coord - o.coord))
                order = nuc_res.get_id()[1] < acid.get_id()[1] < his.get_id()[1]
                cands.append((not order, d1 + d2, his.get_id()[1], acid.get_id()[1], acid.get_resname(), d1, d2, order))
    if cands:
        cands.sort()
        _, _, h, ac, at, d1, d2, order = cands[0]
        out.update(his=h, acid=ac, acid_type=at, d_nuc_his=round(d1, 2), d_his_acid=round(d2, 2), order_ok=order,
                   n_order_ok_triads=len({(c[2], c[3]) for c in cands if c[7]}))
    # oxyanion hole candidate: nearest non-local backbone N to the nucleophile heavy atom
    if names and names[0] in nuc_res:
        nuc = nuc_res[names[0]]
        best = None
        for a in ns.search(nuc.coord, 6.0, level="A"):
            if a.get_name() != "N":
                continue
            rr = a.get_parent()
            if abs(rr.get_id()[1] - nuc_resnum) <= 3:
                continue
            d = float(np.linalg.norm(nuc.coord - a.coord))
            if best is None or d < best[0]:
                best = (d, f"{THREE_TO_ONE.get(rr.get_resname(), 'X')}{rr.get_id()[1]}")
        if best:
            out.update(oxy_dist=round(best[0], 2), oxy_res=best[1])
    return out


def his_near_ca(path, resnum, cutoff=9.0, accession=None):
    """True if a His ring nitrogen lies within `cutoff` of the CA of `resnum` (used to rank elbows whose residue is not
    a nucleophile type, e.g. ABHD5 Asn153, against incidental phi>0 glycines)."""
    structure = PDBParser(QUIET=True).get_structure(accession or path.stem, str(path))
    ns = NeighborSearch(list(structure.get_atoms()))
    res = next((r for r in structure.get_residues() if r.get_id()[1] == resnum), None)
    if res is None or "CA" not in res:
        return False
    return any(a.get_name() in ("NE2", "ND1") and a.get_parent().get_resname() == "HIS" and a.get_parent() is not res
               for a in ns.search(res["CA"].coord, cutoff, level="A"))
