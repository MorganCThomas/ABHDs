"""Nucleophile-elbow detection from backbone dihedrals, plus catalytic-triad geometry."""
import math
import re

import numpy as np
from Bio.PDB import PDBParser, PPBuilder
from Bio.PDB.NeighborSearch import NeighborSearch

from . import config as C

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
_LEGACY = re.compile(C.LEGACY_REGEX)


def _one(res):
    return THREE_TO_ONE.get(res.get_resname(), "X")


def _deg(pair):
    return tuple(None if x is None else math.degrees(x) for x in pair)


def _between(x, lo_hi):
    return x is not None and lo_hi[0] <= x <= lo_hi[1]


def is_elbow_core(phi, psi):
    return _between(phi, C.ELBOW_PHI) and _between(psi, C.ELBOW_PSI)


def is_beta(phi, psi):
    return phi is not None and psi is not None and phi <= C.BETA_PHI_MAX and (psi >= C.BETA_PSI_MIN or psi <= C.BETA_PSI_WRAP)


def is_helix(phi, psi):
    return _between(phi, C.HELIX_PHI) and _between(psi, C.HELIX_PSI)


def find_triad(nuc_res, ns):
    """Best (smallest summed distance) nucleophile-His-acid triad around one nucleophile residue, or None."""
    best = None
    for atom_name in C.NUC_RESIDUES[nuc_res.get_resname()]:
        if atom_name not in nuc_res:
            continue
        nuc = nuc_res[atom_name]
        for n_prox in ns.search(nuc.coord, C.D_NUC_HIS, level="A"):
            if n_prox.get_name() not in ("NE2", "ND1"):
                continue
            his = n_prox.get_parent()
            if his.get_resname() != "HIS" or his is nuc_res:
                continue
            n_dist = his["ND1" if n_prox.get_name() == "NE2" else "NE2"]
            for o in ns.search(n_dist.coord, C.D_HIS_ACID, level="A"):
                if o.get_name() not in ("OD1", "OD2", "OE1", "OE2"):
                    continue
                acid = o.get_parent()
                if acid.get_resname() not in ("ASP", "GLU") or acid is his or acid is nuc_res:
                    continue
                d1 = float(np.linalg.norm(nuc.coord - n_prox.coord))
                d2 = float(np.linalg.norm(n_dist.coord - o.coord))
                if best is None or d1 + d2 < best["d_nuc_his"] + best["d_his_acid"]:
                    best = dict(nuc_atom=atom_name, his_resnum=his.get_id()[1], acid_resnum=acid.get_id()[1],
                                acid_resname=acid.get_resname(), d_nuc_his=round(d1, 2), d_his_acid=round(d2, 2))
    return best


def scan_structure(path, accession=None):
    """Return one dict per candidate residue (elbow-like dihedrals OR legacy-regex match)."""
    accession = accession or path.stem
    structure = PDBParser(QUIET=True).get_structure(accession, str(path))
    ns = NeighborSearch(list(structure.get_atoms()))
    rows = []
    for pp in PPBuilder().build_peptides(structure):
        res = list(pp)
        angles = [_deg(p) for p in pp.get_phi_psi_list()]
        for i in range(2, len(res) - 2):
            r = res[i]
            if r.get_resname() not in C.NUC_RESIDUES:
                continue
            ctx = "".join(_one(x) for x in res[i - 2:i + 3])
            regex = bool(_LEGACY.match(ctx))
            phi, psi = angles[i]
            core = is_elbow_core(phi, psi)
            if not (core or regex):
                continue
            plddt = r["CA"].get_bfactor()
            beta_before = is_beta(*angles[i - 1])
            helix_after = is_helix(*angles[i + 1])
            elbow = core and beta_before and helix_after and plddt >= C.PLDDT_MIN
            triad = find_triad(r, ns) or {}
            rows.append(dict(accession=accession, resnum=r.get_id()[1], aa=_one(r), context=ctx,
                             phi=round(phi, 1) if phi is not None else None,
                             psi=round(psi, 1) if psi is not None else None,
                             plddt=round(plddt, 1), regex=regex, core=core,
                             beta_before=beta_before, helix_after=helix_after, elbow=elbow,
                             has_triad=bool(triad), **triad))
    return rows


def structure_summary(path):
    """Length and mean pLDDT of one AlphaFold model (CA atoms)."""
    structure = PDBParser(QUIET=True).get_structure(path.stem, str(path))
    ca = [a for a in structure.get_atoms() if a.get_name() == "CA"]
    return dict(accession=path.stem, n_residues=len(ca), mean_plddt=round(float(np.mean([a.get_bfactor() for a in ca])), 1))
