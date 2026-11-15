"""SOAP descriptors centred on the elbow nucleophile atom of each structure.

Input : data/derived/nucleophile_elbow.csv (rows with elbow == True)
Output: data/derived/soap.npz (accessions, X) and cosine distance tables for ProtBert and SOAP.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from ase import Atoms
from Bio.PDB import PDBParser
from dscribe.descriptors import SOAP

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.data import cosine_matrix, load_protbert


def to_ase(structure):
    atoms = list(structure.get_atoms())
    sym = [(a.element.strip() or a.get_name()[0]).capitalize() for a in atoms]
    return atoms, Atoms(symbols=sym, positions=np.array([a.coord for a in atoms]))


def main():
    el = pd.read_csv(C.ELBOW_CSV)
    el = el[el.elbow].drop_duplicates("accession")
    parser = PDBParser(QUIET=True)
    loaded = {}
    species = set()
    for r in el.itertuples():
        s = parser.get_structure(r.accession, str(C.STRUCTURES / f"{r.accession}.pdb"))
        atoms, ase = to_ase(s)
        loaded[r.accession] = (atoms, ase)
        species.update(ase.get_chemical_symbols())
    soap = SOAP(species=sorted(species), periodic=False, r_cut=C.SOAP_RCUT, n_max=C.SOAP_NMAX, l_max=C.SOAP_LMAX)
    accs, feats = [], []
    for r in el.itertuples():
        atoms, ase = loaded[r.accession]
        atom_name = r.nuc_atom if isinstance(r.nuc_atom, str) else {"S": "OG", "D": "OD1"}[r.aa]
        idx = next(i for i, a in enumerate(atoms)
                   if a.get_parent().get_id()[1] == r.resnum and a.get_name() == atom_name)
        feats.append(soap.create(ase, centers=[idx])[0])
        accs.append(r.accession)
    X = np.array(feats)
    np.savez(C.SOAP_NPZ, accessions=np.array(accs), X=X)
    cosine_matrix(X, accs).to_csv(C.SOAP_COS)
    pb_accs, pb_X = load_protbert()
    cosine_matrix(pb_X, pb_accs).to_csv(C.PROTBERT_COS)
    print(f"SOAP {X.shape} for {len(accs)} structures; wrote cosine tables")


if __name__ == "__main__":
    main()
