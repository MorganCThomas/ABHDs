#!/usr/bin/env python3
"""
Study stage 2 - structure quality control and trimming to the hydrolase domain.

Domain definition: for proteins in a CATH FunFam, the CATH domain range(s) recorded in the FunFam alignment headers
(`ACC/start-end[_start-end]`; span of all ranges over all FunFams that contain the protein). For proteins not in any FunFam the longest run of residues with
pLDDT >= PLDDT_MIN (gaps of <= 15 low-confidence residues bridged) is used, and flagged `domain_source=plddt`.
Outputs:
  data/derived/structure_qc.csv         one row per universe gene
  data/derived/domains/<acc>.pdb        model trimmed to the domain (original residue numbering kept)
  data/derived/attrition_stage2.csv
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.PDB import PDBIO, PDBParser, Select

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C

DOMAINS = C.DERIVED / "domains"
MIN_DOMAIN_LEN = 150
GAP_BRIDGE = 15


def funfam_ranges():
    rng = {}
    for f in sorted((C.FUNFAM_DIR / "FF_alns").glob("*.faa")):
        for line in open(f):
            if line.startswith(">"):
                m = re.match(r">([A-Z0-9]+)/(\S+)", line.strip())
                if m:
                    segs = [tuple(map(int, s.split("-"))) for s in m.group(2).split("_")]
                    rng.setdefault(m.group(1), []).append((f.stem, segs))
    return rng


def plddt_domain(res_ids, plddt):
    ok = np.where(np.asarray(plddt) >= C.PLDDT_MIN)[0]
    if not len(ok):
        return None
    runs, start = [], ok[0]
    for a, b in zip(ok, ok[1:]):
        if b - a - 1 > GAP_BRIDGE:
            runs.append((start, a)); start = b
    runs.append((start, ok[-1]))
    s, e = max(runs, key=lambda r: r[1] - r[0])
    return res_ids[s], res_ids[e]


class Trim(Select):
    def __init__(self, lo, hi): self.lo, self.hi = lo, hi
    def accept_residue(self, r): return self.lo <= r.get_id()[1] <= self.hi


def feature_span(text, kind):
    return [tuple(map(int, m)) for m in re.findall(rf"{kind} (\d+)\.\.(\d+)", str(text))]


def main():
    DOMAINS.mkdir(parents=True, exist_ok=True)
    uni = pd.read_csv(C.DERIVED / "universe.csv")
    ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
    ff = funfam_ranges()
    parser, io, rows = PDBParser(QUIET=True), PDBIO(), []
    for r in uni.itertuples():
        pdb = C.STRUCTURES / f"{r.accession}.pdb"
        row = dict(gene=r.gene, accession=r.accession, abhd_named=r.abhd_named, source=r.source)
        if not pdb.exists():
            rows.append({**row, "flag": "no_model"}); continue
        s = parser.get_structure(r.accession, pdb)
        res = [x for x in s.get_residues() if "CA" in x]
        ids = [x.get_id()[1] for x in res]
        pl = np.array([x["CA"].get_bfactor() for x in res])
        if r.accession in ff:
            # a protein can be split over several FunFams (e.g. LCAT: 25-93, 98-244_329-354, 355-423): take the union span
            allsegs = [sg for _, segs in ff[r.accession] for sg in segs]
            ff_id = ";".join(sorted({fid for fid, _ in ff[r.accession]}))
            lo, hi, src = min(a for a, b in allsegs), max(b for a, b in allsegs), "cath_funfam"
        else:
            lo, hi = plddt_domain(ids, pl); src, ff_id = "plddt", ""
        m = (np.array(ids) >= lo) & (np.array(ids) <= hi)
        sig = feature_span(ann.loc[r.accession, "Signal peptide"], "SIGNAL") if r.accession in ann.index else []
        tm = feature_span(ann.loc[r.accession, "Transmembrane"], "TRANSMEM") if r.accession in ann.index else []
        io.set_structure(s); io.save(str(DOMAINS / f"{r.accession}.pdb"), Trim(lo, hi))
        dom_mean = float(pl[m].mean())
        flags = []
        if dom_mean < C.PLDDT_MIN: flags.append("low_plddt")
        if m.sum() < MIN_DOMAIN_LEN: flags.append("short_domain")
        if src == "plddt": flags.append("domain_from_plddt")
        rows.append({**row, "n_res": len(res), "mean_plddt": round(float(pl.mean()), 1), "domain_source": src,
                     "funfam": ff_id, "dom_start": lo, "dom_end": hi, "dom_len": int(m.sum()),
                     "dom_mean_plddt": round(dom_mean, 1), "dom_frac_ge70": round(float((pl[m] >= 70).mean()), 2),
                     "signal_peptide": bool(sig), "n_tm": len(tm), "flag": ";".join(flags) or "ok"})
    qc = pd.DataFrame(rows)
    qc.to_csv(C.DERIVED / "structure_qc.csv", index=False)
    att = qc.flag.str.split(";").explode().value_counts().rename_axis("flag").reset_index(name="n_genes")
    att.to_csv(C.DERIVED / "attrition_stage2.csv", index=False)
    print(att.to_string(index=False))
    print(qc.dom_len.describe().round(0).to_dict()); print(qc.dom_mean_plddt.describe().round(1).to_dict())


if __name__ == "__main__":
    main()
