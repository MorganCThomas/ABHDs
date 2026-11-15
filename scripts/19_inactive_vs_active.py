#!/usr/bin/env python3
"""
Follow-up A (docs/followup_design.md): what separates inactive ABHD15 / ABHD5 from their active relatives?

A1 alternative nucleophile at the elbow window; A2 residues at catalytic positions via TM-align; A3 structural divergence in
the pocket shell vs elsewhere (inactive-vs-active against active-vs-active); A4 vertebrate orthologue conservation.
Outputs: data/derived/followup_a1_alt_nucleophile.csv, a2_catalytic_positions.csv, a3_divergence.csv, a3_segments.csv,
a4_orthologue_conservation.csv, a4_elbow_residue_counts.csv; data/raw/orthologues/*.tsv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import requests
from Bio.Align import PairwiseAligner, substitution_matrices
from Bio.PDB import PDBParser
from scipy.spatial import cKDTree

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from abhd import config as C
from abhd.catalytic import triad_and_oxyanion
from abhd.structops import align, load

ORTH = C.RAW / "orthologues"
OUT = C.DERIVED
PAIRS = [("ABHD15", "ABHD3"), ("ABHD15", "ABHD1"), ("ABHD15", "ABHD2"),            # inactive vs active
         ("ABHD5", "ABHD4"),
         ("ABHD1", "ABHD3"), ("ABHD1", "ABHD2"), ("ABHD2", "ABHD3")]               # active vs active (baseline)
DIVERGED = 3.0
SHELL = 6.0
EXTRA = {"ABHD15": {"his_aligned": 391, "acid_aligned": 360, "ser_i+1_check": 238}, "ABHD5": {"his_aligned": 327, "acid_aligned": 177}}


def a1(cs):
    rows = []
    for r in cs.itertuples():
        if pd.isna(r.elbow_resnum):
            continue
        path = C.DERIVED / "domains" / f"{r.accession}.pdb"
        ids, _, seq = load(r.accession)
        pos = {n: s for n, s in zip(ids, seq)}
        base = str(r.category).startswith("no nucleophile")
        for off in range(-3, 4):
            n = int(r.elbow_resnum) + off
            if (off == 0 and not base) or pos.get(n) not in tuple("SCTD"):
                continue
            t = triad_and_oxyanion(path, n, r.accession)
            rows.append(dict(gene=r.gene, accession=r.accession, no_nucleophile_gene=base, offset=off, resnum=n, aa=pos[n],
                             triad=bool(t["n_order_ok_triads"]), his=t["his"], acid=t["acid"]))
    df = pd.DataFrame(rows); df.to_csv(OUT / "followup_a1_alt_nucleophile.csv", index=False)
    nn = df[df.no_nucleophile_gene]
    print("A1 no-nucleophile genes with an off-register S/C/T/D triad:", nn[nn.triad].gene.unique().tolist())
    ctl = df[~df.no_nucleophile_gene]
    print(f"A1 control: off-register S/C/T/D candidates in nucleophile genes: {len(ctl)}; with an ordered triad: {int(ctl.triad.sum())}"
          f" ({ctl.triad.mean():.2f})")
    return df


def a2_a3(cs, uni):
    gacc = dict(zip(uni.gene, uni.accession))
    crow = cs.set_index("gene")
    cat_rows, div_rows, seg_rows = [], [], []
    for g1, g2 in PAIRS:
        a, b = gacc[g1], gacc[g2]
        (ia, ca, sa), (ib, cb, sb), X, pairs, r = align(a, b)
        # reference = g2 (the active/second protein); pocket shell from its cavity points
        cav = np.load(C.DERIVED / "domains" / f"{b}.cavity.npy")
        t = cKDTree(cav)
        dist_pocket = t.query(cb)[0]
        dd = np.full(len(ib), np.nan); partner = {}
        for i, j in pairs:
            dd[j] = np.linalg.norm(X[i] - cb[j]); partner[j] = i
        inshell = dist_pocket <= SHELL
        unal = np.isnan(dd)
        div = unal | (dd > DIVERGED)
        kind = "inactive_vs_active" if g1 in ("ABHD15", "ABHD5") else "active_vs_active"
        div_rows.append(dict(pair=f"{g1} vs {g2}", kind=kind, tm_rmsd=round(r.rmsd, 2), n_aligned=len(pairs),
                             len_a=len(ia), len_ref=len(ib), n_shell=int(inshell.sum()),
                             diverged_in_shell=round(div[inshell].mean(), 3), diverged_outside=round(div[~inshell].mean(), 3),
                             unaligned_in_shell=round(unal[inshell].mean(), 3), unaligned_outside=round(unal[~inshell].mean(), 3),
                             mean_disp_shell=round(float(np.nanmean(dd[inshell])), 2), mean_disp_outside=round(float(np.nanmean(dd[~inshell])), 2)))
        # divergent segments (>=4 consecutive reference residues) overlapping the pocket shell
        j = 0
        while j < len(ib):
            if div[j]:
                k = j
                while k + 1 < len(ib) and div[k + 1]:
                    k += 1
                if k - j + 1 >= 4 and inshell[j:k + 1].any():
                    seg_rows.append(dict(pair=f"{g1} vs {g2}", ref_start=int(ib[j]), ref_end=int(ib[k]), ref_seq=sb[j:k + 1],
                                         n_shell_residues=int(inshell[j:k + 1].sum()), unaligned_frac=round(float(unal[j:k + 1].mean()), 2)))
                j = k + 1
            else:
                j += 1
        # catalytic positions of the reference
        ref = crow.loc[g2]
        for role, num in (("nucleophile", ref.elbow_resnum), ("His", ref.his), ("acid", ref.acid)):
            if pd.isna(num):
                continue
            jj = int(np.where(ib == int(num))[0][0]); ii = partner.get(jj)
            cat_rows.append(dict(pair=f"{g1} vs {g2}", role=role, ref_residue=f"{sb[jj]}{int(num)}",
                                 partner_residue=f"{sa[ii]}{int(ia[ii])}" if ii is not None else "unaligned",
                                 ca_displacement=round(float(dd[jj]), 2) if not np.isnan(dd[jj]) else np.nan))
    D = pd.DataFrame(div_rows); D.to_csv(OUT / "followup_a3_divergence.csv", index=False); print(D.to_string(index=False))
    C2 = pd.DataFrame(cat_rows); C2.to_csv(OUT / "followup_a2_catalytic_positions.csv", index=False); print(C2.to_string(index=False))
    S = pd.DataFrame(seg_rows); S.to_csv(OUT / "followup_a3_segments.csv", index=False); print(S.to_string(index=False))


def fetch_orthologues(gene):
    ORTH.mkdir(parents=True, exist_ok=True)
    p = ORTH / f"{gene}.tsv"
    if not p.exists():
        r = requests.get("https://rest.uniprot.org/uniprotkb/search", timeout=60, params={
            "query": f"gene:{gene} AND taxonomy_id:7742", "format": "tsv", "size": 500,
            "fields": "accession,organism_name,length,reviewed,sequence"})
        r.raise_for_status(); p.write_text(r.text)
    return pd.read_csv(p, sep="\t")


def a4(cs, uni):
    gacc = dict(zip(uni.gene, uni.accession)); crow = cs.set_index("gene")
    al = PairwiseAligner(mode="local", substitution_matrix=substitution_matrices.load("BLOSUM62"), open_gap_score=-10, extend_gap_score=-0.5)
    rows, counts = [], []
    for gene in ("ABHD15", "ABHD3", "ABHD1", "ABHD2", "ABHD4", "ABHD5"):
        ids, _, hseq = load(gacc[gene])
        orth = fetch_orthologues(gene)
        orth = orth[orth["Entry"] != gacc[gene]]
        cols = []                                  # per orthologue: aa aligned to each human position (or '-')
        orgs = set()
        for r in orth.itertuples():
            org = r.Organism
            if org in orgs or not isinstance(r.Sequence, str):
                continue
            a = al.align(hseq, r.Sequence)[0]
            m = ["-"] * len(hseq); nid = 0; nal = 0
            for (s0, e0), (s1, e1) in zip(*a.aligned):
                for k in range(e0 - s0):
                    m[s0 + k] = r.Sequence[s1 + k]; nal += 1; nid += hseq[s0 + k] == r.Sequence[s1 + k]
            if nal / len(hseq) < 0.7 or nid / nal < 0.4:
                continue
            orgs.add(org); cols.append("".join(m))
        n = len(cols)
        if n < 6:
            print(f"A4 {gene}: only {n} orthologues - skipped"); continue
        M = np.array([list(c) for c in cols])
        cons = np.array([np.mean(M[:, k][M[:, k] != "-"] == hseq[k]) if (M[:, k] != "-").any() else np.nan for k in range(len(hseq))])
        el = int(crow.loc[gene, "elbow_resnum"]); ie = int(np.where(ids == el)[0][0])
        keys = {f"elbow{o:+d}": ie + o for o in range(-2, 3)}
        for role in ("his", "acid"):
            v = crow.loc[gene, role]
            if pd.notna(v) and int(v) in ids:
                keys[role] = int(np.where(ids == int(v))[0][0])
        for role, num in EXTRA.get(gene, {}).items():       # His/acid positions that TM-align maps onto the active relatives'
            keys[role] = int(np.where(ids == num)[0][0])
        for name, k in keys.items():
            pct = float((cons[~np.isnan(cons)] < cons[k]).mean())
            rows.append(dict(gene=gene, n_orthologues=n, position=name, human=f"{hseq[k]}{int(ids[k])}", conservation=round(float(cons[k]), 3),
                             domain_median=round(float(np.nanmedian(cons)), 3), percentile_in_domain=round(pct, 3)))
        counts.append(dict(gene=gene, n_orthologues=n, human_elbow=f"{hseq[ie]}{el}",
                           elbow_residues=";".join(f"{a}:{c}" for a, c in zip(*np.unique(M[:, ie][M[:, ie] != '-'], return_counts=True))),
                           next_residue=";".join(f"{a}:{c}" for a, c in zip(*np.unique(M[:, ie + 1][M[:, ie + 1] != '-'], return_counts=True)))))
    df = pd.DataFrame(rows); df.to_csv(OUT / "followup_a4_orthologue_conservation.csv", index=False)
    pd.DataFrame(counts).to_csv(OUT / "followup_a4_elbow_residue_counts.csv", index=False)
    print(df.to_string(index=False)); print(pd.DataFrame(counts).to_string(index=False))


def a5(cs, uni):
    """geometry of the putative alternative nucleophile (ABHD15 Ser238) vs the retained His; Fisher test of shell divergence."""
    from scipy.stats import fisher_exact
    gacc = dict(zip(uni.gene, uni.accession)); P = PDBParser(QUIET=True)
    s15 = P.get_structure("a", str(C.DERIVED / "domains" / f"{gacc['ABHD15']}.pdb"))
    r = {x.get_id()[1]: x for x in s15.get_residues()}
    d = np.linalg.norm(r[238]["OG"].coord - r[391]["NE2"].coord), np.linalg.norm(r[238]["OG"].coord - r[391]["ND1"].coord)
    dg = np.linalg.norm(r[237]["CA"].coord - r[391]["NE2"].coord)
    s3 = P.get_structure("b", str(C.DERIVED / "domains" / f"{gacc['ABHD3']}.pdb")); r3 = {x.get_id()[1]: x for x in s3.get_residues()}
    d3 = np.linalg.norm(r3[220]["OG"].coord - r3[375]["NE2"].coord)
    print(f"A5 ABHD15 Ser238 OG to His391 NE2/ND1: {d[0]:.2f}/{d[1]:.2f} A (ABHD3 Ser220 OG to His375 NE2: {d3:.2f}); Gly237 CA to His391 NE2: {dg:.2f}")
    rng = np.random.default_rng(5)

    def shell_test(g_other, g_ref):
        """odds ratio of divergence inside the reference's pocket shell, Fisher p (residues independent) and a circular-shift p
        that keeps run structure (2000 shifts of the divergence mask along the chain; two-sided)."""
        a, b = gacc[g_other], gacc[g_ref]
        (ia, ca, sa), (ib, cb, sb), X, pairs, rr = align(a, b)
        cav = np.load(C.DERIVED / "domains" / f"{b}.cavity.npy"); sh = cKDTree(cav).query(cb)[0] <= SHELL
        dd = np.full(len(ib), np.inf)
        for i, j in pairs:
            dd[j] = np.linalg.norm(X[i] - cb[j])
        div = dd > DIVERGED
        tab = [[int(div[sh].sum()), int((~div[sh]).sum())], [int(div[~sh].sum()), int((~div[~sh]).sum())]]
        o, p = fisher_exact(tab)
        obs = div[sh].sum()
        null = np.array([np.roll(div, int(k))[sh].sum() for k in rng.integers(1, len(div), 2000)])
        circ = min(1.0, 2 * min((np.sum(null >= obs) + 1) / 2001, (np.sum(null <= obs) + 1) / 2001))
        return tab, o, p, circ

    rows = []
    for g1, g2 in PAIRS:
        for other, ref in ((g1, g2), (g2, g1)):
            tab, o, p, circ = shell_test(other, ref)
            rows.append(dict(pair=f"{g1} vs {g2}", reference=ref, shell_div=tab[0][0], shell_ok=tab[0][1], out_div=tab[1][0], out_ok=tab[1][1],
                             odds_ratio=round(o, 2), fisher_p=round(p, 4), circular_shift_p=round(circ, 4)))
    F = pd.DataFrame(rows); F.to_csv(OUT / "followup_a5_shell_fisher.csv", index=False); print(F.to_string(index=False))


def main():
    cs = pd.read_csv(C.DERIVED / "catalytic_sites.csv"); uni = pd.read_csv(C.DERIVED / "universe.csv")
    a1(cs); a2_a3(cs, uni); a4(cs, uni); a5(cs, uni)


if __name__ == "__main__":
    main()
