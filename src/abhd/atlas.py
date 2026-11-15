"""Per-protein access to everything the study computed: one object that joins the derived tables for any gene.

    from abhd.atlas import Atlas
    A = Atlas(); A.info("ABHD2"); A.neighbours("ABHD2", "structure", k=10)

All inputs are the tables written by scripts 07-29 (`data/derived`, `data/raw/uniprot`); nothing is recomputed here.
"""
import difflib
import re
from functools import lru_cache

import numpy as np
import pandas as pd

from . import config as C

D = C.DERIVED

# view name -> distance table (gene x gene, accession-indexed). Views that are missing on disk are skipped.
VIEWS = {
    "structure": "dist_tmalign.csv",
    "structure_short": "dist_tmshort.csv",
    "structure_core": "dist_tmalign_core200.csv",
    "sequence": "dist_seqid.csv",
    "funfam": "dist_funfam.csv",
    "elbow_window": "dist_anchor_elbow.csv",
    "catalytic_windows": "dist_anchor_cat.csv",
    "local_site": "dist_soap.csv",
    "triad": "dist_triad.csv",
    "language_model": "dist_protbert.csv",
}
CORE_VIEWS = ["structure", "sequence", "funfam", "elbow_window", "local_site"]
IDENTITY_FLOOR = 0.30        # at or below this neighbour identity, EC transfer was no better than guessing (report 13)
IDENTITY_GOOD = 0.40         # the gain over guessing was only detected (FunFam, elbow window) for neighbours up to ~40% identity and above


def family_stem(gene):
    """Paralog group used throughout the study: ABHD17A -> ABHD17; CES1 -> CES."""
    m = re.match(r"ABHD\d+", gene)
    return m.group(0) if m else re.match(r"[A-Z]+", gene).group(0)


def ec_set(x, level=3):
    """Complete EC classes at `level` (incomplete ECs such as 3.-.- are not classes)."""
    if not isinstance(x, str) or not x.strip():
        return set()
    out = set()
    for e in x.split(";"):
        parts = e.strip().split(".")[:level]
        if e.strip() and len(parts) == level and "-" not in parts:
            out.add(".".join(parts))
    return out


class Atlas:
    def __init__(self):
        self.uni = pd.read_csv(D / "universe.csv")
        self.acc_of = dict(zip(self.uni.gene, self.uni.accession))
        self.gene_of = {a: g for g, a in self.acc_of.items()}
        self.qc = pd.read_csv(D / "structure_qc.csv").set_index("accession")
        self.cs = pd.read_csv(D / "catalytic_sites.csv").set_index("accession")
        self.prov = pd.read_csv(D / "label_provenance.csv").set_index("accession")
        self.ann = pd.read_csv(C.RAW / "uniprot" / "universe_annotations.tsv", sep="\t").set_index("Entry")
        self.u = self.uni.set_index("accession")
        self.ffs = {a: set(str(f).split(";")) - {"nan"} for a, f in self.u.funfam_ids.items()}
        self.stem = {a: family_stem(g) for a, g in self.gene_of.items()}
        self.abhd = self.u.abhd_named.to_dict()
        self.ec3 = {a: ec_set(self.ann.loc[a, "EC number"]) if a in self.ann.index else set() for a in self.u.index}
        self._optional = {}

    # ---- lookup ----------------------------------------------------------------------------------------------------------
    def acc(self, gene):
        """Accession for a gene symbol (case-insensitive) or an accession in the universe."""
        g = gene.strip()
        if g in self.gene_of:
            return g
        hits = {k.upper(): v for k, v in self.acc_of.items()}
        if g.upper() in hits:
            return hits[g.upper()]
        close = difflib.get_close_matches(g.upper(), list(hits), n=5, cutoff=0.6)
        raise KeyError(f"'{gene}' is not one of the {len(self.acc_of)} genes in the universe. Similar: {', '.join(close) or 'none'}")

    def gene(self, acc):
        return self.gene_of[acc]

    def same_family(self, a, b):
        """Same paralog group (name stem) or a shared FunFam."""
        return self.stem[a] == self.stem[b] or bool(self.ffs[a] & self.ffs[b])

    # ---- tables ----------------------------------------------------------------------------------------------------------
    def table(self, name, **kw):
        """Optional derived CSV (None if it is not on disk)."""
        if name not in self._optional:
            p = D / name
            self._optional[name] = pd.read_csv(p, **kw) if p.exists() else None
        return self._optional[name]

    @lru_cache(maxsize=None)
    def dist(self, view):
        f = VIEWS[view]
        p = D / f
        if not p.exists():
            raise FileNotFoundError(f"{f} not found (run the study scripts first)")
        M = pd.read_csv(p, index_col=0).reindex(index=self.u.index, columns=self.u.index)
        if M.isna().any().any():                      # partial coverage views (FunFam, triad): symmetrise on the finite value
            v = M.values; M = pd.DataFrame(np.fmin(v, v.T), index=M.index, columns=M.columns)
        return M

    def available_views(self):
        return [v for v, f in VIEWS.items() if (D / f).exists()]

    @lru_cache(maxsize=None)
    def identity_matrix(self):
        I = pd.DataFrame(np.nan, index=self.u.index, columns=self.u.index)
        t = pd.read_csv(D / "identity_pairs.csv")
        for a, b, i in t[["a", "b", "identity"]].itertuples(index=False):
            I.loc[a, b] = I.loc[b, a] = i
        for a in I.index:
            I.loc[a, a] = 1.0
        return I

    def identity(self, a, b):
        return float(self.identity_matrix().loc[a, b])

    # ---- per-protein summaries ---------------------------------------------------------------------------------------------
    def info(self, gene):
        """One dict with the identity, domain and catalytic-site facts for a gene."""
        a = self.acc(gene)
        u, q, c = self.u.loc[a], self.qc.loc[a], self.cs.loc[a] if a in self.cs.index else None
        ann = self.ann.loc[a] if a in self.ann.index else None
        out = dict(gene=self.gene_of[a], accession=a, name=str(u["name"]).split(" (EC")[0], abhd_named=bool(u.abhd_named), source=u.source,
                   funfam_ids=sorted(self.ffs[a]), paralog_group=self.stem[a], length=int(u.length))
        out.update(dom_start=int(q.dom_start), dom_end=int(q.dom_end), dom_len=int(q.dom_len), domain_source=q.domain_source, mean_plddt=float(q.dom_mean_plddt),
                   n_tm=int(q.n_tm), signal_peptide=bool(q.signal_peptide))
        if c is not None:
            out.update({k: (None if pd.isna(c[k]) else c[k]) for k in ("category", "elbow_aa", "elbow_resnum", "elbow_context", "phi", "psi", "plddt", "his", "acid", "acid_type",
                                                                      "d_nuc_his", "d_his_acid", "order_ok", "gold_nuc")})
        out["pdb_entries"] = [p for p in str(ann["PDB"]).split(";") if p.strip() and p != "nan"] if ann is not None else []
        out["ec"] = str(ann["EC number"]) if ann is not None and isinstance(ann["EC number"], str) else ""
        out["ec3"] = sorted(self.ec3[a])
        out["evidence_class"] = self.prov.loc[a, "evidence_class"] if a in self.prov.index else "unknown"
        out["n_rhea"] = int(self.prov.loc[a, "n_rhea"]) if a in self.prov.index else 0
        go = ann["Gene Ontology (molecular function)"] if ann is not None else None
        out["go_mf"] = [g.strip() for g in str(go).split(";") if g.strip() and g != "nan"] if isinstance(go, str) else []
        return out

    def neighbours(self, gene, view, k=10, cross_family=False, labelled_only=False):
        """Nearest other proteins in a view, as a DataFrame (accession, gene, distance, identity, flags)."""
        a = self.acc(gene)
        d = self.dist(view).loc[a].drop(a).dropna()
        rows = []
        for b, v in d.sort_values().items():
            if cross_family and self.same_family(a, b):
                continue
            if labelled_only and not self.ec3[b]:
                continue
            rows.append(dict(accession=b, gene=self.gene_of[b], distance=float(v), identity=self.identity(a, b), abhd_named=self.abhd[b],
                             same_family=self.same_family(a, b), ec3=";".join(sorted(self.ec3[b]))))
            if len(rows) >= k:
                break
        return pd.DataFrame(rows, columns=["accession", "gene", "distance", "identity", "abhd_named", "same_family", "ec3"])

    def rank_in(self, gene, view):
        """Rank (1 = nearest) of every other protein from the point of view of `gene`."""
        a = self.acc(gene)
        d = self.dist(view).loc[a].drop(a).dropna()
        return d.rank(method="min").astype(int)

    def predict_ec(self, gene, view="structure"):
        """EC3 of the nearest labelled cross-family neighbour, with its identity and a reliability flag (report 13)."""
        a = self.acc(gene)
        cat = str(self.cs.loc[a, "category"]) if a in self.cs.index else ""
        out = dict(view=view, gated=cat.startswith("no nucleophile"), category=cat)
        n = self.neighbours(gene, view, k=1, cross_family=True, labelled_only=True)
        if len(n):
            r = n.iloc[0]
            tier = "below_floor" if r.identity <= IDENTITY_FLOOR else ("borderline" if r.identity < IDENTITY_GOOD else "above")
            out.update(neighbour=r.gene, ec3=r.ec3, identity=float(r.identity), tier=tier)
        return out

    def cluster_members(self, gene):
        """Genes sharing a 20-cluster (average linkage) with `gene` in structure and sequence space."""
        t = self.table("global_cluster_assignments.csv")
        a = self.acc(gene)
        if t is None:
            return {}
        row = t[t.gene == self.gene_of[a]]
        out = {}
        for col, name in (("cl20_tmalign", "structure"), ("cl20_seqid", "sequence")):
            if col in t and len(row):
                members = sorted(t[t[col] == row[col].iat[0]].gene)
                out[name] = [m for m in members if m != self.gene_of[a]]
        return out
