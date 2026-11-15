"""Loaders shared by scripts and notebooks. Notebooks should only read through these."""
import json

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform

from . import config as C


def load_entries():
    """UniProt entries as a DataFrame (one row per accession) with EC levels and an ABHD-name flag."""
    meta = pd.DataFrame(json.load(open(C.ENTRIES_JSON)))
    meta["reviewed"] = meta["reviewed"].eq("reviewed")
    for lv in (1, 2, 3):
        meta[f"ec_lv{lv}"] = meta["ec"].apply(
            lambda x: ";".join(sorted({".".join(z.strip().split(".")[:lv]) for z in x.split(";")})) if x else "")
    meta["abhd_named"] = meta["gene_name"].str.startswith("ABHD")
    return meta.set_index("accession", drop=False)


def reviewed_accession(meta, gene):
    return meta.loc[(meta.gene_name == gene) & meta.reviewed, "accession"].iloc[0]


def load_funfam():
    """FunFam metadata annotated with human gene names, and the HHalign E-value matrix."""
    ff = pd.read_csv(C.FUNFAM_META).set_index("funfam_id", drop=False)
    m = json.load(open(C.FF_UNIPROT_MAP))
    ff["human_genes"] = ff.human_uniprots.apply(lambda s: ";".join(m[a]["gene"] for a in s.split(";")))
    ff["human_names"] = ff.human_uniprots.apply(lambda s: ";".join(m[a]["name"] for a in s.split(";")))
    ff["has_abhd_named"] = ff.human_genes.str.contains(r"(?:^|;)ABHD")
    ev = pd.read_csv(C.FUNFAM_HH, index_col=0)
    return ff, ev


def evalue_to_distance(ev):
    """Symmetrise HHalign E-values (keep the better direction) and map to a distance: -log10(E), flipped."""
    e = np.minimum(ev.values, ev.values.T)
    s = -np.log10(np.clip(e, 1e-300, None))
    s = np.clip(s, 0, None)
    np.fill_diagonal(s, s.max())
    d = s.max() - s
    np.fill_diagonal(d, 0)
    return pd.DataFrame(d, index=ev.index, columns=ev.columns)


def load_protbert():
    z = np.load(C.PROTBERT_NPZ, allow_pickle=False)
    return list(z["accessions"]), z["X"]


def cosine_matrix(X, labels):
    return pd.DataFrame(squareform(pdist(X, metric="cosine")), index=labels, columns=labels)


def nearest(dist, acc, meta, n=5, reviewed_only=True):
    """Nearest neighbours of `acc` (excluding itself) from a distance DataFrame, with gene names."""
    d = dist[acc].drop(acc).sort_values()
    out = pd.DataFrame({"distance": d, "gene": meta.reindex(d.index).gene_name.values,
                        "reviewed": meta.reindex(d.index).reviewed.values})
    if reviewed_only:
        out = out[out.reviewed.fillna(False)]
    return out.head(n)


def load_elbow():
    return pd.read_csv(C.ELBOW_CSV)


def umap_coords(X, seed=42, n_neighbors=15):
    import umap
    return umap.UMAP(n_components=2, random_state=seed, n_neighbors=n_neighbors, metric="euclidean").fit_transform(X)
