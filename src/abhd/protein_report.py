"""One-page HTML report for a protein: what the study knows about it and how it relates to the other 109.

    from abhd.atlas import Atlas; from abhd.protein_report import build_report
    build_report(Atlas(), "ABHD2", "out/protein_reports")

Sections: identity, catalytic site, neighbours in several views, function and label quality, position in the family, optional orthologue conservation,
and an optional superposed 3D view. Every number comes from the study's derived tables; a guardrails box states what the study found about their limits.
"""
import html
import json
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs, get_plotlyjs_version

from . import config as C
from .atlas import CORE_VIEWS, IDENTITY_FLOOR, IDENTITY_GOOD, Atlas

INK, ORANGE, BLUE, GREY = "#1B2A41", "#D9731A", "#2F6DA8", "#B8BEC6"
FONT = dict(family="Arial, Helvetica, sans-serif", size=13, color=INK)
PLOT_CFG = {"displaylogo": False, "responsive": True}

CSS = """body{font-family:Arial,Helvetica,sans-serif;color:#1B2A41;max-width:1100px;margin:24px auto;padding:0 16px;line-height:1.4}
h1{margin-bottom:0}h2{border-bottom:2px solid #D9731A;padding-bottom:4px;margin-top:34px}.sub{color:#5B6677;margin-top:2px}
table{border-collapse:collapse;margin:8px 0;font-size:14px}th,td{border-bottom:1px solid #ddd;padding:4px 10px;text-align:left}th{background:#F0F2F5}
.kv td:first-child{font-weight:bold;color:#5B6677;white-space:nowrap}.box{background:#F0F2F5;border-left:5px solid #D9731A;padding:8px 14px;margin:10px 0}
.warn{border-left-color:#B03A2E}.good{border-left-color:#2F6DA8}small{color:#5B6677}img{max-width:100%}th.s{cursor:pointer}th.s:hover{background:#e3e7ee}a{color:#2F6DA8}"""


def _fmt(x, nd=2):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "—"
    return f"{x:.{nd}f}" if isinstance(x, (float, np.floating)) else str(x)


def _table(df, index=False, classes=""):
    return df.to_html(index=index, classes=classes, border=0, na_rep="—", escape=True, float_format=lambda v: f"{v:.2f}")


def _kv(rows):
    return "<table class=kv>" + "".join(f"<tr><td>{html.escape(str(k))}</td><td>{v}</td></tr>" for k, v in rows) + "</table>"


def _plot(fig, post_script=None):
    """Plotly figure as an HTML fragment; plotly.js itself is loaded once per page from plotly.min.js next to the report."""
    return fig.to_html(full_html=False, include_plotlyjs=False, config=PLOT_CFG, post_script=post_script)


SORT_JS = """<script>document.querySelectorAll('table.sortable').forEach(t=>{t.querySelectorAll('th').forEach((th,i)=>{th.classList.add('s');th.title='click to sort';
th.onclick=()=>{const rows=[...t.tBodies[0].rows];const asc=th.dataset.asc!=='1';rows.sort((a,b)=>{const x=a.cells[i].innerText.trim(),y=b.cells[i].innerText.trim();
const nx=parseFloat(x),ny=parseFloat(y);const e=(v,n)=>v===''||v==='—'||isNaN(n);if(e(x,nx)&&e(y,ny))return 0;if(e(x,nx))return 1;if(e(y,ny))return -1;return (nx-ny)*(asc?1:-1)});
th.dataset.asc=asc?'1':'0';rows.forEach(r=>t.tBodies[0].appendChild(r));};});});</script>"""


# ---------- figures ----------------------------------------------------------------------------------------------------------
def fig_relations(A, gene, known=()):
    """Every other protein as a point: identity to `gene` (x) against structural distance (y). Hover for details; click opens that protein's report."""
    a = A.acc(gene)
    d = A.dist("structure").loc[a].drop(a).dropna()
    ident = A.identity_matrix().loc[a, d.index] * 100
    rows = []
    for b in d.index:
        g = A.gene(b)
        rows.append(dict(gene=g, x=ident[b], y=d[b], abhd=A.abhd[b], fam=A.same_family(a, b), ec3=", ".join(sorted(A.ec3[b])) or "none",
                         cat=str(A.cs.loc[b, "category"]) if b in A.cs.index else "—"))
    T = pd.DataFrame(rows)

    def hover(r):
        return (f"<b>{r.gene}</b>{'  (ABHD-named)' if r.abhd else ''}<br>identity {r.x:.0f}% · structural distance {r.y:.2f}<br>EC3: {r.ec3}<br>{r.cat}"
                f"{'<br>same paralog group / FunFam' if r.fam else ''}{'<br><i>click to open its report</i>' if r.gene in known else ''}")
    fig = go.Figure()
    fig.add_vrect(x0=0, x1=IDENTITY_FLOOR * 100, fillcolor="#B8BEC6", opacity=0.18, line_width=0, layer="below")
    fig.add_annotation(x=IDENTITY_FLOOR * 50, y=1.0, yref="paper", text="no better than guessing<br>(EC transfer)", showarrow=False, font=dict(size=11, color="#5B6677"), yanchor="top")
    for flag, name, col in ((False, "other names", BLUE), (True, "ABHD-named", ORANGE)):
        S = T[T.abhd == flag]
        fig.add_trace(go.Scatter(x=S.x, y=S.y, mode="markers", name=name, marker=dict(color=col, size=9, opacity=0.85), customdata=S.gene, hovertext=[hover(r) for r in S.itertuples()], hoverinfo="text"))
    F = T[T.fam]
    fig.add_trace(go.Scatter(x=F.x, y=F.y, mode="markers", name="same paralog group / FunFam", marker=dict(size=17, color="rgba(0,0,0,0)", line=dict(color=INK, width=1.5)), hoverinfo="skip", customdata=F.gene))
    top = T.sort_values("y").head(8)
    fig.add_trace(go.Scatter(x=top.x, y=top.y, mode="text", text=top.gene, textposition="top right", textfont=dict(size=11), showlegend=False, hoverinfo="skip"))
    fig.update_layout(font=FONT, height=480, margin=dict(l=60, r=20, t=50, b=50), plot_bgcolor="white", legend=dict(orientation="h", y=1.12, x=0),
                      xaxis=dict(title=f"sequence identity to {A.gene(a)} (%)", gridcolor="#eee", range=[0, max(60, T.x.max() + 3)]), yaxis=dict(title="structural distance (1 − TM-score)", gridcolor="#eee"))
    post = ("var el=document.getElementById('{plot_id}');var known=" + json.dumps(sorted(known)) + ";el.on('plotly_click',function(d){var p=d.points[0];"
            "if(p.customdata&&known.indexOf(p.customdata)>=0){window.location.href=p.customdata+'.html';}});")
    return _plot(fig, post)


def fig_conservation(cons, info):
    """Conservation across vertebrate orthologues along the domain, catalytic positions marked; hover shows residue and values."""
    p = cons["per_position"]
    fig = go.Figure(go.Scatter(x=p.resnum, y=p.conservation, mode="lines", line=dict(color=BLUE, width=1.5), name="conservation",
                               hovertext=[f"<b>{r.aa}{r.resnum}</b><br>conserved in {r.conservation:.0%} of {r.n_aligned} orthologues<br>percentile within the domain {r.percentile:.0%}" for r in p.itertuples()], hoverinfo="text"))
    for role, key, col in (("nucleophile position", "elbow_resnum", ORANGE), ("His", "his", INK), ("acid", "acid", "#3C8C6E")):
        v = info.get(key)
        if v is not None and not (isinstance(v, float) and np.isnan(v)) and int(v) in set(p.resnum):
            r = p[p.resnum == int(v)].iloc[0]
            fig.add_vline(x=float(v), line=dict(color=col, dash="dash", width=1))
            fig.add_trace(go.Scatter(x=[r.resnum], y=[r.conservation], mode="markers+text", text=[f"{r.aa}{r.resnum}"], textposition="top center", marker=dict(size=10, color=col), name=role,
                                     hovertext=f"<b>{role}: {r.aa}{r.resnum}</b><br>conserved in {r.conservation:.0%}", hoverinfo="text"))
    fig.update_layout(font=FONT, height=300, margin=dict(l=60, r=20, t=20, b=50), plot_bgcolor="white", xaxis=dict(title="residue number", gridcolor="#eee"),
                      yaxis=dict(title="conserved in orthologues", range=[0, 1.12], tickformat=".0%", gridcolor="#eee"), legend=dict(orientation="h", y=1.18, x=0))
    return _plot(fig)


# ---------- sections ---------------------------------------------------------------------------------------------------------
def sec_identity(A, i):
    flags = []
    if i["n_tm"]: flags.append(f"{i['n_tm']} predicted transmembrane segment(s)")
    if i["signal_peptide"]: flags.append("signal peptide")
    rows = [("Gene / accession", f"<b>{i['gene']}</b> / {i['accession']}"), ("Name", html.escape(i["name"])), ("ABHD-named", "yes" if i["abhd_named"] else "no"),
            ("Source", {"both": "UniProt family query and CATH FunFam", "funfam_only": "CATH FunFam only", "uniprot_only": "UniProt family query only"}.get(i["source"], i["source"])),
            ("FunFam(s)", ", ".join(i["funfam_ids"]) or "none"), ("Paralog group", i["paralog_group"]),
            ("Hydrolase domain", f"residues {i['dom_start']}–{i['dom_end']} ({i['dom_len']} aa; {i['domain_source']}); full length {i['length']} aa"),
            ("AlphaFold confidence", f"mean pLDDT {i['mean_plddt']:.1f} over the domain"), ("Flags", "; ".join(flags) or "none"),
            ("PDB entries", f"{len(i['pdb_entries'])}" + (f" ({', '.join(i['pdb_entries'][:8])}{'…' if len(i['pdb_entries']) > 8 else ''})" if i["pdb_entries"] else ""))]
    return "<h2>1. Identity</h2>" + _kv(rows)


def sec_catalytic(A, i):
    a = i["accession"]; g = i["gene"]
    cat = i.get("category") or "no elbow call"
    out = ["<h2>2. Catalytic site (from the AlphaFold model)</h2>"]
    rows = [("Category", f"<b>{html.escape(str(cat))}</b>")]
    if i.get("elbow_resnum") is not None:
        rows += [("Elbow residue", f"{i['elbow_aa']}{int(i['elbow_resnum'])} in …{i['elbow_context']}… (φ {_fmt(i['phi'], 0)}°, ψ {_fmt(i['psi'], 0)}°, pLDDT {_fmt(i['plddt'], 0)})")]
    if i.get("his") is not None:
        rows += [("Triad", f"{i['elbow_aa']}{int(i['elbow_resnum'])} – His{int(i['his'])} – {str(i['acid_type']).title()}{int(i['acid'])} (nucleophile–His {_fmt(i['d_nuc_his'])} Å, His–acid {_fmt(i['d_his_acid'])} Å; "
                               f"sequence order nucleophile, acid, His: {'yes' if i['order_ok'] else 'no'})")]
    d2 = A.table("followup_b1_donor2_all.csv")
    if d2 is not None and (d2.gene == g).any():
        r = d2[d2.gene == g].iloc[0]
        rows += [("Oxyanion donors", f"first: backbone N of the residue after the nucleophile ({_fmt(r.donor1_nu_plus1_dist)} Å); second: {r.donor2_aa}{_fmt(r.donor2_resnum, 0)} in …{r.donor2_context}… ({_fmt(r.donor2_dist)} Å)"
                                    if isinstance(r.donor2_aa, str) else f"first: backbone N after the nucleophile ({_fmt(r.donor1_nu_plus1_dist)} Å); no second donor within 5.5 Å")]
    a1 = A.table("followup_a1_alt_nucleophile.csv")
    if str(cat).startswith("no nucleophile"):
        hit = a1[(a1.gene == g) & (a1.triad)] if a1 is not None else None
        rows += [("Alternative nucleophile", "none: no Ser, Cys, Thr or Asp within three residues of the elbow forms a His/acid triad" if hit is not None and not len(hit) else "see scan table")]
    cs = A.table("cavity_summary.csv")
    if cs is not None and (cs.gene == g).any():
        r = cs[cs.gene == g].iloc[0]
        rows += [("Active-site cavity", f"median {_fmt(r['median_3.2'], 0)} grid units at 3.2 Å clearance (range {_fmt(r['vmin_3.2'], 0)}–{_fmt(r['vmax_3.2'], 0)}), {_fmt(r['median_3.6'], 0)} at 3.6 Å "
                                         f"(range {_fmt(r['vmin_3.6'], 0)}–{_fmt(r['vmax_3.6'], 0)}) across 2 models × 3 orientations; {r.frac_box_truncated:.0%} of variants touch the grid box"
                                         " (then volumes are lower bounds)")]
    pl = A.table("pocket_lining.csv")
    if pl is not None and (pl.gene == g).any():
        res = pl[pl.gene == g]
        rows += [("Pocket-lining residues", f"{len(res)}: " + " ".join(f"{r.aa}{r.resnum}" for r in res.itertuples()))]
    out.append(_kv(rows))
    return "".join(out)


def neighbour_table(A, gene, k=6):
    a = A.acc(gene)
    views = [v for v in CORE_VIEWS if v in A.available_views()]
    ranks = {v: A.rank_in(gene, v) for v in views}
    pool = []
    for v in views:
        pool += list(A.neighbours(gene, v, k=k).accession)
    pool = list(dict.fromkeys(pool))
    rows = []
    for b in pool:
        r = dict(gene=A.gene(b), identity_pct=A.identity(a, b) * 100, abhd=("yes" if A.abhd[b] else ""), same_family=("yes" if A.same_family(a, b) else ""), ec3=";".join(sorted(A.ec3[b])))
        for v in views:
            r[f"rank_{v}"] = int(ranks[v][b]) if b in ranks[v].index else np.nan
        rows.append(r)
    T = pd.DataFrame(rows)
    T["_s"] = T[[f"rank_{v}" for v in views]].min(axis=1)
    T = T.sort_values(["_s", f"rank_{views[0]}"], na_position="last").drop(columns="_s").reset_index(drop=True)
    T["identity_pct"] = T.identity_pct.round(0).astype(int)
    for v in views:
        T[f"rank_{v}"] = T[f"rank_{v}"].astype("Int64")
    return T.rename(columns={"abhd": "ABHD-named", "same_family": "same paralog group/FunFam", "identity_pct": "identity %"})


def sec_neighbours(A, i, k=6, known=()):
    g = i["gene"]
    T = neighbour_table(A, g, k=k)
    Tl = T.copy()
    Tl["gene"] = [f'<a href="{x}.html">{x}</a>' if x in known else x for x in Tl.gene]
    out = ["<h2>3. Neighbours in several views</h2>",
           "<p>Rank 1 = nearest. Views: structure (TM-align on the trimmed domain), sequence (local alignment), FunFam (HMM-to-HMM), elbow window (31 residues around the nucleophile position, ungapped), "
           "local site (SOAP descriptor at the elbow). Identity is a global alignment on the trimmed domains. The study found that function transfer from a neighbour was no better than guessing "
           f"below about {IDENTITY_FLOOR:.0%} identity.</p>", Tl.to_html(index=False, classes="sortable", border=0, na_rep="—", escape=False)]
    cross = []
    for v in ("structure", "sequence"):
        n = A.neighbours(g, v, k=5, cross_family=True)
        cross.append(f"<b>{v}</b>: " + ", ".join(f"{r.gene} ({r.identity:.0%})" for r in n.itertuples()))
    out.append("<p>Nearest from other paralog groups and FunFams — " + "; ".join(cross) + "</p>")
    return "".join(out), T


def sec_function(A, i):
    g, a = i["gene"], i["accession"]
    rows = [("EC", html.escape(i["ec"] or "none") + (f" (EC3: {', '.join(i['ec3'])})" if i["ec3"] else "")),
            ("Evidence behind the EC label", {"experimental_only": "experimentally supported reactions only", "mixed": "mixed: experimental and inferred", "non_experimental": "curator or similarity inference only",
                                              "no_record": "EC without a catalytic-activity record", "no_complete_ec3": "no complete EC3 label"}.get(i["evidence_class"], i["evidence_class"])),
            ("Rhea reactions", i["n_rhea"]), ("GO molecular function", ", ".join(i["go_mf"][:6]) + ("…" if len(i["go_mf"]) > 6 else "") if i["go_mf"] else "none")]
    out = ["<h2>4. Function and label quality</h2>", _kv(rows)]
    notes = []
    for view in ("structure", "sequence", "funfam"):
        if view not in A.available_views():
            continue
        p = A.predict_ec(g, view)
        if p["gated"]:
            notes.append("<div class='box warn'><b>No function transferred:</b> this protein has no nucleophile-type residue at the elbow in the AlphaFold model "
                         f"({html.escape(p['category'])}); an inherited EC would be inconsistent with the site.</div>")
            break
        if "neighbour" in p:
            cls = {"above": "good", "borderline": "", "below_floor": "warn"}[p["tier"]]
            verdict = {"below_floor": f"identity is at or below {IDENTITY_FLOOR:.0%}, where transfer was no better than guessing in our leave-one-out test: a hint only",
                       "borderline": f"identity of {IDENTITY_FLOOR:.0%}–{IDENTITY_GOOD:.0%}: the gain over guessing was unresolved in this range (intervals include zero)",
                       "above": f"identity of {IDENTITY_GOOD:.0%} or more: the range where a gain over guessing was seen (weakly; paralogs dominate it)"}[p["tier"]]
            agree = ""
            if i["ec3"]:
                agree = f" — agrees with its own EC3" if p["ec3"] in i["ec3"] else " — does NOT match its own EC3"
            notes.append(f"<div class='box {cls}'>{view}: nearest labelled cross-family neighbour <b>{p['neighbour']}</b> (EC3 {p['ec3']}, {p['identity']:.0%} identity){agree}; {verdict}.</div>")
    h = A.table("function_hypotheses.csv")
    if h is not None and (h.gene == g).any():
        r = h[h.gene == g].iloc[0]
        notes.append(f"<div class='box'>Hypothesis table (3 nearest structural neighbours): predicted EC3 <b>{r.predicted_ec3}</b> (support {r.support}); {html.escape(str(r.neighbours))}; "
                     f"transfer applicable: {'yes' if r.transfer_applicable else 'no'}.</div>")
    return "".join(out + notes)


def sec_position(A, i, known=()):
    g, a = i["gene"], i["accession"]
    out = ["<h2>5. Position in the family</h2>"]
    rows = []
    cm = A.cluster_members(g)
    for name, members in cm.items():
        rows.append((f"20-cluster ({name} space)", (", ".join(members) or "alone") + f" <small>({len(members)} others)</small>"))
    base = A.uni.abhd_named.mean()
    for view in ("structure", "sequence"):
        if view in A.available_views():
            n = A.neighbours(g, view, k=19)
            rows.append((f"ABHD-named among its 19 nearest ({view})", f"{n.abhd_named.mean():.0%} (expected by chance {base:.0%})"))
    d = A.table("multiview_gene_discordance.csv")
    if d is not None and (d.gene == g).any():
        r = d[d.gene == g].iloc[0]
        parts = [f"{v.replace('_', ' ')}: {r[f'nn_{v}']} ({'agrees' if r[f'agree_{v}'] else 'differs'})" for v in ("structure", "sequence", "funfam") if isinstance(r.get(f"nn_{v}"), str)]
        rows.append(("Nearest cross-family neighbour and EC3 agreement", "; ".join(parts)))
    out.append(_kv(rows))
    out.append("<p><small>Hover for details; click a point to open that protein's report (when generated). Shaded: neighbour identity at which EC transfer was no better than guessing.</small></p>")
    out.append(fig_relations(A, g, known))
    return "".join(out)


def sec_conservation(A, i, cons):
    out = ["<h2>6. Orthologue conservation</h2>"]
    if "error" in cons:
        return "".join(out) + f"<p><i>{html.escape(cons['error'])}</i></p>"
    out.append(f"<p>{cons['n_orthologues']} vertebrate orthologues (best sequence per species, paralogs filtered from {cons['n_candidates']} candidates). "
               "Classes: " + ", ".join(f"{r.taxon_class} {r.n}" for r in cons["classes"].itertuples()) + ". Sets are dominated by birds and mammals.</p>")
    out.append(fig_conservation(cons, i))
    kp = cons["key_positions"]
    if len(kp):
        out.append(_table(kp))
    return "".join(out)


GUARD = """<h2>Guardrails</h2><div class=box><ul>
<li>Everything is from AlphaFold apo models; side-chain positions, cavities and the “no nucleophile” calls are predictions.</li>
<li>Cavity volumes depend on the model (trimmed domain or full length), the anchor and the grid; only contrasts that held in all variants were treated as findings.</li>
<li>Function transfer from a neighbour was no better than guessing for neighbours of 30% identity or less and only weakly better above 40%; the identity tiers mark where a gain was or was not seen, not a validated prediction.</li>
<li>The ABHD name marks a neighbourhood, not a clade: the nearest neighbours of an ABHD protein are usually ABHD, but close paralogs drive that signal.</li>
<li>The universe of 110 genes is a lower bound (11 further reviewed candidates exist), so “nearest” is relative to this set.</li>
</ul></div>"""


def build_report(A, gene, out_dir, orthologues=False, structure=True, k=6, known_reports=None, plotly_cdn=False):
    """Write <out_dir>/<GENE>.html (and, with structure=True, the 3D view files). Returns a dict with paths and key facts."""
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    i = A.info(gene); g = i["gene"]
    if plotly_cdn:                                           # single self-contained page (needs internet): plotly.js from jsDelivr
        plotly_src = f"https://cdn.jsdelivr.net/npm/plotly.js-dist-min@{get_plotlyjs_version()}/plotly.min.js"
    else:                                                    # shared local copy next to the reports (works offline)
        plotly_src = "plotly.min.js"
        (out_dir / "plotly.min.js").exists() or (out_dir / "plotly.min.js").write_text(get_plotlyjs())
    known = set(known_reports) if known_reports is not None else {p.stem for p in out_dir.glob("*.html") if not p.stem.endswith("_view")}
    known.add(g)
    neigh_html, T = sec_neighbours(A, i, k=k, known=known)
    parts = [f"<h1>{g}</h1><p class=sub>{html.escape(i['name'])} · {i['accession']} · generated from the ABHD study tables</p>",
             sec_identity(A, i), sec_catalytic(A, i), neigh_html, sec_function(A, i), sec_position(A, i, known)]
    cons = None
    if orthologues:
        cons = __import__("abhd.orthologs", fromlist=["conservation"]).conservation(A, g)
        parts.append(sec_conservation(A, i, cons))
    view_note = ""
    if structure:
        try:
            from . import viewer
            vp = viewer.build(A, g, list(T.gene[:3]), out_dir)
            view_note = f"<h2>{'7' if orthologues else '6'}. 3D view</h2><p><a href='{vp.name}'>{vp.name}</a> (interactive Mol* viewer: superposed neighbours, catalytic residues, cavity; works offline after <code>scripts/fetch_molstar.sh</code>, otherwise loads Mol* from a CDN); PyMOL: <code>{g}_view.pml</code> with the superposed models in <code>{g}_structures/</code>.</p>"
        except Exception as e:                                    # missing model files must not kill the report
            view_note = f"<h2>3D view</h2><p><i>not generated ({type(e).__name__}: {html.escape(str(e))})</i></p>"
    parts += [view_note, GUARD]
    p = out_dir / f"{g}.html"
    p.write_text(f"<!doctype html><meta charset=utf-8><title>{g}</title><style>{CSS}</style><script src=\"{plotly_src}\"></script>" + "".join(parts) + SORT_JS)
    return dict(path=p, gene=g, info=i, neighbours=T, conservation=cons)
