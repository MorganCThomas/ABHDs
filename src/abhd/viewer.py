"""Superposed 3D view of a protein and its nearest neighbours: multi-structure files for PyMOL and a self-contained 3Dmol.js page."""
import base64
import io
import json
import shutil
from pathlib import Path

import molviewspec as mvs
import numpy as np
from Bio.PDB import PDBIO, PDBParser
from tmtools import tm_align

from . import config as C
from .structops import load


MOLSTAR_VERSION = "5.13.1"
VENDOR = C.ROOT / "data" / "vendor" / f"molstar-{MOLSTAR_VERSION}"
CDN = f"https://cdn.jsdelivr.net/npm/molstar@{MOLSTAR_VERSION}/build/viewer"


def molstar_assets(out_dir):
    """Script and stylesheet paths for the viewer: a local copy (copied next to the reports) if scripts/fetch_molstar.sh has been run, else the CDN."""
    if (VENDOR / "molstar.js").exists():
        d = Path(out_dir) / "assets" / "molstar"; d.mkdir(parents=True, exist_ok=True)
        for f in ("molstar.js", "molstar.css"):
            if not (d / f).exists():
                shutil.copy(VENDOR / f, d / f)
        return "assets/molstar/molstar.js", "assets/molstar/molstar.css"
    return f"{CDN}/molstar.js", f"{CDN}/molstar.css"


def mvs_state(atlas, accs, labels, models, cat, pts):
    """MolViewSpec state: query cartoon blue, neighbours grey; catalytic residues as orange sticks with labels/tooltips; cavity points as spheres."""
    b = mvs.create_builder()
    for k, (a, lab, txt) in enumerate(zip(accs, labels, models)):
        url = "data:chemical/x-pdb;base64," + base64.b64encode(txt.encode()).decode()
        st = b.download(url=url).parse(format="pdb").model_structure()
        rep = st.component(selector="polymer").representation(type="cartoon")
        rep.color(color="#2F6DA8" if k == 0 else "#B8BEC6"); rep.opacity(opacity=1.0 if k == 0 else 0.4)
        if k == 0:
            st.component(selector="polymer").tooltip(text=lab)
        for role, num in cat[a].items():
            comp = st.component(selector=mvs.ComponentExpression(auth_seq_id=int(num)))
            comp.representation(type="ball_and_stick", size_factor=1.6 if k == 0 else 1.0).color(color="#D9731A" if k == 0 else "#F2C79B")
            comp.tooltip(text=f"{lab} {role} {num}")
            if k == 0:
                comp.label(text=f"{role} {num}")
                if role == "nucleophile":
                    comp.focus(radius_extent=14)
    if pts:
        prim = b.primitives(opacity=0.35)
        for x, y, z in pts:
            prim.sphere(center=[x, y, z], radius=0.5, color="#D9731A")
    return b.get_state().dumps()


def model_path(acc):
    p = C.DERIVED / "domains" / f"{acc}.pdb"
    return p if p.exists() else C.STRUCTURES / f"{acc}.pdb"


def superposed_pdb(acc, ref_acc):
    """PDB text of the trimmed domain of `acc` moved onto `ref_acc` by TM-align (reference is returned unchanged)."""
    s = PDBParser(QUIET=True).get_structure(acc, str(model_path(acc)))
    if acc != ref_acc:
        _, ca1, s1 = load(acc); _, ca2, s2 = load(ref_acc)
        r = tm_align(ca1, ca2, s1, s2)
        for atom in s.get_atoms():
            atom.set_coord(atom.coord @ r.u.T + r.t)
    io_ = PDBIO(); io_.set_structure(s); buf = io.StringIO(); io_.save(buf)
    return buf.getvalue()


def catalytic_residues(atlas, acc):
    c = atlas.cs.loc[acc] if acc in atlas.cs.index else None
    out = {}
    if c is not None:
        for role, key in (("nucleophile", "elbow_resnum"), ("His", "his"), ("acid", "acid")):
            if not np.isnan(c[key]):
                out[role] = int(c[key])
    return out


def build(atlas, gene, neighbour_genes, out_dir):
    """Write <GENE>_structures/*.pdb, <GENE>_view.pml and <GENE>_view.html; returns the HTML path."""
    q = atlas.acc(gene); g = atlas.gene(q)
    out_dir = Path(out_dir); sd = out_dir / f"{g}_structures"; sd.mkdir(parents=True, exist_ok=True)
    accs = [q] + [atlas.acc(n) for n in neighbour_genes]
    models, cat = [], {}
    for a in accs:
        txt = superposed_pdb(a, q)
        (sd / f"{atlas.gene(a)}.pdb").write_text(txt)
        models.append(txt); cat[a] = catalytic_residues(atlas, a)
    # PyMOL script
    pml = ["bg_color white", "set cartoon_fancy_helices, 1"]
    for a in accs:
        n = atlas.gene(a)
        pml += [f"load {g}_structures/{n}.pdb, {n}", f"hide everything, {n}", f"show cartoon, {n}", f"color {'marine' if a == q else 'grey70'}, {n}"]
        if cat[a]:
            sel = "+".join(str(v) for v in cat[a].values())
            pml += [f"show sticks, {n} and resi {sel} and not name N+C+O", f"color {'orange' if a == q else 'wheat'}, {n} and resi {sel} and not name N+C+O"]
    pml += ["zoom", f"# query {g} blue; neighbours grey; catalytic residues (nucleophile, His, acid) as sticks"]
    (out_dir / f"{g}_view.pml").write_text("\n".join(pml) + "\n")
    # cavity points of the query (stage-5 grid cavity, if present)
    cav = C.DERIVED / "domains" / f"{q}.cavity.npy"
    pts = np.load(cav).round(2).tolist() if cav.exists() else []
    if len(pts) > 400:
        pts = [pts[i] for i in np.linspace(0, len(pts) - 1, 400).astype(int)]
    labels = [atlas.gene(a) for a in accs]
    state = mvs_state(atlas, accs, labels, models, cat, pts)
    js, css = molstar_assets(out_dir)
    legend = " ".join(f"<span style='color:{'#2F6DA8' if i == 0 else '#5B6677'}'><b>{labels[i]}</b></span>" for i in range(len(labels)))
    html = f"""<!doctype html><meta charset=utf-8><title>{g} and neighbours (Mol*)</title>
<link rel="stylesheet" href="{css}">
<style>body{{font-family:Arial,Helvetica,sans-serif;margin:16px;color:#1B2A41}}#app{{width:100%;height:660px;position:relative;border:1px solid #ccc}}</style>
<h3>{g} (blue) with its nearest structural neighbours superposed (grey); catalytic residues as sticks (orange)</h3>
<p>Models: {legend}. {'Orange translucent spheres: the stage-5 grid cavity (a descriptor, not a binding-site calculation). ' if pts else ''}Hover residues for details; use the Mol* panel (top right) to toggle structures.</p>
<div id="app"></div>
<script src="{js}"></script>
<script>
const mvsj = {json.dumps(state)};
molstar.Viewer.create('app', {{layoutIsExpanded: false, layoutShowControls: true, layoutShowRemoteState: false, layoutShowSequence: true, layoutShowLog: false,
  layoutShowLeftPanel: false, viewportShowExpand: true, viewportShowSelectionMode: false, viewportShowAnimation: false, collapseRightPanel: true}})
  .then(v => v.loadMvsData(mvsj, 'mvsj'));
</script>
<p style="font-size:13px;color:#555">Viewer: Mol* {MOLSTAR_VERSION} (local copy if present in assets/, otherwise loaded from jsDelivr, which needs internet).</p>
"""
    p = out_dir / f"{g}_view.html"; p.write_text(html)
    return p
