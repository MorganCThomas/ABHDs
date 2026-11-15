"""Superposed 3D view of a protein and its nearest neighbours: multi-structure files for PyMOL and a self-contained 3Dmol.js page."""
import io
import json
from pathlib import Path

import numpy as np
from Bio.PDB import PDBIO, PDBParser
from tmtools import tm_align

from . import config as C
from .structops import load


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
    # cavity points of the query (old grid cavity, if present)
    cav = C.DERIVED / "domains" / f"{q}.cavity.npy"
    pts = np.load(cav).round(2).tolist() if cav.exists() else []
    if len(pts) > 400:
        pts = [pts[i] for i in np.linspace(0, len(pts) - 1, 400).astype(int)]
    labels = [atlas.gene(a) for a in accs]
    html = f"""<!doctype html><meta charset=utf-8><title>{g} and neighbours</title>
<style>body{{font-family:Arial,Helvetica,sans-serif;margin:16px}}#v{{width:100%;height:620px;position:relative;border:1px solid #ccc}}label{{margin-right:14px}}</style>
<h3>{g} (blue) with its nearest structural neighbours superposed (grey); catalytic residues as sticks</h3>
<p>{' '.join(f'<label><input type=checkbox checked data-i={i}> {labels[i]}</label>' for i in range(len(labels)))} <label><input type=checkbox id=cav checked> cavity points</label></p>
<div id=v></div>
<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
const models={json.dumps(models)}, cat={json.dumps([list(cat[a].values()) for a in accs])}, pts={json.dumps(pts)};
const v=$3Dmol.createViewer("v",{{backgroundColor:"white"}});
models.forEach((m,i)=>{{v.addModel(m,"pdb");
  v.setStyle({{model:i}},{{cartoon:{{color:i==0?"#2F6DA8":"#b8bec6",opacity:i==0?1:0.7}}}});
  if(cat[i].length) v.addStyle({{model:i,resi:cat[i]}},{{stick:{{colorscheme:i==0?"orangeCarbon":"whiteCarbon",radius:0.25}}}});}});
let spheres=pts.map(p=>v.addSphere({{center:{{x:p[0],y:p[1],z:p[2]}},radius:0.5,color:"#D9731A",alpha:0.25}}));
document.querySelectorAll("input[data-i]").forEach(b=>b.onchange=()=>{{const i=+b.dataset.i;
  v.setStyle({{model:i}},b.checked?{{cartoon:{{color:i==0?"#2F6DA8":"#b8bec6",opacity:i==0?1:0.7}}}}:{{}});
  if(b.checked&&cat[i].length)v.addStyle({{model:i,resi:cat[i]}},{{stick:{{colorscheme:i==0?"orangeCarbon":"whiteCarbon",radius:0.25}}}});v.render();}});
document.getElementById("cav").onchange=e=>{{spheres.forEach(s=>v.removeShape(s));spheres=e.target.checked?pts.map(p=>v.addSphere({{center:{{x:p[0],y:p[1],z:p[2]}},radius:0.5,color:"#D9731A",alpha:0.25}})):[];v.render();}};
v.zoomTo({{model:0}});v.render();
</script>
<p style="font-size:13px;color:#555">Requires an internet connection for 3Dmol.js. Cavity points are the stage-5 grid cavity of the trimmed domain (a descriptor, not a binding-site calculation).</p>
"""
    p = out_dir / f"{g}_view.html"; p.write_text(html)
    return p
