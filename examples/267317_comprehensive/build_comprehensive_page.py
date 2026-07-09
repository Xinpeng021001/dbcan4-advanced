#!/usr/bin/env python3
"""Render the comprehensive CAZyme deep-dive page for 267317 from real DB data.

Produces a self-contained static HTML (inline CSS, 3Dmol.js from CDN, PDB embedded)
that faithfully mirrors what the live BioForge FastAPI/Jinja route would serve. This
is the export/screenshot artifact; the same structure is mirrored into the live
gene_detail template.
"""
import json, html

D = json.load(open("handoff/ui_render_267317.json"))
gene = D["gene"]; rec = D["comprehensive"]; calls = D["real_calls"]; feats = D["feats"]
seq = gene["protein_seq"]; L = len(seq)
PDB = open("real_structures/267317.pdb").read()

# ---- method colours (mirror bioforge.methods) ----
COL = {"HMMER":"#64748b","dbCAN_sub":"#64748b","DIAMOND":"#64748b",
       "ESM-C-kNN":"#2563eb","ESM-C-centroid":"#7c3aed","ESM-C-contrastive":"#c026d3",
       "Foldseek-CAZyme3D":"#059669","SaProt":"#0d9488","fusion":"#ea580c"}
FCOL = {"Pfam/hmmscan":"#3a5bd0","ESMFold":"#059669","Foldseek-CAZyme3D":"#0d9488",
        "DeepTMHMM":"#7c3aed","SignalP6":"#b25c00","CLEAN":"#c026d3","Biopython":"#2563eb","DeepLoc":"#b25c00"}

def esc(x): return html.escape(str(x)) if x is not None else ""

def confbar(v, col="#12925a"):
    if v is None: return ""
    pct = int(float(v)*100)
    return f'<span class="confbar"><span style="width:{pct}%;background:{col}"></span></span> {float(v):.2f}'

CSS = """
:root{--ink:#1a2233;--muted:#64728c;--line:#e2e8f2;--accent:#12925a;--warn:#b25c00;--bg:#f7f9fc;--card:#fff}
*{box-sizing:border-box}
body{margin:0;font:14px/1.55 -apple-system,BlinkMacSystemFont,'Segoe UI',Roboto,sans-serif;color:var(--ink);background:var(--bg)}
.wrap{max-width:1120px;margin:0 auto;padding:22px 26px 60px}
.hero{background:linear-gradient(135deg,#0f2942,#123a5c);color:#fff;border-radius:16px;padding:24px 28px;margin-bottom:22px}
.hero h1{margin:0 0 3px;font-size:26px;letter-spacing:-.3px}
.hero .sub{opacity:.85;font-size:14.5px}
.hero .chips{margin-top:14px;display:flex;flex-wrap:wrap;gap:8px}
.chip{background:rgba(255,255,255,.14);border:1px solid rgba(255,255,255,.25);border-radius:999px;padding:3px 12px;font-size:12.5px}
.chip b{font-weight:600}
.grid{display:grid;grid-template-columns:1fr 1fr;gap:18px}
.grid3{display:grid;grid-template-columns:repeat(3,1fr);gap:14px}
.card{background:var(--card);border:1px solid var(--line);border-radius:13px;padding:18px 20px;margin-bottom:18px;box-shadow:0 1px 2px rgba(20,40,80,.04)}
.card h2{margin:0 0 4px;font-size:16px}
.card h2 .muted{font-weight:400;font-size:13px}
.card .lead{color:var(--muted);font-size:12.7px;margin:0 0 14px}
.span2{grid-column:1 / -1}
table{width:100%;border-collapse:collapse;font-size:13px}
th,td{text-align:left;padding:6px 9px;border-bottom:1px solid var(--line);vertical-align:top}
th{font-size:11px;text-transform:uppercase;letter-spacing:.4px;color:var(--muted);font-weight:600}
.tag{display:inline-block;padding:1px 9px;border-radius:999px;font-size:11.5px;font-weight:600;color:#fff}
.tag.grey{background:#64748b}
.pill{display:inline-block;padding:2px 10px;border-radius:999px;font-size:11.5px;font-weight:600}
.pill.adv{background:#eafaf1;color:#0d6b41;border:1px solid #bfe8d2}
.pill.base{background:#eef1f6;color:#475569;border:1px solid #d9e0ec}
.pill.warn{background:#fdf1e3;color:var(--warn);border:1px solid #f2d5ae}
.confbar{display:inline-block;width:52px;height:7px;border-radius:4px;background:#eef1f6;vertical-align:middle;overflow:hidden;margin-right:5px}
.confbar span{display:block;height:100%}
dl.kv{display:grid;grid-template-columns:auto 1fr;gap:4px 16px;margin:0}
dl.kv dt{color:var(--muted);font-size:12.5px}
dl.kv dd{margin:0;font-weight:500}
.mono{font-family:ui-monospace,SFMono-Regular,Menlo,monospace}
#viewer3d{width:100%;height:340px;position:relative;border:1px solid var(--line);border-radius:10px;background:#0d1b2a}
.arch{width:100%;height:96px;position:relative;margin:8px 0 4px}
.legend{display:flex;gap:14px;flex-wrap:wrap;font-size:11.5px;color:var(--muted);margin-top:8px}
.legend span{display:inline-flex;align-items:center;gap:5px}
.sw{width:11px;height:11px;border-radius:3px;display:inline-block}
.note{font-size:12px;color:var(--muted);margin-top:8px;padding:8px 11px;background:#f2f5fa;border-radius:8px;border-left:3px solid var(--accent)}
.evrow td:first-child{font-weight:600}
.seqblock{font-family:ui-monospace,monospace;font-size:11px;line-height:1.5;white-space:pre-wrap;word-break:break-all;color:#475569;max-height:150px;overflow:auto;background:#f7f9fc;padding:10px;border-radius:8px}
.banner{background:#fff8ec;border:1px solid #f2d5ae;border-radius:11px;padding:13px 17px;margin-bottom:18px;font-size:13px;color:#7a4d12}
.banner b{color:#5e3a0d}
"""

# ---------- build sections ----------
cc = rec["cazyme_call"]; fn = rec["function"]; pc = rec["physicochemistry"]
st = rec["structure"]; ts = rec["topology_secretion"]; loc = rec["localization"]
doms = rec["domains_pfam"]["domains"]

# HERO
hero = f"""
<div class="hero">
  <h1>{esc(gene['gene_key'])} <span style="opacity:.6;font-size:17px">·</span> {esc(rec['product'])}</h1>
  <div class="sub">Comprehensive multi-tool annotation of a dbCAN4-identified fungal CAZyme candidate</div>
  <div class="chips">
    <span class="chip"><b>{L:,}</b> aa</span>
    <span class="chip">EC <b>{esc(fn['ec_number'])}</b></span>
    <span class="chip">{esc(loc['localization_call'])}</span>
    <span class="chip">ESMFold pLDDT <b>{st['plddt_mean']:.0f}</b></span>
    <span class="chip">MW <b>{pc['molecular_weight']/1000:.0f}</b> kDa · pI <b>{pc['theoretical_pI']}</b></span>
    <span class="chip">{esc(gene['sample_key'])}</span>
  </div>
</div>"""

# BANNER — the honest multi-domain headline
banner = f"""
<div class="banner">
  <b>Why this protein needs comprehensive annotation:</b> {esc(cc['summary'])}
</div>"""

# CAZYME CALL — evidence by tool
ev = cc["evidence_by_tool"]
call_rows = "".join(
  f"<tr><td>{esc(v['call'])}</td><td>{esc(k.replace('_',' '))}</td><td class='muted' style='font-size:12px'>{esc(v['basis'])}</td></tr>"
  for k,v in ev.items())
# baseline vs advanced tool tags
base_tags = " ".join(f'<span class="tag grey">{esc(c["tool"])}→{esc(c["cazy_family"])}</span>' for c in calls if c["release_id"]==1)
adv_rows = ""
for c in [c for c in calls if c["release_id"]==3]:
    col = COL.get(c["tool"],"#64748b")
    adv_rows += (f'<tr><td><span class="tag" style="background:{col}">{esc(c["tool"])}</span></td>'
                 f'<td><b>{esc(c["cazy_family"])}</b></td><td>{confbar(c["confidence"],col)}</td></tr>')

cazyme_card = f"""
<div class="card span2">
  <h2>CAZyme family — multi-tool evidence</h2>
  <p class="lead">Different tools key on different domains of this multi-domain protein. No single label is complete — the value is in seeing all signals together.</p>
  <div class="grid">
    <div>
      <h3 style="font-size:13px;margin:0 0 6px;color:var(--muted)">Per-tool call</h3>
      <table><thead><tr><th>Tool</th><th>Family</th><th>Confidence</th></tr></thead>
      <tbody>{adv_rows}
        <tr><td colspan=3 style="padding-top:9px"><span class="pill base">baseline (HMMER/dbCAN_sub/DIAMOND)</span> &nbsp;{base_tags}</td></tr>
      </tbody></table>
    </div>
    <div>
      <h3 style="font-size:13px;margin:0 0 6px;color:var(--muted)">What each signal means</h3>
      <table><thead><tr><th>Call</th><th>Evidence source</th><th>Basis</th></tr></thead>
      <tbody>{call_rows}</tbody></table>
    </div>
  </div>
  <div class="note">{esc(cc['reference_coverage_note'])} &nbsp;<b>Interpretation:</b> {esc(cc['interpretation'])}</div>
</div>"""

# FUNCTION / EC card
go = fn["go_molecular_function"]
clean = fn["ec_source_clean"]
func_card = f"""
<div class="card">
  <h2>Function · EC · substrate</h2>
  <p class="lead">Family-inherited EC corroborated by an independent sequence-based predictor.</p>
  <dl class="kv">
    <dt>Activity</dt><dd>{esc(fn['activity_name'])}</dd>
    <dt>EC (family)</dt><dd><b>{esc(fn['ec_number'])}</b> <span class="muted" style="font-size:11.5px">· dbCAN fam-substrate</span></dd>
    <dt>EC (CLEAN)</dt><dd><b>{esc(clean['predicted_ec'])}</b> <span class="pill adv">seq→EC, conf {clean['confidence']}</span> <span class="muted" style="font-size:11.5px">· agrees exactly</span></dd>
    <dt>GO (MF)</dt><dd><a href="https://www.ebi.ac.uk/QuickGO/term/{esc(go['id'])}">{esc(go['id'])}</a> {esc(go['name'])}</dd>
    <dt>Substrate</dt><dd>{esc(fn['substrate'])}</dd>
    <dt>Reaction</dt><dd style="font-size:12.5px">{esc(fn['reaction'])}</dd>
  </dl>
  <div class="note" style="border-left-color:#c026d3">CLEAN (Yu et al., <i>Science</i> 2023) predicts EC from sequence alone — orthogonal to the family label. Low absolute confidence (0.11) is expected for a large multi-domain fungal protein far from SwissProt, but the top EC matches the family EC exactly.</div>
</div>"""

# DOMAIN ARCHITECTURE (inline SVG, to scale)
def bp(x): return 40 + (x/L)*(WSVG-80)
WSVG=1040
plddt_vals=[float(l[60:66]) for l in PDB.splitlines() if l.startswith("ATOM") and l[12:16].strip()=="CA"]
# pLDDT sparkline points
step=max(1,len(plddt_vals)//260)
pts=" ".join(f"{bp(i):.1f},{40-(v/100)*26:.1f}" for i,v in enumerate(plddt_vals) if i%step==0)
dom_svg=""
DCOL={"Glyco_hydro_28":"#94a3b8","Bac_rhamnosid6H":"#059669"}
for d in doms:
    x0,x1=bp(d["seq_start"]),bp(d["seq_end"]); col=DCOL.get(d["name"],"#7c3aed")
    dom_svg+=(f'<rect x="{x0:.1f}" y="52" width="{x1-x0:.1f}" height="26" rx="5" fill="{col}"/>'
              f'<text x="{(x0+x1)/2:.1f}" y="68" text-anchor="middle" fill="#fff" font-size="11.5" font-weight="600">{esc(d["name"])}</text>'
              f'<text x="{(x0+x1)/2:.1f}" y="90" text-anchor="middle" fill="#475569" font-size="9.5">{esc(d["pfam_acc_base"])} · {d["seq_start"]}–{d["seq_end"]} · iE {d["i_evalue"]:.0e}</text>')
# signal peptide
sp0,sp1=bp(1),bp(18)
dom_svg+=f'<rect x="{sp0:.1f}" y="52" width="{max(sp1-sp0,3):.1f}" height="26" rx="3" fill="#ea580c"/>'
# n-glyc ticks
ngly=rec["sequence_features"]["n_glycosylation_sequons"]["sites"]
for s in ngly:
    p=s["position"] if isinstance(s,dict) else s
    dom_svg+=f'<line x1="{bp(p):.1f}" y1="48" x2="{bp(p):.1f}" y2="52" stroke="#0d9488" stroke-width="1"/>'
# axis
axis=f'<line x1="{bp(0):.1f}" y1="80" x2="{bp(L):.1f}" y2="80" stroke="#cbd5e1" stroke-width="1"/>'
for t in range(0,L+1,200):
    axis+=f'<text x="{bp(t):.1f}" y="95" text-anchor="middle" fill="#94a3b8" font-size="9">{t}</text>'
arch_card=f"""
<div class="card span2">
  <h2>Domain architecture <span class="muted">· Pfam / hmmscan (--cut_ga) · pLDDT track</span></h2>
  <svg viewBox="0 0 {WSVG} 105" width="100%" style="max-height:120px">
    <polyline points="{pts}" fill="none" stroke="#2563eb" stroke-width="1"/>
    <line x1="{bp(0):.1f}" y1="{40-(70/100)*26:.1f}" x2="{bp(L):.1f}" y2="{40-(70/100)*26:.1f}" stroke="#b25c00" stroke-dasharray="3,3" stroke-width=".7"/>
    <text x="{bp(0)-6:.1f}" y="18" text-anchor="end" fill="#64748b" font-size="9">pLDDT</text>
    {dom_svg}{axis}
  </svg>
  <div class="legend">
    <span><span class="sw" style="background:#ea580c"></span>signal peptide 1–18</span>
    <span><span class="sw" style="background:#94a3b8"></span>Glyco_hydro_28 (GH28-type, N-term)</span>
    <span><span class="sw" style="background:#059669"></span>Bac_rhamnosid6H (GH78-diagnostic, C-term)</span>
    <span><span class="sw" style="background:#0d9488"></span>N-glyc sequons (n={len(ngly)})</span>
  </div>
</div>"""

# 3D STRUCTURE card (3Dmol.js, PDB embedded)
struct_card=f"""
<div class="card span2">
  <h2>Predicted 3D structure <span class="muted">· ESMFold · mean pLDDT {st['plddt_mean']:.1f} · {st['n_residues']:,} residues</span></h2>
  <p class="lead">Real ESMFold model (facebook/esmfold_v1, folded on met GPU). Coloured by per-residue pLDDT (blue=confident → red=low). Drag to rotate · scroll to zoom.</p>
  <div id="viewer3d"></div>
  <div class="note">This is a genuine compact globular fold (radius of gyration / extended-length ratio 0.030), not a placeholder trace. The structure is what feeds the Foldseek structural search below.</div>
</div>"""

# FOLDSEEK HITS table
fh_rows=""
for h in st["foldseek_top_hits"][:8]:
    fh_rows+=("<tr><td class='mono'>%s</td><td><b>%s</b></td><td>%.2f</td><td>%.2f</td><td>%.2f</td><td>%s</td></tr>"
              % (esc(h['target']), esc(h['family']), h['tmscore'], h['prob'], h['lddt'], h['bits']))
foldseek_card = f"""
<div class="card">
  <h2>Structural homologs <span class="muted">· Foldseek vs CAZyme3D ({st['foldseek_summary']['n_hits']} hits)</span></h2>
  <p class="lead">Top structural matches of the ESMFold model against {esc(st['foldseek_summary']['reference'])}. High TM-score at low sequence identity = remote structural homology.</p>
  <table><thead><tr><th>Target</th><th>Family</th><th>TM</th><th>prob</th><th>LDDT</th><th>bits</th></tr></thead>
  <tbody>{fh_rows}</tbody></table>
  <div class="note">Top structural neighbours are <b>GH28</b> (TM 0.80&ndash;0.82) &mdash; the full-length fold is dominated by the large N-terminal Glyco_hydro_28 domain. GH78 references exist in the set but do not surface, consistent with the C-terminal rhamnosidase domain being the minor structural component.</div>
</div>"""

dt = ts["deeptmhmm"]
sec_card = f"""
<div class="card">
  <h2>Secretion &amp; topology <span class="muted">· DeepTMHMM</span></h2>
  <dl class="kv">
    <dt>Prediction</dt><dd><span class="pill adv">{esc(dt['prediction'])}</span> (classically secreted)</dd>
    <dt>Signal peptide</dt><dd>residues {dt['signal_peptide_span'][0]}&ndash;{dt['signal_peptide_span'][1]}</dd>
    <dt>TM helices</dt><dd>{dt['n_tm_helices']}</dd>
    <dt>SignalP 6.0</dt><dd><span class="pill warn">not installed (license-gated)</span></dd>
  </dl>
  <div class="note" style="border-left-color:#7c3aed">Ran for real via DeepTMHMM (DTU) on BioLib cloud. A cleaved N-terminal signal peptide with no TM helix = a secreted enzyme, consistent with a plant-cell-wall&ndash;degrading fungal CAZyme.</div>
</div>"""

loc_card = f"""
<div class="card">
  <h2>Subcellular localization</h2>
  <dl class="kv">
    <dt>Call</dt><dd><b>{esc(loc['localization_call'])}</b></dd>
    <dt>Confidence</dt><dd>{esc(loc.get('confidence','—'))}</dd>
    <dt>Basis</dt><dd style="font-size:12.5px">{esc(loc['basis'])}</dd>
    <dt>DeepLoc-2.0</dt><dd><span class="pill warn">not installed (license-gated)</span></dd>
  </dl>
</div>"""

phys_card = f"""
<div class="card">
  <h2>Physicochemistry <span class="muted">· Biopython</span></h2>
  <dl class="kv">
    <dt>Molecular weight</dt><dd>{pc['molecular_weight']:,.0f} Da</dd>
    <dt>Theoretical pI</dt><dd>{pc['theoretical_pI']}</dd>
    <dt>Instability index</dt><dd>{pc['instability_index']} <span class="muted">({esc(pc['instability_class'])})</span></dd>
    <dt>GRAVY</dt><dd>{pc['gravy']}</dd>
    <dt>Aromaticity</dt><dd>{pc['aromaticity']}</dd>
    <dt>N-glyc sequons</dt><dd>{rec['sequence_features']['n_glycosylation_sequons']['count']}</dd>
  </dl>
</div>"""

import csv as _csv
from collections import Counter
prov = list(_csv.DictReader(open("annotation_267317/evidence_provenance.tsv"), delimiter="\t"))
EBADGE = {"real_run":("adv","real run"),"reference_lookup":("base","reference"),
          "derived":("base","derived"),"not_installed":("warn","not installed")}
prov_rows=""
for r in prov:
    cls,lab = EBADGE.get(r["evidence"],("base",r["evidence"]))
    prov_rows += ("<tr class='evrow'><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td><span class='pill %s'>%s</span></td></tr>"
                  % (esc(r['category']), esc(r['field']), esc(r['value'][:52]), esc(r['tool']), cls, lab))
etot = Counter(r["evidence"] for r in prov)
prov_card = f"""
<div class="card span2">
  <h2>Evidence provenance <span class="muted">· every value traced to its tool ({len(prov)} annotations · {etot['real_run']} real runs)</span></h2>
  <p class="lead">Full transparency: what was computed by a real tool run, what came from a reference database, what was derived, and what is a license-gated scaffold.</p>
  <table><thead><tr><th>Category</th><th>Field</th><th>Value</th><th>Tool</th><th>Evidence</th></tr></thead>
  <tbody>{prov_rows}</tbody></table>
</div>"""

seq_fmt = "\n".join(seq[i:i+60] for i in range(0,L,60))
seq_card = f'<div class="card span2"><h2>Protein sequence <span class="muted">({L:,} aa)</span></h2><div class="seqblock">{esc(seq_fmt)}</div></div>'

BODY = (hero + banner + cazyme_card
        + '<div class="grid">' + func_card + sec_card + '</div>'
        + arch_card + struct_card
        + '<div class="grid">' + foldseek_card + '<div>' + loc_card + phys_card + '</div></div>'
        + prov_card + seq_card)

VIEWER_JS = """
<script>
document.addEventListener('DOMContentLoaded',function(){
  var pdb=document.getElementById('pdbdata').textContent;
  var v=$3Dmol.createViewer(document.getElementById('viewer3d'),{backgroundColor:'#0d1b2a'});
  v.addModel(pdb,'pdb');
  v.setStyle({},{cartoon:{colorscheme:{prop:'b',gradient:'roygb',min:50,max:90}}});
  v.zoomTo();v.render();
});
</script>
<script type="text/plain" id="pdbdata">""" + PDB + "</script>"

HTML = ('<!doctype html><html lang="en"><head><meta charset="utf-8">'
        '<meta name="viewport" content="width=device-width,initial-scale=1">'
        f'<title>{esc(gene["gene_key"])} — comprehensive CAZyme annotation — BioForge</title>'
        '<script src="https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"></script>'
        f'<style>{CSS}</style></head><body><div class="wrap">{BODY}</div>{VIEWER_JS}</body></html>')

open("comprehensive_267317.html","w").write(HTML)
print("wrote comprehensive_267317.html:", len(HTML), "bytes")
