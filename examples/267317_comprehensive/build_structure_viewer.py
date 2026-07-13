#!/usr/bin/env python3
"""
build_structure_viewer.py — generate a self-contained INTERACTIVE 3D structure viewer with a
colour-scheme drop-down the user drives themselves (no need to show every option at once).

The output HTML embeds the structure (PDB) and 3Dmol.js inline, so it works offline in any browser.
A <select> switches the live cartoon between four colouring schemes on the fly:
    domain    Pfam / CAZyme domains  (each domain a distinct colour, linker grey)   [default]
    plddt     per-residue ESMFold confidence  (blue confident -> red low)
    spectrum  sequence position  (N-terminus blue -> C-terminus red)
    sstruc    secondary structure  (helix / sheet / coil)
A live legend updates with the selected mode. Drag to rotate, scroll to zoom.

Inputs:
    --pdb PATH        the structure. --domains FILE (JSON list of {name,family,start,end,color});
                      defaults to the built-in 267317 spec.
    --threedmol PATH  local 3Dmol-min.js (else fetched from cdnjs and vendored inline).
    --title STR, --subtitle STR, --id STR
Output (‑‑out FILE, default structure_viewer_<id>.html): one self-contained HTML file.

Usage:
    python build_structure_viewer.py --pdb 267317.pdb --out structure_viewer_267317.html
"""
import argparse, json, os, urllib.request

CDN = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"
DEMO_DOMAINS_267317 = [
    {"name": "Signal peptide",  "family": "SP (DeepTMHMM)", "start": 1,   "end": 18,  "color": "#f59e0b"},
    {"name": "Glyco_hydro_28",  "family": "GH28 (PF00295)", "start": 62,  "end": 412, "color": "#3b82f6"},
    {"name": "Bac_rhamnosid6H", "family": "GH78 (PF17389)", "start": 667, "end": 886, "color": "#10b981"},
]

TEMPLATE = """<!doctype html><html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
  body{margin:0;font-family:ui-sans-serif,system-ui,-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;
       background:#f7f9fc;color:#0f172a}
  .wrap{max-width:960px;margin:24px auto;padding:0 16px}
  h1{font-size:20px;margin:0 0 2px}
  .sub{color:#475569;font-size:13px;margin:0 0 14px}
  .card{background:#fff;border:1px solid #e2e8f0;border-radius:12px;padding:16px;
        box-shadow:0 1px 3px rgba(15,23,42,.06)}
  .ctl{display:flex;align-items:center;gap:10px;margin-bottom:12px;flex-wrap:wrap}
  .ctl label{font-weight:600;font-size:13.5px}
  select{font-size:14px;padding:6px 10px;border:1px solid #cbd5e1;border-radius:8px;background:#fff;cursor:pointer}
  #viewer3d{position:relative;width:100%;height:560px;border-radius:10px;background:#0d1b2a}
  #legend{margin-top:12px;font-size:12.7px;color:#334155;line-height:1.9}
  .sw{display:inline-block;width:13px;height:13px;border-radius:3px;vertical-align:-2px;margin:0 5px 0 12px}
  .sw:first-child{margin-left:0}
  .hint{color:#64748b;font-size:12px;margin-top:8px}
</style>
<script>__THREEDMOL__</script>
</head><body>
<div class="wrap">
  <h1>__TITLE__</h1>
  <p class="sub">__SUBTITLE__</p>
  <div class="card">
    <div class="ctl">
      <label for="colormode">Colour by:</label>
      <select id="colormode">
        <option value="domain">Pfam / CAZyme domains</option>
        <option value="plddt">pLDDT confidence</option>
        <option value="spectrum">Sequence position (N&#8594;C)</option>
        <option value="sstruc">Secondary structure</option>
      </select>
    </div>
    <div id="viewer3d"></div>
    <div id="legend"></div>
    <p class="hint">Drag to rotate &middot; scroll to zoom &middot; pick a colour scheme above.</p>
  </div>
</div>
<script id="pdbdata" type="text/plain">__PDB__</script>
<script>
var DOMAINS = __DOMAINS__;
var BASE = "#9ca3af";
var viewer = null;
function sw(c){return '<span class="sw" style="background:'+c+'"></span>';}
function applyMode(mode){
  var lg = document.getElementById('legend');
  if(mode==='plddt'){
    viewer.setStyle({},{cartoon:{colorscheme:{prop:'b',gradient:'roygb',min:50,max:90}}});
    lg.innerHTML = 'Per-residue ESMFold pLDDT: '+sw('#2166ac')+'confident (\u226590) \u2192 '+sw('#b2182b')+'low (\u226450)';
  } else if(mode==='spectrum'){
    viewer.setStyle({},{cartoon:{color:'spectrum'}});
    lg.innerHTML = 'Sequence position: '+sw('#2166ac')+'N-terminus \u2192 '+sw('#b2182b')+'C-terminus';
  } else if(mode==='sstruc'){
    viewer.setStyle({},{cartoon:{colorscheme:'ssJmol'}});
    lg.innerHTML = 'Secondary structure: '+sw('#c8102e')+'helix '+sw('#ffd700')+'sheet '+sw('#e2e8f0')+'coil';
  } else { // domain
    viewer.setStyle({},{cartoon:{color:BASE}});
    var html = '';
    for(var i=0;i<DOMAINS.length;i++){
      var d = DOMAINS[i];
      viewer.setStyle({resi: d.start+'-'+d.end},{cartoon:{color:d.color}});
      html += sw(d.color)+d.family+' \u00b7 '+d.start+'\u2013'+d.end+'  ';
    }
    html += sw(BASE)+'inter-domain / linker';
    lg.innerHTML = html;
  }
  viewer.render();
}
document.addEventListener('DOMContentLoaded', function(){
  var pdb = document.getElementById('pdbdata').textContent;
  viewer = $3Dmol.createViewer(document.getElementById('viewer3d'), {backgroundColor:'#0d1b2a'});
  viewer.addModel(pdb,'pdb');
  var sel = document.getElementById('colormode');
  applyMode(sel.value);
  viewer.zoomTo();
  sel.addEventListener('change', function(e){ applyMode(e.target.value); });
});
</script></body></html>"""

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--domains", default="")
    ap.add_argument("--threedmol", default="")
    ap.add_argument("--title", default="267317 \u00b7 ESMFold structure")
    ap.add_argument("--subtitle", default="Multi-domain fungal \u03b1-L-rhamnosidase (GH28 + GH78) \u00b7 1,088 residues \u00b7 mean pLDDT 76.6")
    ap.add_argument("--id", default="267317")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    pdb = open(a.pdb, encoding="utf-8", errors="replace").read().strip()
    n_ca = sum(1 for l in pdb.splitlines() if l.startswith("ATOM") and l[12:16].strip() == "CA")
    domains = json.load(open(a.domains)) if a.domains else DEMO_DOMAINS_267317
    js = (open(a.threedmol, encoding="utf-8", errors="replace").read() if a.threedmol
          else urllib.request.urlopen(CDN, timeout=60).read().decode("utf-8", "replace"))
    out = a.out or ("structure_viewer_%s.html" % a.id)
    html = (TEMPLATE.replace("__THREEDMOL__", js)
                    .replace("__PDB__", pdb)
                    .replace("__DOMAINS__", json.dumps(domains))
                    .replace("__TITLE__", a.title)
                    .replace("__SUBTITLE__", a.subtitle))
    open(out, "w", encoding="utf-8").write(html)
    print("wrote", out, os.path.getsize(out), "B | CA atoms:", n_ca, "| domains:", len(domains),
          "| 3Dmol inline:", "3Dmol" in js[:5000])

if __name__ == "__main__":
    main()
