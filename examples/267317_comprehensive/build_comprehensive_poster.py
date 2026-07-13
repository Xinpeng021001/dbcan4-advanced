#!/usr/bin/env python3
"""
build_comprehensive_poster.py — render the full comprehensive 267317 CAZyme page to a single PNG,
INCLUDING the 3D structure, via headless Google Chrome.

The comprehensive page (comprehensive_267317.html) already carries every card (multi-tool CAZyme
evidence, function/EC, secretion, domain architecture, structural homologs, localization,
physicochemistry, evidence provenance) AND a full protein-sequence block AND a live 3Dmol WebGL
structure viewer. A naive full-page headless screenshot renders everything EXCEPT the WebGL viewer,
which paints on a canvas that the screenshot races. This script fixes that in two phases:

  Phase A  — render the ESMFold structure to a static PNG using 3Dmol's pngURI() framebuffer
             readback (needs preserveDrawingBuffer + a synchronous pngURI right after render()).
  Phase B  — splice that PNG into the page as a static <img> in place of the live viewer, then
             capture the whole page. The structure now shows reliably because no WebGL is needed
             at full-page-capture time.

Inputs  (‑‑html DIR, default the script's own directory):
    comprehensive_267317.html   the comprehensive page with cards+sequence+viewer and the ESMFold
                                PDB embedded inline in <script id="pdbdata">.
Requires: google-chrome-stable (or google-chrome) with SwiftShader WebGL, and Pillow.
Network:  fetches 3Dmol.js 2.1.0 once from cdnjs (met reaches it; if offline, pass ‑‑threedmol PATH).

Output  (‑‑outdir DIR, default '.'):
    comprehensive_267317_full.png       the merged poster (cards + real 3D structure + sequence)
    comprehensive_267317_static.html     the self-contained page with the structure baked in as <img>
    structure_267317.png                 the standalone structure render (kept for reuse)

Usage:
    python build_comprehensive_poster.py --html . --outdir ./out
"""
import argparse, base64, os, re, shutil, subprocess, sys, urllib.request

CDN = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"

def chrome():
    exe = shutil.which("google-chrome-stable") or shutil.which("google-chrome") or shutil.which("chromium")
    if not exe: sys.exit("no chrome/chromium on PATH")
    return exe

def get_3dmol(path):
    if path and os.path.exists(path):
        return open(path, encoding="utf-8", errors="replace").read()
    return urllib.request.urlopen(CDN, timeout=60).read().decode("utf-8", "replace")

def render_structure(pdb, js, chrome_exe, out_png):
    """Phase A: 3Dmol cartoon coloured by pLDDT -> static PNG via pngURI framebuffer readback."""
    page = ('<!doctype html><html><head><meta charset="utf-8"><script>' + js + '</script></head>'
            '<body style="margin:0"><div id="viewer" style="width:1400px;height:1150px;position:relative"></div>'
            '<div id="out">PENDING</div>'
            '<script id="pdb" type="text/plain">' + pdb + '</script>'
            '<script>function go(){var out=document.getElementById("out");try{'
            'var pdb=document.getElementById("pdb").textContent;'
            'var v=$3Dmol.createViewer(document.getElementById("viewer"),'
            '{backgroundColor:"#0d1b2a",preserveDrawingBuffer:true,antialias:true});'
            'v.addModel(pdb,"pdb");'
            'v.setStyle({},{cartoon:{colorscheme:{prop:"b",gradient:"roygb",min:50,max:90}}});'
            'v.zoomTo();v.render();out.textContent=v.pngURI();'   # synchronous readback
            '}catch(e){out.textContent="ERR:"+e;}}'
            'if(document.readyState!=="loading")go();else document.addEventListener("DOMContentLoaded",go);'
            '</script></body></html>')
    tmp = os.path.abspath("_render.html"); open(tmp, "w", encoding="utf-8").write(page)
    cmd = [chrome_exe, "--headless=new", "--no-sandbox", "--in-process-gpu", "--use-gl=angle",
           "--use-angle=swiftshader-webgl", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist",
           "--virtual-time-budget=25000", "--run-all-compositor-stages-before-draw",
           "--dump-dom", "file://" + tmp]
    dom = subprocess.run(cmd, capture_output=True, text=True, timeout=180).stdout
    m = re.search(r'id="out">(data:image/png;base64,[A-Za-z0-9+/=]+)<', dom)
    if not m:
        m2 = re.search(r'id="out">([^<]*)<', dom); sys.exit("structure render failed: " + (m2.group(1)[:200] if m2 else "?"))
    open(out_png, "wb").write(base64.b64decode(m.group(1).split(",", 1)[1]))
    # tighten crop around the molecule (navy bg #0d1b2a)
    from PIL import Image, ImageChops
    im = Image.open(out_png).convert("RGB")
    bbox = ImageChops.difference(im, Image.new("RGB", im.size, (13, 27, 42))).getbbox()
    if bbox and (bbox[2]-bbox[0]) < im.size[0]:
        pad = 50; l, t, rr, bb = bbox
        im.crop((max(0,l-pad), max(0,t-pad), min(im.size[0],rr+pad), min(im.size[1],bb+pad))).save(out_png)
    return out_png

def build(html_dir, outdir, threedmol_path):
    os.makedirs(outdir, exist_ok=True)
    html_in = os.path.join(html_dir, "comprehensive_267317.html")
    html = open(html_in, encoding="utf-8", errors="replace").read()
    js = get_3dmol(threedmol_path)
    ce = chrome()

    # extract the ESMFold PDB embedded in the page
    mpdb = re.search(r'<script[^>]*id="pdbdata"[^>]*>(.*?)</script>', html, re.S)
    if not mpdb: sys.exit("no embedded #pdbdata in HTML")
    pdb = mpdb.group(1).strip()
    n_ca = sum(1 for ln in pdb.splitlines() if ln.startswith("ATOM") and ln[12:16].strip() == "CA")
    print("embedded PDB CA atoms:", n_ca)

    # Phase A
    struct_png = render_structure(pdb, js, ce, os.path.join(outdir, "structure_267317.png"))
    from PIL import Image
    print("structure png:", os.path.getsize(struct_png), "B", Image.open(struct_png).size)
    b64 = base64.b64encode(open(struct_png, "rb").read()).decode("ascii")

    # Phase B: splice static <img> in place of the live WebGL viewer; drop CDN + init script
    img_tag = ('<img alt="ESMFold model of 267317 coloured by pLDDT" '
               'style="width:100%;border-radius:10px;background:#0d1b2a;display:block" '
               'src="data:image/png;base64,' + b64 + '">')
    div = '<div id="viewer3d"></div>'
    assert div in html, "viewer3d div not found"
    html = html.replace(div, img_tag, 1)
    html = re.sub(r'<script[^>]*src="[^"]*3Dmol[^"]*"[^>]*>\s*</script>', '', html, count=1, flags=re.I)
    html = re.sub(r"<script>\s*document\.addEventListener\('DOMContentLoaded'.*?</script>", '', html, count=1, flags=re.S)
    assert "cdnjs" not in html and "createViewer" not in html, "viewer refs not fully removed"
    html_static = os.path.join(outdir, "comprehensive_267317_static.html")
    open(html_static, "w", encoding="utf-8").write(html)
    print("static html:", os.path.getsize(html_static), "B")

    # Phase C: full-page capture (no WebGL needed now) + autocrop the page bg (#f7f9fc)
    out_png = os.path.join(outdir, "comprehensive_267317_full.png")
    cmd = [ce, "--headless=new", "--no-sandbox", "--hide-scrollbars",
           "--window-size=1200,7600", "--force-device-scale-factor=2",
           "--run-all-compositor-stages-before-draw", "--virtual-time-budget=15000",
           "--screenshot=" + out_png, "file://" + os.path.abspath(html_static)]
    subprocess.run(cmd, capture_output=True, text=True, timeout=180)
    if not os.path.exists(out_png): sys.exit("full-page capture failed")
    from PIL import Image, ImageChops
    im = Image.open(out_png).convert("RGB")
    bbox = ImageChops.difference(im, Image.new("RGB", im.size, (247, 249, 252))).getbbox()
    if bbox:
        pad = 48; l, t, rr, bb = bbox
        im.crop((max(0,l-pad), max(0,t-pad), min(im.size[0],rr+pad), min(im.size[1],bb+pad))).save(out_png)
    print("FINAL merged poster:", os.path.getsize(out_png), "B", Image.open(out_png).size)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--html", default=os.path.dirname(os.path.abspath(__file__)))
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--threedmol", default="", help="local 3Dmol-min.js (else fetched from cdnjs)")
    a = ap.parse_args()
    build(a.html, a.outdir, a.threedmol)
