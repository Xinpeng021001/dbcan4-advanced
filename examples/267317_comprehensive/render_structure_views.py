#!/usr/bin/env python3
"""
render_structure_views.py — render a protein structure in MULTIPLE colouring schemes, including a
reusable "colour by domain" mode driven by Pfam/CAZyme domain boundaries.

Colouring modes (‑‑modes, default all four):
    plddt     per-residue ESMFold confidence   (3Dmol 'b'-factor gradient, blue=confident→red=low)
    spectrum  N→C sequence position            (rainbow, blue N-terminus → red C-terminus)
    sstruc    secondary structure              (helix / sheet / coil, 3Dmol ssJmol scheme)
    domain    Pfam / CAZyme domains            (each domain a distinct colour, linker grey) <-- the
                                                requested capability; see domain_cartoon_style_js()

The domain spec is a JSON list of {name, family, start, end, color}; pass it with ‑‑domains FILE
or rely on the built-in 267317 spec. `domain_cartoon_style_js()` is the reusable primitive: hand it
your hmmscan hits (each mapped to a CAZy family + a colour) and it emits the 3Dmol setStyle() calls
that light up each domain on the fold.

Pipeline: 3Dmol.js renders each mode's cartoon offscreen in headless Chrome and reads the framebuffer
back to a PNG via pngURI() (preserveDrawingBuffer + synchronous readback). matplotlib then composes a
labelled multi-panel figure and a large standalone domain figure with a colour→domain legend.

Inputs:
    --pdb PATH        ESMFold (or any) PDB. Default: the 267317 PDB embedded in comprehensive HTML if
                      --html-with-pdb is given, else --pdb is required.
    --domains FILE    JSON list of domain dicts (see above). Default: built-in 267317 spec.
    --threedmol PATH  local 3Dmol-min.js (else fetched from cdnjs; met can reach it).
Outputs (‑‑outdir, default '.'):
    structure_views_<id>.png            2x(N/2) labelled panel of all modes
    structure_domain_<id>.png           large domain-coloured figure with legend
    view_<id>_<mode>.png                each mode's raw cartoon
Requires: google-chrome-stable/-chrome (SwiftShader WebGL) + Pillow + matplotlib.

Usage:
    python render_structure_views.py --pdb 267317.pdb --outdir out
    python render_structure_views.py --pdb X.pdb --domains X_domains.json --modes plddt,domain
"""
import argparse, base64, json, os, re, shutil, subprocess, sys, urllib.request

CDN = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"

# --- built-in demo spec: 267317's two Pfam domains mapped to CAZy families + the signal peptide ---
DEMO_DOMAINS_267317 = [
    {"name": "Signal peptide", "family": "SP (DeepTMHMM)", "start": 1,   "end": 18,  "color": "#f59e0b"},
    {"name": "Glyco_hydro_28", "family": "GH28 (PF00295)", "start": 62,  "end": 412, "color": "#3b82f6"},
    {"name": "Bac_rhamnosid6H","family": "GH78 (PF17389)", "start": 667, "end": 886, "color": "#10b981"},
]
BASE_DOMAIN_COLOR = "#9ca3af"   # residues in no annotated domain

MODE_LABELS = {
    "plddt":    "pLDDT confidence  (blue = confident \u2192 red = low)",
    "spectrum": "Sequence position  (N-terminus blue \u2192 C-terminus red)",
    "sstruc":   "Secondary structure  (helix / sheet / coil)",
    "domain":   "Pfam / CAZyme domains",
}


def domain_cartoon_style_js(domains, base_color=BASE_DOMAIN_COLOR):
    """Reusable primitive: emit 3Dmol.js setStyle() statements colouring a cartoon by sequence domains.

    Each `domains` entry needs 1-based integer 'start'/'end' residue numbers and a hex 'color'.
    Residues outside every domain get `base_color`. Pass Pfam hmmscan hits (each mapped to a CAZy
    family + colour) to colour a structure by its CAZyme/Pfam domain architecture.
    """
    stmts = ["v.setStyle({},{cartoon:{color:'%s'}});" % base_color]
    for d in domains:
        stmts.append("v.setStyle({resi:'%d-%d'},{cartoon:{color:'%s'}});"
                     % (int(d["start"]), int(d["end"]), d["color"]))
    return "".join(stmts)


def _mode_style_js(mode, domains):
    if mode == "plddt":
        return "v.setStyle({},{cartoon:{colorscheme:{prop:'b',gradient:'roygb',min:50,max:90}}});"
    if mode == "spectrum":
        return "v.setStyle({},{cartoon:{color:'spectrum'}});"
    if mode == "sstruc":
        return "v.setStyle({},{cartoon:{colorscheme:'ssJmol'}});"
    if mode == "domain":
        return domain_cartoon_style_js(domains)
    raise ValueError("unknown mode: " + mode)


def chrome():
    exe = shutil.which("google-chrome-stable") or shutil.which("google-chrome") or shutil.which("chromium")
    if not exe:
        sys.exit("no chrome/chromium on PATH")
    return exe


def render_modes(pdb, modes, domains, js, chrome_exe, outdir, wh=(1200, 1000)):
    """Render each mode's cartoon to a raw PNG via one headless-Chrome pass (pngURI readback)."""
    w, h = wh
    viewers, outs = [], []
    for i, m in enumerate(modes):
        viewers.append(
            '<div id="v%d" style="width:%dpx;height:%dpx;position:relative"></div>'
            '<div id="out%d">PENDING</div>' % (i, w, h, i))
        style = _mode_style_js(m, domains)
        outs.append(
            "try{var v=$3Dmol.createViewer(document.getElementById('v%d'),"
            "{backgroundColor:'#0d1b2a',preserveDrawingBuffer:true,antialias:true});"
            "v.addModel(pdb,'pdb');%s v.zoomTo();v.render();"
            "document.getElementById('out%d').textContent=v.pngURI();}"
            "catch(e){document.getElementById('out%d').textContent='ERR:'+e;}" % (i, style, i, i))
    page = ('<!doctype html><html><head><meta charset="utf-8"><script>' + js + '</script></head>'
            '<body style="margin:0">' + "".join(viewers) +
            '<script id="pdb" type="text/plain">' + pdb + '</script>'
            '<script>function go(){var pdb=document.getElementById("pdb").textContent;' + "".join(outs) +
            '}if(document.readyState!=="loading")go();else document.addEventListener("DOMContentLoaded",go);'
            '</script></body></html>')
    tmp = os.path.abspath(os.path.join(outdir, "_views.html"))
    open(tmp, "w", encoding="utf-8").write(page)
    cmd = [chrome_exe, "--headless=new", "--no-sandbox", "--in-process-gpu", "--use-gl=angle",
           "--use-angle=swiftshader-webgl", "--enable-unsafe-swiftshader", "--ignore-gpu-blocklist",
           "--virtual-time-budget=30000", "--run-all-compositor-stages-before-draw",
           "--dump-dom", "file://" + tmp]
    dom = subprocess.run(cmd, capture_output=True, text=True, timeout=240).stdout
    from PIL import Image, ImageChops
    paths = {}
    for i, m in enumerate(modes):
        mm = re.search(r'id="out%d">(data:image/png;base64,[A-Za-z0-9+/=]+)<' % i, dom)
        if not mm:
            m2 = re.search(r'id="out%d">([^<]*)<' % i, dom)
            sys.exit("mode %s render failed: %s" % (m, (m2.group(1)[:160] if m2 else "?")))
        p = os.path.join(outdir, "view_%s.png" % m)
        open(p, "wb").write(base64.b64decode(mm.group(1).split(",", 1)[1]))
        im = Image.open(p).convert("RGB")
        bbox = ImageChops.difference(im, Image.new("RGB", im.size, (13, 27, 42))).getbbox()
        if bbox and (bbox[2] - bbox[0]) < im.size[0]:
            pad = 40; l, t, rr, bb = bbox
            im.crop((max(0, l - pad), max(0, t - pad), min(im.size[0], rr + pad), min(im.size[1], bb + pad))).save(p)
        paths[m] = p
        print("rendered mode %-9s -> %s (%d B)" % (m, p, os.path.getsize(p)))
    return paths


def compose(paths, modes, domains, outdir, sid):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch
    from PIL import Image

    # (1) multi-panel figure of all modes
    n = len(modes); ncol = 2 if n > 1 else 1; nrow = (n + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(6.2 * ncol, 5.4 * nrow), facecolor="white")
    axes = (axes.ravel() if hasattr(axes, "ravel") else [axes])
    for ax, m in zip(axes, modes):
        ax.imshow(Image.open(paths[m])); ax.axis("off")
        ax.set_title(MODE_LABELS.get(m, m), fontsize=12, fontweight="bold", pad=8)
        if m == "domain":
            handles = [Patch(facecolor=d["color"], edgecolor="none", label="%s \u00b7 %d\u2013%d" %
                             (d["family"], d["start"], d["end"])) for d in domains]
            handles.append(Patch(facecolor=BASE_DOMAIN_COLOR, edgecolor="none", label="inter-domain / linker"))
            ax.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.16),
                      ncol=2, fontsize=8.5, frameon=False)
    for ax in axes[len(modes):]:
        ax.axis("off")
    fig.suptitle("267317 \u00b7 ESMFold structure \u00b7 colouring schemes", fontsize=15, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    panel = os.path.join(outdir, "structure_views_%s.png" % sid)
    fig.savefig(panel, dpi=150, bbox_inches="tight", facecolor="white"); plt.close(fig)
    print("panel figure ->", panel)

    # (2) large standalone domain figure with legend
    if "domain" in paths:
        fig2, ax2 = plt.subplots(figsize=(9, 9), facecolor="white")
        ax2.imshow(Image.open(paths["domain"])); ax2.axis("off")
        ax2.set_title("267317 \u00b7 ESMFold model coloured by Pfam / CAZyme domain",
                      fontsize=14, fontweight="bold", pad=10)
        handles = [Patch(facecolor=d["color"], edgecolor="none",
                         label="%s  \u00b7  %s  \u00b7  res %d\u2013%d" % (d["name"], d["family"], d["start"], d["end"]))
                   for d in domains]
        handles.append(Patch(facecolor=BASE_DOMAIN_COLOR, edgecolor="none", label="inter-domain / linker"))
        ax2.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.10),
                   ncol=1, fontsize=10, frameon=False)
        dom_fig = os.path.join(outdir, "structure_domain_%s.png" % sid)
        fig2.savefig(dom_fig, dpi=150, bbox_inches="tight", facecolor="white"); plt.close(fig2)
        print("domain figure ->", dom_fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--pdb", required=True)
    ap.add_argument("--domains", default="")
    ap.add_argument("--modes", default="plddt,spectrum,sstruc,domain")
    ap.add_argument("--threedmol", default="")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--id", default="267317")
    a = ap.parse_args()
    os.makedirs(a.outdir, exist_ok=True)
    pdb = open(a.pdb, encoding="utf-8", errors="replace").read()
    n_ca = sum(1 for l in pdb.splitlines() if l.startswith("ATOM") and l[12:16].strip() == "CA")
    print("PDB CA atoms:", n_ca)
    domains = json.load(open(a.domains)) if a.domains else DEMO_DOMAINS_267317
    js = (open(a.threedmol, encoding="utf-8", errors="replace").read() if a.threedmol
          else urllib.request.urlopen(CDN, timeout=60).read().decode("utf-8", "replace"))
    modes = [m.strip() for m in a.modes.split(",") if m.strip()]
    paths = render_modes(pdb, modes, domains, js, chrome(), a.outdir)
    compose(paths, modes, domains, a.outdir, a.id)


if __name__ == "__main__":
    main()
