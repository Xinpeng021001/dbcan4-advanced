#!/usr/bin/env python3
"""Export the BioForge advanced-vs-baseline demo as a self-contained, offline
HTML bundle — no running server needed.

For each of the demo genes plus a browse landing page it renders the real Jinja
templates through the FastAPI app, then rewrites the output so the bundle opens
straight from disk:

  * the shared stylesheet is inlined into every page,
  * internal /genes/<id> links become gene_<key>.html,
  * the 3D structure viewer reads its PDB from an embedded block (so it works
    from file://, where fetch() is blocked), and the PDB files are also copied
    in so the download link works,
  * 3Dmol.js is inlined when reachable, else left as a CDN <script> (works
    online).

Usage:
  DATABASE_URL=sqlite:///demo.db BIOFORGE_TRACKS_DIR=.../tracks \\
    python demo_export.py --structures-dir .../web_static/structures \\
      --out demo_export
"""
from __future__ import annotations

import argparse
import re
import shutil
from pathlib import Path

from sqlalchemy import create_engine, select
from sqlalchemy.orm import Session
from starlette.testclient import TestClient

from bioforge.api.main import create_app
from bioforge.models import Gene

CDN_3DMOL = "https://cdnjs.cloudflare.com/ajax/libs/3Dmol/2.1.0/3Dmol-min.js"

# viewer init: swap the fetch() for reading an embedded <script> block so the
# structure renders from file:// (browsers block fetch of local files).
_FETCH_SRC = ("var url = el.dataset.pdb;\n"
              "  var viewer = $3Dmol.createViewer(el, { backgroundColor: 'white' });\n"
              "  fetch(url).then(function (r) { return r.text(); }).then(function (pdb) {")
_FETCH_OFFLINE = ("var emb = document.getElementById('pdb-embed');\n"
                  "  var viewer = $3Dmol.createViewer(el, { backgroundColor: 'white' });\n"
                  "  Promise.resolve(emb ? emb.textContent : '').then(function (pdb) {")


def _inline_css(html: str, css: str) -> str:
    return html.replace('<link rel="stylesheet" href="/static/app.css">',
                        f"<style>\n{css}\n</style>")


def _rewrite_links(html: str, id2key: dict[int, str]) -> str:
    # /genes/<id> -> gene_<key>.html
    def gene_sub(m):
        gid = int(m.group(1))
        return f'href="gene_{id2key.get(gid, gid)}.html"'
    html = re.sub(r'href="/genes/(\d+)"', gene_sub, html)
    # browse / home / root-relative app links -> index (offline landing)
    html = re.sub(r'href="/browse(\?[^"]*)?"', 'href="index.html"', html)
    html = re.sub(r'href="/"(?=[ >])', 'href="index.html"', html)
    # structures served relatively (download link)
    html = html.replace('href="/structures/', 'href="structures/')
    html = html.replace('data-pdb="/structures/', 'data-pdb="structures/')
    # remaining app routes we don't export -> neutralise (keep external http links)
    html = re.sub(r'href="/(samples|compare|search|releases|docs|tracks|browse\.csv|browse\.faa)[^"]*"',
                  'href="#" class="disabled-link"', html)
    return html


def _embed_structure(html: str, client: TestClient) -> str:
    """Insert the gene's PDB text as an embedded block and switch the viewer to
    read it offline. Runs on RAW html (path is '/structures/...')."""
    m = re.search(r'<div id="viewer3d" data-pdb="(/structures/[^"]+)"></div>', html)
    if not m:
        return html
    path = m.group(1)  # e.g. /structures/demo_fungal/891310.pdb
    pdb_txt = client.get(path).text
    block = (f'<div id="viewer3d" data-pdb="{path}"></div>\n'
             f'<script type="text/plain" id="pdb-embed">{pdb_txt}</script>')
    html = html.replace(m.group(0), block, 1)
    if _FETCH_SRC not in html:
        raise RuntimeError("viewer init script did not match _FETCH_SRC — "
                           "template changed; update demo_export.py")
    html = html.replace(_FETCH_SRC, _FETCH_OFFLINE)
    return html


def export_demo(out_dir: str | Path, structures_dir: str | Path) -> dict:
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    app = create_app()
    client = TestClient(app)

    css = client.get("/static/app.css").text
    # best-effort inline of 3Dmol.js (works offline if we can fetch it)
    dmol_tag = f'<script src="{CDN_3DMOL}"></script>'
    dmol_inline = dmol_tag
    try:
        import urllib.request
        js = urllib.request.urlopen(CDN_3DMOL, timeout=8).read().decode("utf-8")
        (out / "3Dmol-min.js").write_text(js)
        dmol_inline = '<script src="3Dmol-min.js"></script>'
    except Exception:
        pass  # leave CDN reference; viewer works when opened online

    eng = create_engine(__import__("os").environ["DATABASE_URL"])
    with Session(eng) as s:
        genes = list(s.scalars(select(Gene).order_by(Gene.gene_key)))
        id2key = {g.id: g.gene_key for g in genes}
        gene_ids = [(g.id, g.gene_key) for g in genes]

    # copy served structures into the bundle (for the download links)
    sdir = Path(structures_dir)
    if sdir.exists():
        shutil.copytree(sdir, out / "structures", dirs_exist_ok=True)

    def finish(html: str) -> str:
        html = _inline_css(html, css)
        html = html.replace(dmol_tag, dmol_inline)
        html = _rewrite_links(html, id2key)
        return html

    written = []
    # landing page = browse filtered to advanced-only (the headline view)
    idx = finish(client.get("/browse?advanced_only=1").text)
    idx = idx.replace("<body>", '<body>\n<div class="export-banner">Static export — '
                      'BioForge advanced-vs-baseline fungal-CAZyme demo · '
                      'advanced-only genes shown first</div>')
    (out / "index.html").write_text(idx)
    written.append("index.html")
    # full browse too
    (out / "browse_all.html").write_text(finish(client.get("/browse?sample=demo_fungal").text))
    written.append("browse_all.html")

    for gid, key in gene_ids:
        html = client.get(f"/genes/{gid}").text
        html = _embed_structure(html, client)
        (out / f"gene_{key}.html").write_text(finish(html))
        written.append(f"gene_{key}.html")

    return {"out_dir": str(out), "files": written, "n_genes": len(gene_ids),
            "dmol_inlined": dmol_inline != dmol_tag}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="demo_export")
    ap.add_argument("--structures-dir", required=True)
    a = ap.parse_args()
    rep = export_demo(a.out, a.structures_dir)
    print(f"[export] {len(rep['files'])} pages -> {rep['out_dir']} "
          f"(3Dmol inlined: {rep['dmol_inlined']})")


if __name__ == "__main__":
    main()
