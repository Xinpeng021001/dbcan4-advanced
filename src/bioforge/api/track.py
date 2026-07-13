"""Server-rendered genomic context track (pure SVG, no JS, no external deps).

This is the always-available genome viewer: given a set of features on one
contig and a base-pair region, it draws strand-aware gene arrows positioned by
coordinate, colour-coded by feature type / CAZyme content, on greedily packed
lanes so overlapping features stay readable. Each arrow is a link to its gene
detail page and carries a hover ``<title>``.

Unlike the embedded JBrowse view (which needs bgzip/tabix track files from pysam
and loads its React app from a CDN at view time), this renders inline from the
DB alone, works offline, and themes with the page via CSS classes in app.css.
"""
from __future__ import annotations

from html import escape

# Logical drawing space; the <svg> scales to its container via viewBox + width:100%.
_W = 1000
_PAD_X = 14
_ROW_H = 26      # vertical pitch between lanes
_GENE_H = 15     # arrow body height
_TOP = 48        # headroom for the contig header + focal-gene label rows
_BAND_H = 22     # extra headroom + height for a CGC cluster band, when shown
_AXIS_H = 26     # baseline + coordinate labels below the lanes
_MIN_W = 3.0     # never draw a gene narrower than this (logical units)
_HEAD = 7.0      # arrowhead length (logical units)


def _fmt_bp(n: int) -> str:
    if n >= 1000:
        return f"{n / 1000:.1f} kb"
    return f"{n} bp"


def _gene_class(gene) -> str:
    ft = (gene.feature_type or "").upper()
    if getattr(gene, "cazymes", None):
        return "cazyme"
    if ft in ("CDS", "GENE"):
        return "cds"
    return "rna"


def _pack_lanes(spans: list[tuple[float, float]], gap: float = 8.0) -> list[int]:
    """Greedy interval packing → lane index per span (input order preserved)."""
    order = sorted(range(len(spans)), key=lambda i: spans[i][0])
    lanes: list[float] = []  # rightmost occupied x per lane
    lane_of = [0] * len(spans)
    for i in order:
        x1, x2 = spans[i]
        placed = False
        for li in range(len(lanes)):
            if x1 > lanes[li] + gap:
                lanes[li] = x2
                lane_of[i] = li
                placed = True
                break
        if not placed:
            lane_of[i] = len(lanes)
            lanes.append(x2)
    return lane_of


def render_track(
    genes,
    contig: str,
    region_start: int,
    region_end: int,
    focal_id: int | None = None,
    clusters=None,
) -> str:
    """Return an SVG string for one contig region. ``genes`` are ORM Gene rows.

    ``clusters`` (optional) are CazymeCluster rows on this contig; each is drawn as
    a spanning band above the gene lanes, tying the CGC to its member genes.
    """
    span = max(region_end - region_start, 1)
    inner = _W - 2 * _PAD_X
    clusters = list(clusters or [])
    band_h = _BAND_H if clusters else 0
    top = _TOP + band_h

    def x_of(bp: int) -> float:
        return _PAD_X + (bp - region_start) / span * inner

    spans = [(x_of(g.start), x_of(g.end)) for g in genes]
    lane_of = _pack_lanes(spans) if genes else []
    n_lanes = (max(lane_of) + 1) if lane_of else 1
    lanes_h = n_lanes * _ROW_H
    height = top + lanes_h + _AXIS_H
    axis_y = top + lanes_h + 6

    parts: list[str] = [
        f'<svg class="track" viewBox="0 0 {_W} {height}" role="img" '
        f'preserveAspectRatio="xMidYMid meet" '
        f'aria-label="Genomic context track for {escape(contig)}">',
        f'<text class="contig-label" x="{_PAD_X}" y="18">{escape(contig)}</text>',
        f'<text class="tick-label" x="{_W - _PAD_X}" y="18" text-anchor="end">'
        f'{_fmt_bp(span)} shown</text>',
    ]

    # Cluster bands (CGC): a rounded rect spanning the cluster, above the genes.
    for cl in clusters:
        cx1, cx2 = x_of(cl.start), x_of(cl.end)
        if cx2 - cx1 < 2:
            cx2 = cx1 + 2
        by = _TOP - 4
        label = cl.cluster_key + (f" · {cl.substrate}" if cl.substrate else "")
        title = escape(
            f"{cl.cluster_key} · {cl.contig}:{cl.start}-{cl.end}"
            f"{' · ' + cl.composition if cl.composition else ''}"
            f"{' · substrate ' + cl.substrate if cl.substrate else ''}"
        )
        parts.append(
            f'<rect class="cluster-band" x="{cx1:.1f}" y="{by}" '
            f'width="{cx2 - cx1:.1f}" height="{_BAND_H - 6}" rx="4">'
            f'<title>{title}</title></rect>'
        )
        if cx2 - cx1 > 60:
            parts.append(
                f'<text class="cluster-label" x="{(cx1 + cx2) / 2:.1f}" '
                f'y="{by + _BAND_H - 12:.1f}" text-anchor="middle">{escape(label)}</text>'
            )

    for gene, (x1, x2), lane in zip(genes, spans, lane_of):
        if x2 - x1 < _MIN_W:
            mid = (x1 + x2) / 2
            x1, x2 = mid - _MIN_W / 2, mid + _MIN_W / 2
        y = top + lane * _ROW_H
        yc = y + _GENE_H / 2
        yb = y + _GENE_H
        head = min(_HEAD, (x2 - x1) / 2)
        if gene.strand == "-":
            pts = f"{x2:.1f},{y} {x1 + head:.1f},{y} {x1:.1f},{yc:.1f} " \
                  f"{x1 + head:.1f},{yb} {x2:.1f},{yb}"
        elif gene.strand == "+":
            pts = f"{x1:.1f},{y} {x2 - head:.1f},{y} {x2:.1f},{yc:.1f} " \
                  f"{x2 - head:.1f},{yb} {x1:.1f},{yb}"
        else:
            pts = f"{x1:.1f},{y} {x2:.1f},{y} {x2:.1f},{yb} {x1:.1f},{yb}"

        cls = _gene_class(gene)
        focal = " focal" if focal_id is not None and gene.id == focal_id else ""
        fam = ""
        if getattr(gene, "cazymes", None):
            fam = " · " + ",".join(sorted({c.cazy_family for c in gene.cazymes}))
        title = escape(
            f"{gene.gene_key} · {gene.feature_type} · "
            f"{gene.contig}:{gene.start}-{gene.end}({gene.strand or '.'})"
            f"{' · ' + gene.product if gene.product else ''}{fam}"
        )
        parts.append(
            f'<a href="/genes/{gene.id}">'
            f'<polygon class="gene {cls}{focal}" points="{pts}">'
            f'<title>{title}</title></polygon></a>'
        )
        # Label the focal gene above its lane; label wide neighbours inside.
        if focal:
            parts.append(
                f'<text class="gene-label focal-label" x="{(x1 + x2) / 2:.1f}" '
                f'y="{y - 6}" text-anchor="middle">{escape(gene.gene_key)}</text>'
            )
        elif (x2 - x1) > 90:
            parts.append(
                f'<text class="gene-label" x="{(x1 + x2) / 2:.1f}" '
                f'y="{yc + 4:.1f}" text-anchor="middle">{escape(gene.gene_key)}</text>'
            )

    # Baseline axis with start/end coordinates.
    parts.append(
        f'<line class="axis" x1="{_PAD_X}" y1="{axis_y}" '
        f'x2="{_W - _PAD_X}" y2="{axis_y}"/>'
    )
    parts.append(
        f'<text class="tick-label" x="{_PAD_X}" y="{axis_y + 15}">'
        f'{region_start:,}</text>'
    )
    parts.append(
        f'<text class="tick-label" x="{_W - _PAD_X}" y="{axis_y + 15}" '
        f'text-anchor="end">{region_end:,}</text>'
    )
    parts.append("</svg>")
    return "".join(parts)
