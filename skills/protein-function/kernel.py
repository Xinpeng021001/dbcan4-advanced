"""kernel.py for the protein-function skill.

Helpers that assemble REAL functional evidence for a protein given its CAZy
family. Nothing here invents annotation values: GO terms and domains come only
from evidence dicts you pass in (fetched from the protein-annotation MCP
connector in the repl tool); EC/substrate come only from a dbCAN reference you
supply; physicochemistry/sequence features are computed from the sequence.

Top level defines only functions + literal constants (sidecar rules).
Third-party imports (Bio) are deferred into function bodies.
"""

# Kyte-Doolittle hydropathy scale (Kyte & Doolittle, 1982) — literal constant.
KD_SCALE = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
}
STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
# GO single-letter category codes -> aspect names.
GO_ASPECT = {"F": "molecular_function", "P": "biological_process", "C": "cellular_component"}


def clean_sequence(seq):
    """Uppercase, strip whitespace/gaps/stops. Return (standard_only_seq, sorted_nonstandard_chars)."""
    s = "".join(str(seq).split()).upper().replace("*", "").replace("-", "")
    nonstd = sorted(set(c for c in s if c not in STANDARD_AA))
    std = "".join(c for c in s if c in STANDARD_AA)
    return std, nonstd


def physicochem(seq):
    """Biopython ProteinAnalysis physicochemistry of a protein sequence.

    Returns MW (Da), theoretical pI, aromaticity, instability index (+stable/unstable
    class at the 40 cutoff), GRAVY, and per-residue aa composition (%). Non-standard
    residues (X/B/Z/U/O) are dropped before analysis and reported separately.
    """
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    std, nonstd = clean_sequence(seq)
    pa = ProteinAnalysis(std)
    ii = pa.instability_index()
    # Composition from raw counts -> percent (version-robust across Biopython releases;
    # the amino_acids_percent property changed name/scale between versions).
    counts = pa.count_amino_acids()
    total = sum(counts.values()) or 1
    comp = {k: round(100.0 * v / total, 2) for k, v in sorted(counts.items())}
    return {
        "length": len(std),
        "nonstandard_residues_removed": nonstd,
        "molecular_weight": round(pa.molecular_weight(), 2),
        "theoretical_pI": round(pa.isoelectric_point(), 2),
        "aromaticity": round(pa.aromaticity(), 4),
        "instability_index": round(ii, 2),
        "instability_class": "stable" if ii < 40 else "unstable",
        "gravy": round(pa.gravy(), 4),
        "aa_composition_pct": comp,
    }


def n_glyc_sites(seq):
    """N-linked glycosylation sequons N-X-[S/T], X != P (overlapping matches).

    Returns a list of {position (1-based N), sequon} dicts.
    """
    import re
    s = "".join(str(seq).split()).upper().replace("*", "").replace("-", "")
    return [{"position": m.start() + 1, "sequon": m.group(1)}
            for m in re.finditer(r"(?=(N[^P][ST]))", s)]


def signal_region(seq, window=None, threshold=None):
    """Transparent N-terminal secretory signal-peptide heuristic (NOT SignalP/DeepLoc).

    Scans the first 25 residues for the most hydrophobic window (Kyte-Doolittle);
    a signal peptide is called when a strongly hydrophobic h-region sits near the
    N-terminus. Returns the N-terminal residues, n-region positive charges, the best
    hydrophobic window, and a boolean call. Use this only as a heuristic prior — for
    a committed prediction run SignalP/DeepLoc (see the deeploc skill).
    """
    if window is None:
        window = 9
    if threshold is None:
        threshold = 1.6
    std, _nonstd = clean_sequence(seq)
    nterm = std[:30]
    scan = std[:25]
    best = {"start": None, "mean_kd": None, "peptide": None}
    for i in range(0, max(0, len(scan) - window + 1)):
        w = scan[i:i + window]
        mk = sum(KD_SCALE[c] for c in w) / len(w)
        if best["mean_kd"] is None or mk > best["mean_kd"]:
            best = {"start": i + 1, "mean_kd": round(mk, 3), "peptide": w}
    ncharge = sum(1 for c in std[:5] if c in "KR")
    has_signal = bool(best["mean_kd"] is not None and best["mean_kd"] >= threshold
                      and best["start"] is not None and best["start"] <= 15)
    return {
        "nterminal_1_20": nterm[:20],
        "n_region_pos_charges_1_5": ncharge,
        "h_region_best_window": best,
        "hydrophobicity_threshold": threshold,
        "predicted_signal_peptide": has_signal,
        "method": "Kyte-Doolittle N-terminal hydropathy heuristic (not SignalP/DeepLoc)",
    }


def load_famref(famref):
    """Accept a dbCAN fam->EC/substrate reference as a dict, a JSON path, or None. Return dict ({} if unavailable)."""
    import json
    import os
    if famref is None:
        return {}
    if isinstance(famref, dict):
        return famref
    if isinstance(famref, str) and os.path.exists(famref):
        with open(famref) as fh:
            return json.load(fh)
    return {}


def map_family(cazy_family, famref):
    """Map a CAZy family (e.g. 'GH78', or subfamily 'GH5_4') to EC + substrate via a dbCAN reference.

    famref: dict or JSON path keyed by family. Tries the exact label, then the base
    family (subfamily suffix stripped). Returns {family, ec, substrate, substrate_high,
    name, alternatives, mapped, note}. Returns mapped=False (never a guessed EC) when
    the family is absent from the reference.
    """
    ref = load_famref(famref)
    fam = str(cazy_family).strip()
    entry = ref.get(fam)
    if entry is None and "_" in fam:
        entry = ref.get(fam.split("_")[0])
    if not entry:
        return {"family": cazy_family, "ec": None, "substrate": None, "substrate_high": None,
                "name": None, "alternatives": [], "mapped": False,
                "note": "family not found in provided reference"}
    return {
        "family": cazy_family,
        "ec": entry.get("ec") or None,
        "substrate": entry.get("substrate") or entry.get("substrate_high") or None,
        "substrate_high": entry.get("substrate_high") or None,
        "name": entry.get("name") or None,
        "alternatives": entry.get("all", []),
        "mapped": True,
        "note": None,
    }


def collect_go(entries):
    """Flatten + dedup GO terms from InterPro/Pfam entry dicts.

    Each entry may carry a 'go_terms' or 'go' list of GO dicts (as returned by the
    protein-annotation connector). Returns a list of {identifier, name, category}
    sorted by GO id. Category is normalised to the aspect name
    (molecular_function/biological_process/cellular_component).
    """
    seen = {}
    for e in entries or []:
        if not isinstance(e, dict):
            continue
        for key in ("go_terms", "go"):
            for g in (e.get(key) or []):
                if not isinstance(g, dict):
                    continue
                ident = g.get("identifier") or g.get("id") or g.get("accession")
                if not ident:
                    continue
                cat = g.get("category")
                if isinstance(cat, dict):
                    cat = cat.get("name") or GO_ASPECT.get(cat.get("code"), cat.get("code"))
                elif isinstance(cat, str) and cat in GO_ASPECT:
                    cat = GO_ASPECT[cat]
                seen[ident] = {"identifier": ident, "name": g.get("name"), "category": cat}
    return sorted(seen.values(), key=lambda x: x["identifier"])


def collect_domains(entries):
    """Normalise InterPro/Pfam entry dicts to {accession, name, type, source_database}, deduped by accession."""
    out, seen = [], set()
    for e in entries or []:
        if not isinstance(e, dict):
            continue
        acc = e.get("accession") or e.get("acc")
        if not acc or acc in seen:
            continue
        name = e.get("name")
        if isinstance(name, dict):
            name = name.get("name") or name.get("short")
        seen.add(acc)
        out.append({
            "accession": acc,
            "name": name,
            "type": e.get("type"),
            "source_database": e.get("source_database") or e.get("database"),
        })
    return out


def group_go_by_aspect(go):
    """Split a normalised GO list into molecular_function / biological_process / cellular_component buckets."""
    buckets = {"molecular_function": [], "biological_process": [], "cellular_component": []}
    for g in go or []:
        cat = (g.get("category") or "").lower()
        cat = GO_ASPECT.get(cat.upper(), cat)
        if cat in buckets:
            buckets[cat].append(g)
    return buckets


def annotate_function(seq, cazy_family, interpro_accessions=None, famref=None,
                      go_terms=None, domains=None, interpro_entries=None, protein_id=None):
    """Assemble real functional evidence for a protein given its CAZy family.

    Args:
        seq: protein sequence (str).
        cazy_family: CAZy family label, e.g. 'GH78'.
        interpro_accessions: optional list of accessions to record verbatim.
        famref: dbCAN fam->EC/substrate reference (dict or JSON path).
        go_terms: pre-fetched GO dicts (from the protein-annotation connector, repl tool).
        domains: pre-fetched domain dicts (e.g. get_domain_architecture entries).
        interpro_entries: raw InterPro/Pfam entry dicts carrying go_terms — GO and
            domains are extracted from them when go_terms/domains are not given directly.
        protein_id: optional identifier to stamp on the record.

    Returns a structured dict combining EC/substrate (from famref), GO (grouped by
    aspect), domains, physicochemistry, N-glyc sequons, and the signal-region heuristic.
    GO/domains are echoed ONLY from evidence you pass in — this helper never invents them.
    """
    entries = list(interpro_entries) if interpro_entries else []
    raw_go = list(go_terms) if go_terms else []
    if not raw_go and entries:
        go = collect_go(entries)
    else:
        go = collect_go([{"go_terms": raw_go}]) if raw_go else []
    dom = list(domains) if domains else (collect_domains(entries) if entries else [])
    ec_sub = map_family(cazy_family, famref)
    pchem = physicochem(seq)
    glyc = n_glyc_sites(seq)
    sig = signal_region(seq)
    return {
        "protein_id": protein_id,
        "length_aa": pchem["length"],
        "cazy_family": cazy_family,
        "ec_number": ec_sub["ec"],
        "substrate": ec_sub["substrate"],
        "substrate_high": ec_sub["substrate_high"],
        "activity_name": ec_sub["name"],
        "ec_substrate_alternatives": ec_sub["alternatives"],
        "ec_substrate_mapped": ec_sub["mapped"],
        "interpro_accessions": list(interpro_accessions) if interpro_accessions else [],
        "domains": dom,
        "n_domains": len(dom),
        "go_terms": go,
        "go_by_category": group_go_by_aspect(go),
        "n_go_terms": len(go),
        "physicochemistry": pchem,
        "n_glycosylation_sites": glyc,
        "n_glyc_site_count": len(glyc),
        "signal_region": sig,
    }
