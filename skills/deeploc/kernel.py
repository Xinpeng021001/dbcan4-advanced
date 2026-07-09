"""kernel.py for the deeploc skill (subcellular localization).

DeepLoc-2.0 (DTU Health Tech) is license-gated and not pip/conda installable, so
this is a SCAFFOLD: run_deeploc() wraps the CLI when it is present on PATH and
parses its output; when it is absent it returns status 'not_installed' with the
real install instructions rather than any invented localization.

derive_localization_from_evidence() provides a transparent fallback call built
ONLY from evidence you already have (a signal-peptide flag + GO cellular-component
terms). It is explicitly labelled as derived, never presented as a DeepLoc result.

Top level defines only functions + literal constants (sidecar rules). Third-party
imports are deferred into function bodies (none are needed here — stdlib only).
"""

DEEPLOC_INSTALL_DOC = (
    "DeepLoc-2.0 is distributed by DTU Health Tech under an academic license and is "
    "NOT on PyPI/conda. Install: (1) register and download from "
    "https://services.healthtech.dtu.dk/services/DeepLoc-2.0/ (academic use); "
    "(2) `pip install deeploc-2.0.tar.gz` (pulls torch, fair-esm, etc.); "
    "(3) verify with `deeploc2 --help`. First run downloads the ESM-1b/ProtT5 weights. "
    "CLI: `deeploc2 -f input.fasta -o outdir [-m Accurate|Fast]` writes a results_*.csv."
)
# DeepLoc-2.0's 10 subcellular-localization classes.
DEEPLOC_CLASSES = [
    "Cytoplasm", "Nucleus", "Extracellular", "Cell membrane", "Mitochondrion",
    "Plastid", "Endoplasmic reticulum", "Lysosome/Vacuole", "Golgi apparatus", "Peroxisome",
]
# GO cellular-component identifiers -> a coarse localization label, for the derived call.
GO_CC_TO_LOCALIZATION = {
    "GO:0005576": "Extracellular",              # extracellular region
    "GO:0005615": "Extracellular",              # extracellular space
    "GO:0005737": "Cytoplasm",                  # cytoplasm
    "GO:0005634": "Nucleus",                    # nucleus
    "GO:0016020": "Cell membrane",              # membrane
    "GO:0005886": "Cell membrane",              # plasma membrane
    "GO:0005739": "Mitochondrion",              # mitochondrion
    "GO:0005783": "Endoplasmic reticulum",      # ER
    "GO:0005794": "Golgi apparatus",            # Golgi
    "GO:0005764": "Lysosome/Vacuole",           # lysosome
    "GO:0005773": "Lysosome/Vacuole",           # vacuole
    "GO:0005777": "Peroxisome",                 # peroxisome
    "GO:0009986": "Cell membrane",              # cell surface
    "GO:0031012": "Extracellular",              # extracellular matrix
}


def deeploc_available():
    """Return the deeploc2 CLI path if on PATH, else None."""
    import shutil
    return shutil.which("deeploc2") or shutil.which("deeploc")


def run_deeploc(fasta_path, out_dir=None, model=None):
    """Run DeepLoc-2.0 on a FASTA if the CLI is installed; else report not_installed.

    Returns a dict with a truthful 'status':
      - 'not_installed' + install_doc when the CLI is absent (NO invented output);
      - 'failed' + real stderr when the CLI ran but errored;
      - 'success' + parsed per-protein predictions when it worked.
    Never fabricates localization values.
    """
    import os
    import glob
    import subprocess
    if out_dir is None:
        out_dir = "deeploc_out"
    if model is None:
        model = "Accurate"
    exe = deeploc_available()
    if not exe:
        return {"tool": "DeepLoc-2.0", "status": "not_installed",
                "error": "deeploc2 CLI not found on PATH",
                "install_doc": DEEPLOC_INSTALL_DOC, "predictions": []}
    os.makedirs(out_dir, exist_ok=True)
    try:
        proc = subprocess.run([exe, "-f", fasta_path, "-o", out_dir, "-m", model],
                              capture_output=True, text=True)
    except Exception as e:
        return {"tool": "DeepLoc-2.0", "status": "failed", "error": repr(e),
                "install_doc": DEEPLOC_INSTALL_DOC, "predictions": []}
    if proc.returncode != 0:
        return {"tool": "DeepLoc-2.0", "status": "failed",
                "error": (proc.stderr or proc.stdout or "").strip()[:2000],
                "predictions": []}
    csvs = sorted(glob.glob(os.path.join(out_dir, "results*.csv")) +
                  glob.glob(os.path.join(out_dir, "*.csv")))
    return {"tool": "DeepLoc-2.0", "status": "success", "output_dir": out_dir,
            "result_files": csvs, "predictions": parse_deeploc_csv(csvs[0]) if csvs else [],
            "stdout_tail": (proc.stdout or "").strip()[-500:]}


def parse_deeploc_csv(csv_path):
    """Parse a DeepLoc-2.0 results CSV into a list of per-protein dicts (stdlib csv)."""
    import csv
    rows = []
    with open(csv_path, newline="") as fh:
        for row in csv.DictReader(fh):
            rows.append(dict(row))
    return rows


def derive_localization_from_evidence(signal_region=None, go_terms=None, cazy_family=None):
    """Derive a localization call from evidence you ALREADY have — NOT a DeepLoc prediction.

    Combines an N-terminal signal-peptide flag (secretory pathway prior) with any GO
    cellular-component terms. Returns the call, a confidence label, and the exact
    evidence + reasoning used, so the provenance is explicit ('derived from
    signal-peptide + GO-CC', never attributed to DeepLoc).

    Args:
        signal_region: dict from protein-function's signal_region() (uses
            'predicted_signal_peptide' and 'h_region_best_window').
        go_terms: normalised GO dicts ({identifier,name,category}); only
            cellular_component terms are used.
        cazy_family: optional family label, recorded as context only (not decisive).
    """
    evidence, cc_locs = [], []
    go_cc = []
    for g in (go_terms or []):
        cat = (g.get("category") or "")
        if isinstance(cat, str) and cat.lower() == "cellular_component":
            go_cc.append(g)
            loc = GO_CC_TO_LOCALIZATION.get(g.get("identifier"))
            if loc:
                cc_locs.append(loc)
                evidence.append("GO-CC %s (%s) -> %s" % (g.get("identifier"), g.get("name"), loc))

    has_signal = bool(signal_region and signal_region.get("predicted_signal_peptide"))
    if has_signal:
        hw = (signal_region or {}).get("h_region_best_window", {})
        evidence.append("N-terminal signal peptide predicted (h-region %s, mean KD %s) -> secretory pathway"
                        % (hw.get("peptide"), hw.get("mean_kd")))

    # Decision: GO-CC evidence wins when present; otherwise a signal peptide implies
    # secretory routing (Extracellular for a fungal CAZyme unless CC says otherwise).
    if cc_locs:
        call = max(set(cc_locs), key=cc_locs.count)
        confidence = "medium" if has_signal or len(cc_locs) > 1 else "low"
        basis = "GO cellular-component terms" + (" + signal peptide" if has_signal else "")
    elif has_signal:
        call = "Extracellular (secretory pathway)"
        confidence = "medium"
        basis = "N-terminal signal peptide (no GO-CC available)"
    else:
        call = "Undetermined"
        confidence = "none"
        basis = "no signal peptide and no GO cellular-component evidence"

    return {
        "method": "derived from signal-peptide + GO-CC (NOT DeepLoc)",
        "localization_call": call,
        "confidence": confidence,
        "basis": basis,
        "signal_peptide": has_signal,
        "go_cellular_component_terms": go_cc,
        "evidence": evidence,
        "cazy_family_context": cazy_family,
    }
