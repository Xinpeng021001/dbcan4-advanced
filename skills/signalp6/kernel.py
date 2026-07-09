"""kernel.py for the signalp6 skill: locate + run + parse the licensed SignalP-6.0 CLI.

Scaffold: wraps `signalp6` only when present on PATH. Never fabricates output.
"""

SIGNALP6_RESULTS = "prediction_results.txt"
SIGNALP6_INSTALL = ("SignalP-6.0 is license-gated (academic download). Request it at "
                    "https://services.healthtech.dtu.dk/services/SignalP-6.0/ , then "
                    "`pip install signalp-6-package/` and `signalp6-register <weights>`. "
                    "It is NOT pip-installable from PyPI without the license file.")


def find_signalp6():
    """Return the path to the signalp6 binary on PATH, or None if not installed."""
    import shutil
    return shutil.which("signalp6")


def run_signalp6(fasta, out_dir=None, organism=None, mode=None, fmt=None):
    """Run the SignalP-6.0 CLI if present. Returns a status dict; never raises for the common cases.

    status is one of: "not_installed" (binary absent - see result["install"]),
    "ok" (returncode 0), "failed" (nonzero exit or exception).
    """
    import os, subprocess, shutil
    if out_dir is None:
        out_dir = "signalp6_out"
    if organism is None:
        organism = "eukarya"
    if mode is None:
        mode = "fast"
    if fmt is None:
        fmt = "txt"
    result = {"tool": "signalp6", "out_dir": out_dir,
              "fasta": os.path.abspath(fasta) if os.path.exists(fasta) else fasta}
    binp = shutil.which("signalp6")
    if binp is None:
        result["status"] = "not_installed"
        result["error"] = "signalp6 not found on PATH."
        result["install"] = SIGNALP6_INSTALL
        return result
    result["binary"] = binp
    os.makedirs(out_dir, exist_ok=True)
    cmd = [binp, "--fastafile", fasta, "--output_dir", out_dir,
           "--organism", organism, "--mode", mode, "--format", fmt]
    result["command"] = " ".join(cmd)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        result["returncode"] = proc.returncode
        result["stdout_tail"] = (proc.stdout or "")[-2000:]
        result["stderr_tail"] = (proc.stderr or "")[-2000:]
        result["status"] = "ok" if proc.returncode == 0 else "failed"
    except Exception as e:
        import traceback
        result["status"] = "failed"
        result["error"] = repr(e)
        result["traceback"] = traceback.format_exc()
    return result


def parse_signalp6(out_dir):
    """Parse SignalP-6.0 prediction_results.txt (TSV) into a dict keyed by protein id.

    See SKILL.md 'Output schema'. Returns status="no_output_found" if the TSV is absent.
    """
    import os, glob, re
    pred_norm = {"OTHER": "NO_SP", "SP": "SP", "LIPO": "LIPO", "TAT": "TAT",
                 "TATLIPO": "TATLIPO", "PILIN": "PILIN"}
    result = {"status": "ok", "n_proteins": 0, "proteins": {}, "source_files": {}}
    tsvs = glob.glob(os.path.join(out_dir, "**", SIGNALP6_RESULTS), recursive=True)
    if not tsvs:
        tsvs = glob.glob(os.path.join(out_dir, "**", "*.txt"), recursive=True)
    result["source_files"]["tsv"] = tsvs
    if not tsvs:
        result["status"] = "no_output_found"
        return result

    header = None
    with open(tsvs[0]) as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                # a header comment carries the column names (contains "Prediction")
                if "Prediction" in line and ("\t" in line):
                    header = line.lstrip("#").strip().split("\t")
                continue
            cols = line.split("\t")
            if not cols or not cols[0]:
                continue
            pid = cols[0]
            row = dict(zip(header, cols)) if (header and len(header) == len(cols)) else {}
            ptype = cols[1].strip() if len(cols) > 1 else ""
            # probability columns: any numeric-looking header that is a class prob
            all_probs = {}
            for k, v in row.items():
                kk = k.strip()
                if kk in ("ID", "Prediction", "CS Position"):
                    continue
                try:
                    all_probs[kk] = float(v)
                except Exception:
                    pass
            # winning probability: prob column whose name starts with the predicted type
            sp_prob = None
            for k, v in all_probs.items():
                base = k.split("(")[0].strip().upper()
                if ptype and base == ptype.upper():
                    sp_prob = v
                    break
            if sp_prob is None and ptype.upper() != "OTHER" and all_probs:
                # fall back to the max non-OTHER probability
                cand = {k: v for k, v in all_probs.items() if not k.upper().startswith("OTHER")}
                if cand:
                    sp_prob = max(cand.values())
            # cleavage site from "CS Position" e.g. "CS pos: 23-24. Pr: 0.87"
            cs = row.get("CS Position", "") if row else ""
            if not cs:
                m0 = re.search(r"CS pos:\s*([0-9]+-[0-9]+)", line)
                cs = m0.group(0) if m0 else ""
            cleavage_site = None
            cleavage_prob = None
            m = re.search(r"([0-9]+-[0-9]+)", cs)
            if m:
                cleavage_site = m.group(1)
            mp = re.search(r"Pr:\s*([0-9.]+)", cs)
            if mp:
                try:
                    cleavage_prob = float(mp.group(1))
                except Exception:
                    pass
            result["proteins"][pid] = {
                "prediction": pred_norm.get(ptype.upper(), ptype),
                "prediction_raw": ptype,
                "sp_probability": sp_prob,
                "cleavage_site": cleavage_site,
                "cleavage_prob": cleavage_prob,
                "all_probs": all_probs,
            }
    result["n_proteins"] = len(result["proteins"])
    return result
