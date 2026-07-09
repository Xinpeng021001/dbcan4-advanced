"""kernel.py for the deeptmhmm skill: run + parse DeepTMHMM (DTU/DeepTMHMM via pybiolib)."""

DEEPTMHMM_APP = "DTU/DeepTMHMM"


def run_deeptmhmm(fasta, workdir=None):
    """Run DeepTMHMM on a FASTA via pybiolib (BioLib cloud) and save outputs to workdir.

    Returns a status dict; never raises for the common failure modes. Check
    result["status"] in {"ok","not_installed","failed"} before calling parse_deeptmhmm.
    """
    import os, traceback
    if workdir is None:
        workdir = "deeptmhmm_out"
    os.makedirs(workdir, exist_ok=True)
    result = {"status": "unknown", "out_dir": workdir, "app": DEEPTMHMM_APP,
              "fasta": os.path.abspath(fasta) if os.path.exists(fasta) else fasta}
    try:
        import biolib
    except Exception as e:
        result["status"] = "not_installed"
        result["error"] = "pybiolib not importable: %r. Install with: pip install pybiolib" % (e,)
        return result
    result["pybiolib_version"] = getattr(biolib, "__version__", "?")
    try:
        app = biolib.load(DEEPTMHMM_APP)
        result["app_loaded"] = True
        job = app.cli(args="--fasta %s" % fasta)
        try:
            result["exit_code"] = job.get_exit_code()
        except Exception as e:
            result["exit_code_err"] = repr(e)
        job.save_files(workdir)
        try:
            result["stdout_tail"] = job.get_stdout().decode("utf-8", "replace")[-2000:]
        except Exception as e:
            result["stdout_err"] = repr(e)
        result["status"] = "ok"
    except Exception as e:
        result["status"] = "failed"
        result["error"] = repr(e)
        result["traceback"] = traceback.format_exc()
    return result


def parse_deeptmhmm(out_dir):
    """Parse DeepTMHMM output (predicted_topologies.3line + TMRs.gff3) into a dict keyed by protein id.

    See the skill's SKILL.md 'Output schema' section for the full return shape.
    """
    import os, glob
    label_names = {"S": "signal_peptide", "I": "inside", "O": "outside",
                   "M": "TMhelix", "B": "beta_strand", "P": "periplasm"}
    pred_norm = {"SP": "SP", "TM": "TM", "SP+TM": "SP+TM", "GLOB": "Globular",
                 "BETA": "TM_beta", "SP+BETA": "SP+TM_beta"}
    result = {"status": "ok", "n_proteins": 0, "proteins": {}, "source_files": {}}
    three = glob.glob(os.path.join(out_dir, "**", "*.3line"), recursive=True)
    gff = glob.glob(os.path.join(out_dir, "**", "*.gff3"), recursive=True)
    result["source_files"]["three_line"] = three
    result["source_files"]["gff3"] = gff
    if not three:
        result["status"] = "no_3line_found"
        return result

    # GFF3: authoritative "Number of predicted TMRs" per id (comment lines)
    gff_tmr_count = {}
    if gff:
        with open(gff[0]) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if "Number of predicted TMRs" in line:
                    body = line.lstrip("#").strip()
                    try:
                        pid = body.split("Number of predicted TMRs")[0].strip()
                        cnt = int(body.rsplit(":", 1)[1].strip())
                        gff_tmr_count[pid] = cnt
                    except Exception:
                        pass

    with open(three[0]) as fh:
        lines = [l.rstrip("\n") for l in fh]
    i = 0
    while i < len(lines):
        if lines[i].startswith(">"):
            header = lines[i][1:]
            seq = lines[i + 1] if i + 1 < len(lines) else ""
            topo = lines[i + 2] if i + 2 < len(lines) else ""
            if "|" in header:
                pid, ptype = [x.strip() for x in header.split("|", 1)]
            else:
                pid, ptype = header.strip(), ""
            regions = []
            if topo:
                run_start = 0
                cur = topo[0]
                for j in range(1, len(topo)):
                    if topo[j] != cur:
                        regions.append({"label": cur, "type": label_names.get(cur, cur),
                                        "start": run_start + 1, "end": j})
                        run_start = j
                        cur = topo[j]
                regions.append({"label": cur, "type": label_names.get(cur, cur),
                                "start": run_start + 1, "end": len(topo)})
            m_runs = sum(1 for r in regions if r["label"] == "M")
            n_tm = gff_tmr_count.get(pid, m_runs)
            sp_regions = [r for r in regions if r["label"] == "S"]
            has_sp = bool(sp_regions) or ("SP" in ptype.upper())
            sp_span = [sp_regions[0]["start"], sp_regions[0]["end"]] if sp_regions else None
            result["proteins"][pid] = {
                "prediction": pred_norm.get(ptype, ptype),
                "prediction_raw": ptype,
                "n_tm_helices": n_tm,
                "topology_string": topo,
                "has_signal_peptide": has_sp,
                "signal_peptide_span": sp_span,
                "length": len(seq) if seq else len(topo),
                "regions": [{"type": r["type"], "start": r["start"], "end": r["end"]} for r in regions],
            }
            i += 3
        else:
            i += 1
    result["n_proteins"] = len(result["proteins"])
    return result
