# kernel.py — helpers for the `pfam-scan` skill.
# Runs hmmscan against Pfam-A and parses the --domtblout into an ordered
# domain-architecture list. Pure stdlib; no third-party imports needed.
import os
import shutil
import subprocess

PFAM_HMM_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
PFAM_CLANS_URL = "https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.clans.tsv.gz"
STD_AA = "ACDEFGHIKLMNPQRSTVWY"


def find_hmmer_binary(name="hmmscan", hint=None):
    """Locate an HMMER binary. Checks `hint`, then PATH, then common install dirs."""
    if hint and os.path.exists(hint):
        return hint
    found = shutil.which(name)
    if found:
        return found
    for cand in ("/usr/bin/" + name, "/usr/local/bin/" + name):
        if os.path.exists(cand):
            return cand
    raise FileNotFoundError(
        f"{name} not found on PATH. Install HMMER (conda: `bioconda::hmmer`, "
        f"or apt `hmmer`) or pass an explicit path via the *_bin argument."
    )


def fetch_url(url, dest):
    """Download `url` -> `dest`. Prefers curl (robust retries), falls back to urllib."""
    import urllib.request
    curl = shutil.which("curl")
    if curl:
        r = subprocess.run(
            [curl, "-fSL", "--retry", "3", "-o", dest, url],
            capture_output=True, text=True,
        )
        if r.returncode == 0 and os.path.exists(dest):
            return dest
    urllib.request.urlretrieve(url, dest)
    return dest


def hmmpress_db(pfam_hmm, hmmpress_bin=None, force=False):
    """Build the hmmscan binary index (.h3m/.h3i/.h3f/.h3p) for Pfam-A.hmm.
    Idempotent: skips if the index already exists unless force=True."""
    idx = pfam_hmm + ".h3m"
    if os.path.exists(idx) and not force:
        return {"pressed": False, "index": idx, "reason": "already pressed"}
    if hmmpress_bin is None:
        hmmpress_bin = find_hmmer_binary("hmmpress")
    cmd = [hmmpress_bin] + (["-f"] if force else []) + [pfam_hmm]
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"hmmpress failed (rc={r.returncode}): {r.stderr[-2000:]}")
    return {"pressed": True, "index": idx, "cmd": " ".join(cmd)}


def run_pfam_hmmscan(fasta, pfam_hmm, out, cut_ga=True, cpu=4, evalue=None,
                     hmmscan_bin=None, press_if_needed=True, extra_args=None):
    """Run `hmmscan --domtblout <out>` for `fasta` against `pfam_hmm` (Pfam-A.hmm).

    cut_ga=True applies Pfam per-family gathering thresholds (--cut_ga), the
    curated cutoff that defines Pfam family membership — recommended over a
    flat E-value for domain annotation. Set cut_ga=False and pass `evalue` to
    use a flat -E threshold instead.

    Returns a dict: {returncode, domtblout, cmd, stdout_tail, stderr_tail, ok}.
    Does not raise on a non-zero hmmscan exit — inspect `ok`/`returncode`."""
    if hmmscan_bin is None:
        hmmscan_bin = find_hmmer_binary("hmmscan")
    if press_if_needed and not os.path.exists(pfam_hmm + ".h3m"):
        hmmpress_db(pfam_hmm)
    cmd = [hmmscan_bin, "--domtblout", out]
    if cpu:
        cmd += ["--cpu", str(cpu)]
    if cut_ga:
        cmd += ["--cut_ga"]
    elif evalue is not None:
        cmd += ["-E", str(evalue)]
    if extra_args:
        cmd += list(extra_args)
    cmd += [pfam_hmm, fasta]
    r = subprocess.run(cmd, capture_output=True, text=True)
    return {
        "returncode": r.returncode,
        "domtblout": out,
        "cmd": " ".join(cmd),
        "stdout_tail": r.stdout[-3000:],
        "stderr_tail": r.stderr[-3000:],
        "ok": r.returncode == 0 and os.path.exists(out),
    }


def load_pfam_clans(clans_tsv):
    """Load Pfam-A.clans.tsv -> {pfam_acc_base: {clan_acc, clan_id, pfam_id, desc}}.
    Columns: pfam_acc, clan_acc, clan_id, pfam_id, description (clan_acc may be blank)."""
    m = {}
    with open(clans_tsv) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            acc = parts[0].split(".")[0]
            m[acc] = {
                "clan_acc": parts[1] or None,
                "clan_id": parts[2] or None,
                "pfam_id": parts[3] or None,
                "desc": parts[4] or None,
            }
    return m


def parse_domtblout(path, clan_map=None, order_by="seq_start"):
    """Parse an hmmscan --domtblout file into a domain-architecture list.

    NOTE on hmmscan orientation: each query protein is scanned against the
    profile database, so in the output the TARGET is the Pfam profile and the
    QUERY is the protein. Sequence coordinates therefore come from the ali/env
    (alignment/envelope) columns; HMM coordinates give profile coverage.

    Returns list[dict], one per domain hit, ordered along the sequence
    (by query_name then seq_start) by default. For a multi-sequence FASTA,
    split per protein with group_domains_by_query()."""
    rows = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            f = line.split()
            if len(f) < 22:
                continue
            hmm_len = int(f[2])
            hmm_from, hmm_to = int(f[15]), int(f[16])
            ali_from, ali_to = int(f[17]), int(f[18])
            env_from, env_to = int(f[19]), int(f[20])
            pfam_acc = f[1]
            base = pfam_acc.split(".")[0]
            d = {
                "pfam_acc": pfam_acc,
                "pfam_acc_base": base,
                "name": f[0],
                "query_name": f[3],
                "query_len": int(f[5]),
                "seq_start": env_from,        # envelope coords = domain boundary on protein
                "seq_end": env_to,
                "ali_start": ali_from,        # aligned region (tighter than envelope)
                "ali_end": ali_to,
                "hmm_start": hmm_from,
                "hmm_end": hmm_to,
                "hmm_len": hmm_len,
                "hmm_coverage": round((hmm_to - hmm_from + 1) / hmm_len, 4) if hmm_len else None,
                "full_evalue": float(f[6]),
                "full_score": float(f[7]),
                "c_evalue": float(f[11]),
                "i_evalue": float(f[12]),     # independent E-value — use for per-domain significance
                "bitscore": float(f[13]),
                "bias": float(f[14]),
                "acc_posterior": float(f[21]),
                "domain_num": int(f[9]),
                "domain_total": int(f[10]),
                "description": " ".join(f[22:]) if len(f) > 22 else "",
                "clan": None,
            }
            if clan_map and base in clan_map:
                d["clan"] = clan_map[base].get("clan_acc")
            rows.append(d)
    if order_by:
        rows.sort(key=lambda x: (x["query_name"], x[order_by]))
    return rows


def group_domains_by_query(domains):
    """Group a parse_domtblout() list into {query_name: [domains ordered by seq_start]}."""
    out = {}
    for d in domains:
        out.setdefault(d["query_name"], []).append(d)
    for q in out:
        out[q].sort(key=lambda x: x["seq_start"])
    return out


def architecture_string(domains, sep=" - "):
    """Render an ordered domain list as a compact N->C architecture string,
    e.g. 'Bac_rhamnosid_N - Bac_rhamnosid - Bac_rhamnosid6H'."""
    return sep.join(d["name"] for d in sorted(domains, key=lambda x: x["seq_start"]))


def download_pfam_db(dest_dir, url=None, clans=True, press=True):
    """Fetch Pfam-A.hmm(.gz) into dest_dir, gunzip, hmmpress, and (optionally) the
    clans TSV. Large download (~1-2 GB) and slow; idempotent — skips existing files.
    Returns {pfam_hmm, [clans_tsv], [pressed]}."""
    import gzip
    if url is None:
        url = PFAM_HMM_URL
    os.makedirs(dest_dir, exist_ok=True)
    hmm = os.path.join(dest_dir, "Pfam-A.hmm")
    gz = hmm + ".gz"
    if not os.path.exists(hmm):
        if not os.path.exists(gz):
            fetch_url(url, gz)
        with gzip.open(gz, "rb") as fi, open(hmm, "wb") as fo:
            shutil.copyfileobj(fi, fo)
    result = {"pfam_hmm": hmm}
    if press:
        hmmpress_db(hmm)
        result["pressed"] = True
    if clans:
        ctsv = os.path.join(dest_dir, "Pfam-A.clans.tsv")
        cgz = ctsv + ".gz"
        if not os.path.exists(ctsv):
            fetch_url(PFAM_CLANS_URL, cgz)
            with gzip.open(cgz, "rb") as fi, open(ctsv, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        result["clans_tsv"] = ctsv
    return result
