---
name: pfam-scan
description: >
  Annotate the Pfam domain architecture of a protein by running hmmscan
  (HMMER3) against the Pfam-A profile database and parsing the --domtblout into
  an ordered, per-domain list (Pfam accession, name, sequence start/end, HMM
  coverage, E-value, bitscore, clan). Use this skill whenever the task involves
  Pfam domains, hmmscan/hmmsearch, HMMER domain tables, "what domains does this
  protein have", domain architecture, multi-domain / modular protein layout, or
  cross-checking a family call (e.g. a CAZyme GH/GT/PL/CE assignment) against
  independent HMM evidence — even if the user does not say "Pfam" or "hmmscan"
  by name. Works for any protein FASTA (single or multi-sequence).
---

# Pfam domain annotation via hmmscan

Pfam domains are the standard vocabulary for describing what a protein is built
from. Running `hmmscan` against the Pfam-A profile library and reading the
per-domain table gives you the **domain architecture** — which conserved
modules occur, where they sit on the sequence, and in what N→C order. That
architecture is often more informative than a single top-hit family label: it
distinguishes catalytic from accessory (carbohydrate-binding, structural)
modules, and it is independent HMM evidence you can use to confirm or challenge
a call made by another method (DIAMOND/BLAST best-hit, a CAZyme family
assignment, etc.).

This skill runs the scan and turns the raw output into a clean list you can
reason over. Helper functions are pre-loaded from `kernel.py` when the skill
loads — you do not need to re-implement HMMER parsing.

## When to reach for this

- "What Pfam domains / what domain architecture does this protein have?"
- Cross-checking a family/function call against independent evidence.
- Building a domain-architecture column for a set of proteins.
- Any mention of hmmscan, HMMER, Pfam-A, `--domtblout`, or domain tables.

## Prerequisites

- **HMMER3** (`hmmscan`, `hmmpress`) on PATH. Check with
  `find_hmmer_binary("hmmscan")`. Install via conda (`bioconda::hmmer`) or the
  system package manager (`apt install hmmer`) if missing.
- **Pfam-A.hmm**, decompressed and `hmmpress`-ed (produces the
  `.h3m/.h3i/.h3f/.h3p` binary index hmmscan actually reads). This is a large
  (~1.5–2 GB decompressed), one-time download; reuse it across runs. Optional
  but recommended: **Pfam-A.clans.tsv** to attach clan membership.

## Workflow

The kernel plugin exposes these helpers (all pure-stdlib, no heavy deps):

| helper | purpose |
|---|---|
| `find_hmmer_binary(name, hint=None)` | locate `hmmscan`/`hmmpress`/`hmmpress` |
| `download_pfam_db(dest_dir)` | fetch + gunzip + hmmpress Pfam-A (+ clans) |
| `hmmpress_db(pfam_hmm)` | build the binary index (idempotent) |
| `run_pfam_hmmscan(fasta, pfam_hmm, out, ...)` | run the scan → `--domtblout` |
| `parse_domtblout(path, clan_map=None)` | → ordered `list[dict]`, one per domain |
| `load_pfam_clans(clans_tsv)` | `{pfam_acc: clan}` map for clan annotation |
| `group_domains_by_query(domains)` | split a multi-FASTA parse per protein |
| `architecture_string(domains)` | render a compact `N→C` architecture string |

### 1. Make sure the database exists (one-time, slow)

```python
db = download_pfam_db("/path/to/pfam")   # ~1-2 GB, minutes; idempotent
# db == {"pfam_hmm": ".../Pfam-A.hmm", "clans_tsv": ".../Pfam-A.clans.tsv", "pressed": True}
```

If you already have `Pfam-A.hmm`, just ensure it is pressed:
`hmmpress_db("/path/to/pfam/Pfam-A.hmm")`.

On a remote host with a package manager, prefer downloading directly there
(EBI: `https://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz`),
`gunzip`, then `hmmpress Pfam-A.hmm`, so the multi-GB file never crosses the
wire twice.

### 2. Run the scan

```python
res = run_pfam_hmmscan("protein.fasta", db["pfam_hmm"], "protein.domtblout")
assert res["ok"], res["stderr_tail"]
```

**Thresholding — use `--cut_ga` (the default here).** Pfam curates a per-family
*gathering threshold* (GA) that defines membership in that family; applying it
(`cut_ga=True`) is how Pfam/InterPro themselves decide what counts as a real
domain, and it avoids the arbitrariness of a flat E-value across families of
very different lengths and conservation. Use a flat E-value
(`cut_ga=False, evalue=1e-5`) only when you deliberately want to see
sub-threshold/borderline hits.

### 3. Parse into a domain-architecture list

```python
clans = load_pfam_clans(db["clans_tsv"])          # optional
domains = parse_domtblout("protein.domtblout", clan_map=clans)
```

Each element is a dict with the fields the task asks for:
`pfam_acc`, `pfam_acc_base`, `name`, `query_name`, `seq_start`, `seq_end`
(protein envelope coordinates), `ali_start/ali_end`, `hmm_start/hmm_end`,
`hmm_len`, `hmm_coverage`, `full_evalue`, `i_evalue` (per-domain independent
E-value — the right significance number for an individual domain), `bitscore`,
`acc_posterior`, `clan`, `description`. The list is ordered along the sequence.

For a multi-sequence FASTA: `by_prot = group_domains_by_query(domains)` →
`{protein_id: [domains...]}`, and `architecture_string(doms)` renders e.g.
`Bac_rhamnosid_N - Bac_rhamnosid - Bac_rhamnosid6H`.

### 4. Interpret / cross-check

- Read the architecture N→C: catalytic module(s) plus any accessory domains.
- Overlapping hits from the **same clan** are usually the same region matched
  by related profiles — keep the higher-bitscore one when reporting a single
  architecture, or note the redundancy.
- When confirming another method's call, check whether the Pfam domains are
  *consistent with* that family's known architecture. Agreement across
  independent evidence (Pfam + InterPro + the original caller) is a strong
  signal; a mismatch is worth flagging, not silently resolving.

## Output shape

Persist the parsed result as JSON — a list of domain dicts (or
`{protein_id: [domains]}` for a batch) — plus keep the raw `.domtblout` so the
annotation is reproducible and auditable. Report the domain **count**, the
ordered **architecture string**, and the top domains with their coordinates and
E-values.

## Honesty

`hmmscan` output is empirical. Report exactly what the scan returns — real
accessions, coordinates, E-values, bitscores. If the database is missing or the
scan fails, say so with the real error; never invent domains, coordinates, or
scores to fill a gap.
