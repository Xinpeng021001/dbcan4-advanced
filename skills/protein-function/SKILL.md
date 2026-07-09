---
name: protein-function
description: >-
  Assemble real functional annotation for a protein from its sequence and CAZy
  family — GO terms and InterPro/Pfam domains, EC number and substrate, and
  physicochemistry (MW, pI, aromaticity, instability, GRAVY, aa composition)
  plus sequence features (N-glycosylation sequons, N-terminal signal region).
  Use this skill whenever the task is to functionally annotate or characterize a
  protein or CAZyme, to get GO/EC/substrate/domain evidence for a FASTA sequence,
  to map a CAZy family (GH/GT/PL/CE/CBM/AA) to its EC and substrate, or to compute
  protein physicochemical properties and glycosylation sites — even when the user
  does not name a specific tool. Pairs with the dbcan-annotation and deeploc skills.
---

# protein-function

Assemble **real, sourced** functional evidence for a protein given its sequence
and its CAZy family. This skill never invents annotation values: GO terms and
domains are echoed only from evidence you fetch (InterPro/Pfam via the
`protein-annotation` MCP connector), EC and substrate come only from a dbCAN
family reference you supply, and physicochemistry / sequence features are
computed deterministically from the sequence with Biopython. When a piece of
evidence is genuinely absent (e.g. a family InterPro does not annotate with GO),
report that honestly — do not backfill it.

The `kernel.py` sidecar loads automatically when this skill is loaded; its
functions are then callable in the `python` kernel. Install Biopython first if
missing: `manage_packages(mode="install", environment=<env>, packages=["biopython"])`.

## What you get

`annotate_function(seq, cazy_family, ...)` returns one dict combining:

- **EC number + substrate** — mapped from the CAZy family through a dbCAN
  reference (`famref`). Subfamily labels like `GH5_4` fall back to the base
  family. Alternatives (a family with >1 activity) are listed; `ec_substrate_mapped`
  is `False` (never a guessed EC) when the family is not in the reference.
- **GO terms** — flattened, deduped, grouped by aspect (molecular_function /
  biological_process / cellular_component) from the InterPro/Pfam entries you pass in.
- **Domains** — InterPro/Pfam accessions, names, types.
- **Physicochemistry** — MW, theoretical pI, aromaticity, instability index
  (+stable/unstable class at the 40 cutoff), GRAVY, aa composition (%).
- **N-glycosylation sequons** — N-X-[S/T], X≠P, with 1-based positions.
- **Signal region** — a transparent Kyte-Doolittle N-terminal hydropathy heuristic
  (a prior, explicitly *not* SignalP/DeepLoc).

## Workflow

### 1. Get the sequence and its CAZy family

From a FASTA, a dbCAN `overview.tsv` "Recommend Results" call (see the
`dbcan-annotation` skill), or a record dict. You need the amino-acid sequence
and a family label such as `GH78`.

### 2. Fetch GO + domains from the protein-annotation connector (repl tool)

MCP calls run **only in the `repl` tool**. Two directions are useful:

```python
# repl tool — protein -> domains (needs a UniProt accession):
arch = host.mcp("protein-annotation", "get_domain_architecture", accessions=["P0DTE7"])
# entry detail (InterPro IPRxxxxxx or Pfam PFxxxxx) carries go_terms + signatures:
entry = host.mcp("protein-annotation", "get_interpro_entry", accession="IPR016007")
# find the InterPro entries for a family by keyword:
hits  = host.mcp("protein-annotation", "search_interpro_entries", query="alpha-L-rhamnosidase")
```

Write the collected entry dicts to `./handoff/evidence.json`, then load them in
the `python` kernel. **Do not assume a family carries GO** — many CAZy-family
InterPro entries (e.g. the GH78 α-L-rhamnosidase families) have none. If the
domain route returns no GO, record `n_go_terms: 0` and say so; that is the real
answer, not a gap to fill with a remembered GO id.

### 3. Map the family to EC + substrate

Pass a dbCAN family→EC/substrate reference — a dict or a JSON path — as `famref`.
Its entries look like `{"ec": "3.2.1.40", "substrate": "...", "name": "...",
"substrate_high": "...", "all": [ ... ]}`. `map_family(fam, famref)` handles the
lookup and the subfamily fallback.

### 4. Compute physicochemistry + features (python kernel)

```python
result = annotate_function(
    seq, "GH78",
    interpro_accessions=[d["accession"] for d in domain_entries],
    famref="/path/to/cazy_famref.json",
    interpro_entries=detail_entries,   # dicts carrying go_terms
    domains=collect_domains(domain_entries),
    protein_id="267317",
)
```

`physicochem(seq)`, `n_glyc_sites(seq)`, and `signal_region(seq)` are also
callable standalone. Non-standard residues (X/B/Z/U/O) are dropped before
ProteinAnalysis and reported under `nonstandard_residues_removed`.

### 5. Save the result

`json.dump(result, open("<id>_function.json", "w"), indent=2)` and
`save_artifacts([...], language="python")`. Stamp provenance: which source each
field came from, and any honest gaps (e.g. "InterPro attaches no GO to GH78").

## Honesty

If a live fetch fails, record the real error and status — never substitute a
plausible value. GO/EC/substrate/domains must trace to a real source
(connector response or the supplied reference); physicochemistry/features are
computed from the actual sequence. A result with `n_go_terms: 0` and a clear
note is correct; a result with an invented GO id is not.
