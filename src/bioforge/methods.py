"""The CAZyme-prediction method registry — one source of truth for both the
ingester and the web UI.

BioForge's baseline ingester loads dbCAN `overview.tsv` calls keyed by a ``tool``
field (``HMMER``/``dbCAN_sub``/``DIAMOND``). The dbCAN4-advanced module adds
protein-language-model and structure-based predictors. Rather than scatter the
"which tools are advanced / what colour / how to label them" knowledge across
parsers and templates, every method is described exactly once here.

Each entry:
    key          canonical value stored in ``CazymeAnnotation.tool``
    display      human-readable label for the UI
    family       'baseline' | 'advanced'      (the method_family column)
    kind         hmm|diamond|dbcan_sub|sequence-plm|structure|fusion (method_kind)
    signal       'sequence' | 'structure' | 'fusion'  (what evidence it uses)
    colour       hex for UI tags/badges (baseline greys, advanced accent hues)
    blurb        one-line description surfaced in tooltips / the /advanced page
    order        stable sort order in side-by-side tables

The design rationale (design_dbcan4_advanced.md §6): baseline detects homology in
*sequence* space (HMMER profiles, DIAMOND alignment); the advanced tier adds
*embedding* similarity (ESM-C kNN / nearest-centroid / contrastive head) and
*structure* similarity (Foldseek vs CAZyme3D, SaProt structure-aware embedding),
then a *fusion* call that combines the orthogonal signals. Agreement across
orthogonal signals is the strongest evidence for a genuine remote-homolog CAZyme.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Method:
    key: str
    display: str
    family: str          # 'baseline' | 'advanced'
    kind: str            # hmm|diamond|dbcan_sub|sequence-plm|structure|fusion
    signal: str          # 'sequence' | 'structure' | 'fusion'
    colour: str          # hex
    blurb: str
    order: int

    @property
    def is_advanced(self) -> bool:
        return self.family == "advanced"


# --- baseline tools (dbCAN3/4 overview.tsv) -------------------------------
_BASELINE = [
    Method("HMMER", "HMMER", "baseline", "hmm", "sequence", "#64748b",
           "Profile-HMM search vs dbCAN CAZyme HMMs (sequence homology).", 10),
    Method("dbCAN_sub", "dbCAN_sub", "baseline", "dbcan_sub", "sequence", "#64748b",
           "Subfamily-level HMM assignment (dbCAN-sub).", 11),
    Method("DIAMOND", "DIAMOND", "baseline", "diamond", "sequence", "#64748b",
           "DIAMOND BLASTp vs the CAZy sequence database (sequence homology).", 12),
]

# --- advanced tools (dbCAN4-advanced module) ------------------------------
_ADVANCED = [
    Method("ESM-C-kNN", "ESM-C · kNN", "advanced", "sequence-plm", "sequence",
           "#2563eb",
           "ESM-C 600M embedding, k-nearest-reference majority vote; "
           "confidence = top-1 cosine similarity.", 20),
    Method("ESM-C-centroid", "ESM-C · centroid", "advanced", "sequence-plm",
           "sequence", "#7c3aed",
           "ESM-C embedding, nearest per-family centroid (CLEAN-style prototype); "
           "confidence = best centroid cosine.", 21),
    Method("ESM-C-contrastive", "ESM-C · contrastive", "advanced", "sequence-plm",
           "sequence", "#c026d3",
           "Supervised-contrastive projection of ESM-C embeddings, then "
           "centroid/kNN in the learned space (CLEAN mechanism).", 22),
    Method("Foldseek-CAZyme3D", "Foldseek · CAZyme3D", "advanced", "structure",
           "structure", "#059669",
           "Foldseek 3Di structural search vs the CAZyme3D structure database; "
           "family label from best structural hit (resolves the sub-25%-identity "
           "twilight zone).", 23),
    Method("SaProt", "SaProt", "advanced", "structure", "structure", "#0d9488",
           "Structure-aware pLM (residue + Foldseek 3Di token) embedding "
           "retrieval — a second, orthogonal structure signal.", 24),
    Method("fusion", "Fusion (consensus)", "advanced", "fusion", "fusion",
           "#ea580c",
           "Calibrated combination of sequence-pLM and structure evidence into "
           "one call + confidence; agreement across orthogonal signals is the "
           "strongest remote-homolog evidence.", 25),
]

REGISTRY: dict[str, Method] = {m.key: m for m in (_BASELINE + _ADVANCED)}

BASELINE_TOOLS = [m.key for m in _BASELINE]
ADVANCED_TOOLS = [m.key for m in _ADVANCED]


# --- feature-annotation tools (produce ProteinFeature rows, NOT family calls) -
# These are the comprehensive per-protein annotation tools (OUTPUT_CONTRACT §2.2–2.9).
# Kept separate from the CAZyme-prediction REGISTRY so they never pollute the
# advanced-vs-baseline family comparison, while still giving the UI a real
# label / colour / blurb / status for each tool.
@dataclass(frozen=True)
class FeatureTool:
    key: str             # value stored in ProteinFeature.tool
    display: str
    feature_type: str    # signal_peptide|tm_topology|structure|domain|structure_hit|localization|physicochem|ec_prediction
    status: str          # 'real_tool' | 'scaffold' (license-gated, not installed here)
    colour: str
    blurb: str


_FEATURE_TOOLS = [
    FeatureTool("Pfam/hmmscan", "Pfam · hmmscan", "domain", "real_tool", "#3a5bd0",
                "HMMER hmmscan vs Pfam-A with curated gathering thresholds; full domain architecture (coords + E-values)."),
    FeatureTool("ESMFold", "ESMFold", "structure", "real_tool", "#059669",
                "Single-sequence 3D structure prediction (facebook/esmfold_v1); per-residue pLDDT."),
    FeatureTool("Foldseek-CAZyme3D", "Foldseek · CAZyme3D", "structure_hit", "real_tool", "#0d9488",
                "Structural-homology search of the predicted fold vs the CAZyme3D reference set (TM-score / LDDT)."),
    FeatureTool("DeepTMHMM", "DeepTMHMM", "tm_topology", "real_tool", "#7c3aed",
                "Transmembrane topology + signal-peptide prediction (DTU DeepTMHMM via BioLib)."),
    FeatureTool("SignalP6", "SignalP 6.0", "signal_peptide", "scaffold", "#b25c00",
                "Signal-peptide type + cleavage site (DTU SignalP-6.0). License-gated; scaffold runs it if installed."),
    FeatureTool("CLEAN", "CLEAN", "ec_prediction", "real_tool", "#c026d3",
                "Sequence-based EC-number prediction via contrastive learning (Yu et al., Science 2023) — independent of family-inherited EC."),
    FeatureTool("Biopython", "Biopython", "physicochem", "real_tool", "#2563eb",
                "Physicochemical properties (MW, pI, instability, GRAVY, aromaticity, aa composition)."),
    FeatureTool("DeepLoc", "DeepLoc-2.0", "localization", "scaffold", "#b25c00",
                "Subcellular localization (DTU DeepLoc-2.0). License-gated; localization otherwise derived from signal peptide + GO-CC."),
]

FEATURE_TOOLS: dict[str, FeatureTool] = {m.key: m for m in _FEATURE_TOOLS}


def feature_tool(key: str | None) -> FeatureTool | None:
    """Look up a FeatureTool by its stored ``ProteinFeature.tool`` key."""
    if key is None:
        return None
    return FEATURE_TOOLS.get(key)


def feature_colour_of(key: str | None) -> str:
    m = feature_tool(key)
    return m.colour if m else "#64748b"


def feature_display_of(key: str | None) -> str:
    m = feature_tool(key)
    return m.display if m else (key or "—")


def get(tool: str | None) -> Method | None:
    """Look up a Method by its stored ``tool`` key (exact, case-sensitive)."""
    if tool is None:
        return None
    return REGISTRY.get(tool)


def family_of(tool: str | None) -> str:
    """'baseline' | 'advanced' | 'unknown' for a tool key (never raises)."""
    m = get(tool)
    return m.family if m else "unknown"


def kind_of(tool: str | None) -> str | None:
    m = get(tool)
    return m.kind if m else None


def is_advanced(tool: str | None) -> bool:
    m = get(tool)
    return bool(m and m.is_advanced)


def display_of(tool: str | None) -> str:
    m = get(tool)
    return m.display if m else (tool or "—")


def colour_of(tool: str | None) -> str:
    m = get(tool)
    return m.colour if m else "#64748b"


def as_dict(m: Method) -> dict:
    return {
        "key": m.key, "display": m.display, "family": m.family, "kind": m.kind,
        "signal": m.signal, "colour": m.colour, "blurb": m.blurb,
        "order": m.order, "is_advanced": m.is_advanced,
    }
