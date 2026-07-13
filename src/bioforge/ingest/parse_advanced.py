"""Parser for the dbCAN4-advanced STANDARD OUTPUT CONTRACT.

The advanced annotation pipeline (github.com/Xinpeng021001/dbcan4-advanced,
`nf/OUTPUT_CONTRACT.md`) publishes a `cazyme_advanced/manifest.json` that points
at normalized per-method prediction TSVs and per-protein feature TSVs. This
module turns that manifest into typed records the loader inserts:

    parse_manifest(path)   -> AdvancedManifest (samples, tool versions, release meta)
    read_predictions(...)  -> Iterator[AdvancedCall]     (one per protein per tool)
    read_features(...)     -> Iterator[AdvancedFeature]  (SignalP6/DeepTMHMM/structure)

Two TSV modes, both supported:
  * normalized §2.1  — columns: protein_id, family, confidence, ec, all_families, extra
  * raw/legacy       — the manifest entry carries id_col/family_col/confidence_col
                       overrides so an arbitrary wide TSV maps onto the standard
                       fields without new code (e.g. the prototype's
                       esmc_retrieval_pred.tsv with cent_pred/cent_conf).

Every `tool` in the manifest must exist in ``bioforge.methods.REGISTRY`` — an
unknown tool raises, exactly like the baseline dbCAN parser rejects an
unrecognised column layout, so a contract drift fails loudly.
"""
from __future__ import annotations

import csv
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

from ..methods import REGISTRY


@dataclass
class PredFile:
    tool: str
    path: Path
    id_col: str | None = None
    family_col: str | None = None
    confidence_col: str | None = None


@dataclass
class FeatFile:
    feature_type: str
    tool: str
    path: Path


@dataclass
class SampleManifest:
    sample_key: str
    predictions: list[PredFile] = field(default_factory=list)
    features: list[FeatFile] = field(default_factory=list)


@dataclass
class AdvancedManifest:
    contract_version: str
    pipeline: str
    pipeline_version: str | None
    release_label: str
    release_notes: str | None
    tool_versions: dict
    samples: list[SampleManifest]
    base_dir: Path              # dir the relative paths resolve against
    manifest_path: Path


@dataclass
class AdvancedCall:
    protein_id: str
    cazy_family: str
    confidence: float | None
    ec: str | None
    tool: str
    method_kind: str | None
    all_families: str | None
    extra: dict


@dataclass
class AdvancedFeature:
    protein_id: str
    feature_type: str
    tool: str
    label: str | None
    score: float | None
    start: int | None
    end: int | None
    structure_rel_path: str | None
    attributes: dict


_NULL = {"", "-", "na", "none", "nan", None}


def _s(v):
    if v is None:
        return None
    v = str(v).strip()
    return None if v.lower() in {"", "-", "na", "none", "nan"} else v


def _f(v):
    v = _s(v)
    if v is None:
        return None
    try:
        return float(v)
    except ValueError:
        return None


def _i(v):
    v = _s(v)
    if v is None:
        return None
    try:
        return int(float(v))
    except ValueError:
        return None


def _jsonobj(v) -> dict:
    v = _s(v)
    if v is None:
        return {}
    try:
        obj = json.loads(v)
        return obj if isinstance(obj, dict) else {"value": obj}
    except (json.JSONDecodeError, ValueError):
        return {}


def parse_manifest(path: str | Path) -> AdvancedManifest:
    path = Path(path)
    if path.is_dir():  # tolerate being handed the cazyme_advanced/ dir
        path = path / "manifest.json"
    data = json.loads(path.read_text())
    base = path.parent
    samples: list[SampleManifest] = []
    for s in data.get("samples", []):
        preds = []
        for p in s.get("cazyme_predictions", []):
            tool = p["tool"]
            if tool not in REGISTRY:
                raise ValueError(
                    f"manifest {path}: unknown tool {tool!r} not in "
                    f"bioforge.methods.REGISTRY ({list(REGISTRY)})"
                )
            preds.append(PredFile(
                tool=tool, path=base / p["path"],
                id_col=p.get("id_col"), family_col=p.get("family_col"),
                confidence_col=p.get("confidence_col")))
        feats = [
            FeatFile(feature_type=f["feature_type"], tool=f["tool"],
                     path=base / f["path"])
            for f in s.get("protein_features", [])
        ]
        samples.append(SampleManifest(
            sample_key=s["sample_key"], predictions=preds, features=feats))
    return AdvancedManifest(
        contract_version=str(data.get("contract_version", "?")),
        pipeline=data.get("pipeline", "dbcan4-advanced"),
        pipeline_version=data.get("pipeline_version"),
        release_label=data.get("release_label")
        or f"advanced-{path.parent.name}",
        release_notes=data.get("release_notes"),
        tool_versions=data.get("tool_versions", {}) or {},
        samples=samples,
        base_dir=base,
        manifest_path=path,
    )


def read_predictions(pf: PredFile) -> Iterator[AdvancedCall]:
    """Yield one AdvancedCall per non-abstaining row of a prediction TSV.

    Rows whose family is null/`-` (an abstention) are skipped — a CAZyme
    annotation needs a family — but callers may count them separately."""
    meth = REGISTRY[pf.tool]
    id_col = pf.id_col or "protein_id"
    fam_col = pf.family_col or "family"
    conf_col = pf.confidence_col or "confidence"
    with pf.path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if id_col not in (reader.fieldnames or []):
            raise ValueError(
                f"{pf.path}: id column {id_col!r} not found; "
                f"columns={reader.fieldnames}")
        for row in reader:
            pid = _s(row.get(id_col))
            fam = _s(row.get(fam_col))
            if pid is None or fam is None:
                continue  # abstain / missing
            extra = _jsonobj(row.get("extra"))
            # In raw mode there's no `extra`; capture any non-standard columns.
            if not extra and (pf.family_col or pf.confidence_col):
                std = {id_col, fam_col, conf_col, "ec", "all_families", "extra"}
                extra = {k: v for k, v in row.items()
                         if k not in std and _s(v) is not None}
            yield AdvancedCall(
                protein_id=pid,
                cazy_family=fam,
                confidence=_f(row.get(conf_col)),
                ec=_s(row.get("ec")),
                tool=pf.tool,
                method_kind=meth.kind,
                all_families=_s(row.get("all_families")),
                extra=extra,
            )


def read_features(ff: FeatFile) -> Iterator[AdvancedFeature]:
    """Yield AdvancedFeature rows from a §2.2/2.3/2.4 feature TSV.

    The feature TSVs have per-type column names (see OUTPUT_CONTRACT §2); we map
    each type's columns onto the generic ProteinFeature shape."""
    with ff.path.open(newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = set(reader.fieldnames or [])
        for row in reader:
            pid = _s(row.get("protein_id"))
            if pid is None:
                continue
            ft = ff.feature_type
            if ft == "signal_peptide":
                yield AdvancedFeature(
                    pid, ft, ff.tool,
                    label=_s(row.get("prediction")),
                    score=_f(row.get("sp_prob")),
                    start=None, end=_i(row.get("cs_position")),
                    structure_rel_path=None,
                    attributes=_jsonobj(row.get("extra")))
            elif ft == "tm_topology":
                attrs = _jsonobj(row.get("extra"))
                if _s(row.get("topology")):
                    attrs["topology"] = _s(row.get("topology"))
                yield AdvancedFeature(
                    pid, ft, ff.tool,
                    label=_s(row.get("prediction")),
                    score=_f(row.get("n_tm")),
                    start=None, end=None,
                    structure_rel_path=None, attributes=attrs)
            elif ft == "structure":
                attrs = _jsonobj(row.get("extra"))
                if _s(row.get("length")):
                    attrs["length"] = _i(row.get("length"))
                if _s(row.get("source")):
                    attrs["source"] = _s(row.get("source"))
                yield AdvancedFeature(
                    pid, ft, ff.tool,
                    label=_s(row.get("source")),
                    score=_f(row.get("plddt")),
                    start=None, end=None,
                    structure_rel_path=_s(row.get("path")),
                    attributes=attrs)
            elif ft == "domain":
                # §2.5 Pfam/hmmscan — one row per domain occurrence.
                attrs = _jsonobj(row.get("extra"))
                for k in ("acc", "evalue"):
                    if _s(row.get(k)):
                        attrs[k] = _s(row.get(k))
                yield AdvancedFeature(
                    pid, ft, ff.tool,
                    label=_s(row.get("name")),
                    score=_f(row.get("score")),
                    start=_i(row.get("start")), end=_i(row.get("end")),
                    structure_rel_path=None, attributes=attrs)
            elif ft == "structure_hit":
                # §2.6 Foldseek — one row per structural-homology hit.
                attrs = _jsonobj(row.get("extra"))
                for k in ("target", "prob", "lddt", "evalue"):
                    if _s(row.get(k)):
                        attrs[k] = _s(row.get(k))
                yield AdvancedFeature(
                    pid, ft, ff.tool,
                    label=_s(row.get("target_family")),
                    score=_f(row.get("tmscore")),
                    start=None, end=None,
                    structure_rel_path=None, attributes=attrs)
            elif ft == "localization":
                # §2.7 DeepLoc / derived.
                attrs = _jsonobj(row.get("extra"))
                if _s(row.get("method")):
                    attrs["method"] = _s(row.get("method"))
                yield AdvancedFeature(
                    pid, ft, ff.tool,
                    label=_s(row.get("localization")),
                    score=_f(row.get("confidence")),
                    start=None, end=None,
                    structure_rel_path=None, attributes=attrs)
            elif ft == "physicochem":
                # §2.8 Biopython — one row per protein (summary).
                attrs = _jsonobj(row.get("extra"))
                for k in ("pi", "instability", "gravy", "aromaticity"):
                    if _s(row.get(k)):
                        attrs[k] = _f(row.get(k))
                yield AdvancedFeature(
                    pid, ft, ff.tool,
                    label="physicochem",
                    score=_f(row.get("mw_da")),
                    start=None, end=None,
                    structure_rel_path=None, attributes=attrs)
            elif ft == "ec_prediction":
                # §2.9 CLEAN — independent sequence-based EC prediction.
                attrs = _jsonobj(row.get("extra"))
                if _s(row.get("rank")):
                    attrs["rank"] = _i(row.get("rank"))
                yield AdvancedFeature(
                    pid, ft, ff.tool,
                    label=_s(row.get("ec_number")),
                    score=_f(row.get("confidence")),
                    start=None, end=None,
                    structure_rel_path=None, attributes=attrs)
            else:
                # Unknown feature type: keep the whole row as attributes.
                yield AdvancedFeature(
                    pid, ft, ff.tool, label=None, score=None,
                    start=None, end=None, structure_rel_path=None,
                    attributes={k: v for k, v in row.items() if _s(v)})
