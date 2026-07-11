"""Resolve the on-disk assets the dbcan4 engine needs (pipeline, scripts, index).

Resolution order for the asset root:
  1. --assets CLI flag / DBCAN4_ASSETS env var
  2. a bundled nf/ inside the installed package (packaged pipeline)
  3. the dbcan4-advanced repo checkout inferred from this file's location
  4. the known met working dir (/array1/xinpeng/dbcan4-advanced)
Heavy data assets (reference embedding index, trained heads) are looked up under
the asset root's results/ + emb/, and are overridable via env vars.
"""
from __future__ import annotations

import os
import sys
from dataclasses import dataclass
from pathlib import Path

_MET_ROOT = Path("/array1/xinpeng/dbcan4-advanced")


@dataclass
class Assets:
    root: Path
    nf_dir: Path
    scripts_dir: Path
    ref_emb_prefix: str
    heads_pt: str
    proj_ref_npz: str
    engine_python: str

    def script(self, name: str) -> str:
        """Absolute path to a stage script (prefers pipeline bin/, then scripts/)."""
        for cand in (self.nf_dir / "bin" / name, self.scripts_dir / name):
            if cand.exists():
                return str(cand)
        # fall back to PATH (Nextflow puts bin/ on PATH inside tasks)
        return name


def _first_existing(*paths: Path) -> Path | None:
    for p in paths:
        if p and p.exists():
            return p
    return None


def resolve_assets(override: str | None = None) -> Assets:
    env_root = os.environ.get("DBCAN4_ASSETS")
    here = Path(__file__).resolve()

    # candidate roots, in priority order
    candidates = []
    if override:
        candidates.append(Path(override))
    if env_root:
        candidates.append(Path(env_root))
    # packaged pipeline: <pkg>/nf
    candidates.append(here.parent)
    # repo checkout: src/dbcan4_advanced/config.py -> repo root (parents[2])
    candidates.append(here.parents[2])
    candidates.append(_MET_ROOT / "repo_clone")
    candidates.append(_MET_ROOT)

    nf_dir = None
    root = None
    for c in candidates:
        cand_nf = _first_existing(c / "nf", c)  # packaged as <pkg>/nf, or c IS the nf dir
        if cand_nf and (cand_nf / "main.nf").exists():
            nf_dir = cand_nf
            root = c
            break
    if nf_dir is None:
        # last resort: met repo_clone/nf
        root = _MET_ROOT / "repo_clone"
        nf_dir = root / "nf"

    scripts_dir = _first_existing(root / "scripts", _MET_ROOT / "scripts") or (root / "scripts")

    ref_emb_prefix = os.environ.get("DBCAN4_REF_EMB", str(_MET_ROOT / "emb" / "ref2024"))
    heads_pt = os.environ.get("DBCAN4_HEADS", str(_MET_ROOT / "results" / "heads" / "heads.pt"))
    proj_ref = os.environ.get("DBCAN4_PROJ_REF", str(_MET_ROOT / "results" / "heads" / "proj_ref.npz"))
    # The GPU engine (embed/infer) needs torch+esm — this lives in the project
    # venv, NOT necessarily the venv the CLI is installed in. Resolve it explicitly.
    _venv_py = _MET_ROOT / "venv" / "bin" / "python"
    engine_python = os.environ.get("DBCAN4_ENGINE_PYTHON",
                                   str(_venv_py) if _venv_py.exists() else sys.executable)

    return Assets(root=Path(root), nf_dir=Path(nf_dir), scripts_dir=Path(scripts_dir),
                  ref_emb_prefix=ref_emb_prefix, heads_pt=heads_pt, proj_ref_npz=proj_ref,
                  engine_python=engine_python)
