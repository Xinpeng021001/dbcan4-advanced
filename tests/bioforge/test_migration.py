"""`alembic upgrade head` must produce the same schema as Base.metadata.create_all.

Tests build the schema with create_all (conftest), while production uses Alembic;
this guards against the two drifting apart.
"""
from __future__ import annotations

from pathlib import Path

import pytest
from sqlalchemy import create_engine, inspect

from bioforge.models import Base


def _find_repo_root() -> Path:
    """Walk up to the dir holding alembic.ini (robust to vendored nesting)."""
    here = Path(__file__).resolve()
    for parent in here.parents:
        if (parent / "alembic.ini").is_file():
            return parent
    return here.parents[1]


REPO_ROOT = _find_repo_root()


def _schema(url: str) -> dict[str, tuple[set, set]]:
    insp = inspect(create_engine(url))
    out: dict[str, tuple[set, set]] = {}
    for t in insp.get_table_names():
        if t == "alembic_version":
            continue
        cols = {c["name"] for c in insp.get_columns(t)}
        idx = {i["name"] for i in insp.get_indexes(t)}
        out[t] = (cols, idx)
    return out


def test_alembic_matches_models(tmp_path, monkeypatch):
    pytest.importorskip("alembic")
    from alembic import command
    from alembic.config import Config

    ca_url = f"sqlite:///{tmp_path / 'ca.db'}"
    Base.metadata.create_all(create_engine(ca_url))

    al_url = f"sqlite:///{tmp_path / 'al.db'}"
    monkeypatch.setenv("DATABASE_URL", al_url)
    cfg = Config(str(REPO_ROOT / "alembic.ini"))
    command.upgrade(cfg, "head")

    ca, al = _schema(ca_url), _schema(al_url)
    assert set(ca) == set(al), f"table mismatch: {set(ca) ^ set(al)}"
    for t in ca:
        assert ca[t][0] == al[t][0], f"{t} column mismatch: {ca[t][0] ^ al[t][0]}"
        assert ca[t][1] == al[t][1], f"{t} index mismatch: {ca[t][1] ^ al[t][1]}"
