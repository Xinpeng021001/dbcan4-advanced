"""Shared pytest fixtures: a fresh in-file SQLite DB with the schema applied."""
from __future__ import annotations

from pathlib import Path

import pytest
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from bioforge.models import Base


def _find_sample_data() -> Path:
    """Locate the vendored sample_data dir by walking up from this file.

    Robust to how deep these tests are nested (upstream they sit at
    <biodb>/tests/; vendored here they sit at <repo>/tests/bioforge/).
    """
    here = Path(__file__).resolve()
    for parent in here.parents:
        cand = parent / "sample_data"
        if cand.is_dir():
            return cand
    return here.parents[1] / "sample_data"


SAMPLE_DATA = _find_sample_data()


@pytest.fixture()
def engine(tmp_path):
    url = f"sqlite:///{tmp_path / 'test.db'}"
    eng = create_engine(url, future=True, connect_args={"check_same_thread": False})
    Base.metadata.create_all(eng)  # schema == models (kept in sync with Alembic)
    return eng


@pytest.fixture()
def Session(engine):
    return sessionmaker(bind=engine, future=True, expire_on_commit=False)


@pytest.fixture()
def sample_data() -> Path:
    return SAMPLE_DATA
