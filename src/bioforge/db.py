"""SQLAlchemy engine/session helpers shared by the CLI and the API."""
from __future__ import annotations

from sqlalchemy import create_engine, event
from sqlalchemy.engine import Engine
from sqlalchemy.orm import sessionmaker

from .config import database_url, is_sqlite


def make_engine(url: str | None = None) -> Engine:
    url = url or database_url()
    connect_args = {"check_same_thread": False} if url.startswith("sqlite") else {}
    engine = create_engine(url, future=True, connect_args=connect_args)
    if url.startswith("sqlite"):
        # Enforce foreign keys on SQLite (off by default) so cascade/versioning
        # behaves the same as on Postgres.
        @event.listens_for(engine, "connect")
        def _fk_pragma(dbapi_con, _):  # pragma: no cover - trivial
            dbapi_con.execute("PRAGMA foreign_keys=ON")
    return engine


def make_session_factory(engine: Engine | None = None):
    engine = engine or make_engine()
    return sessionmaker(bind=engine, future=True, expire_on_commit=False)
