"""Runtime configuration, resolved from environment variables.

Postgres is the primary target; SQLite is the zero-setup dev fallback used when
DATABASE_URL is unset. Everything else (paths for served JBrowse tracks) is
derived here so the ingester, API, and Alembic all agree on one source of truth.
"""
from __future__ import annotations

import os
from pathlib import Path

# Repo root = two levels up from this file's package (src/bioforge/config.py -> repo).
REPO_ROOT = Path(__file__).resolve().parents[2]


def database_url() -> str:
    """Return the SQLAlchemy URL. Defaults to a local SQLite dev DB."""
    url = os.environ.get("DATABASE_URL")
    if url:
        return url
    # SQLite dev fallback — file lives at repo root so CLI + API share it.
    return f"sqlite:///{REPO_ROOT / 'bioforge_dev.db'}"


def tracks_dir() -> Path:
    """Directory where per-sample JBrowse track files + config.json are written.

    Served as static files by the web layer; JBrowse reads bgzip/tabix over HTTP
    range requests, so this must be web-served, not DB-stored.
    """
    d = Path(os.environ.get("BIOFORGE_TRACKS_DIR", REPO_ROOT / "web_static" / "tracks"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def is_sqlite(url: str | None = None) -> bool:
    return (url or database_url()).startswith("sqlite")
