#!/usr/bin/env bash
# Daily-update helper. Usage: GITHUB_TOKEN=... ./scripts/daily_push.sh "Day N: summary"
# Commits the day's work and pushes. Assumes docs/daily-log.md already has today's entry.
set -euo pipefail
MSG="${1:?provide a commit message, e.g. 'Day 2: baseline + splits'}"
cd "$(dirname "$0")/.."
git add -A
if git diff --cached --quiet; then echo "nothing to commit"; exit 0; fi
git commit -m "$MSG"
if [ -n "${GITHUB_TOKEN:-}" ]; then
  git push "https://x-access-token:${GITHUB_TOKEN}@github.com/Xinpeng021001/dbcan4-advanced.git" main
else
  git push origin main
fi
echo "pushed: $(git rev-parse --short HEAD)"
