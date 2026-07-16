#!/usr/bin/env bash
# ============================================================================
# seed-series.sh — Reproducible seed for the series-builder heatmap harness.
#
# Adds ONE indexed, ordered SERIES to a writable copy of the Himalaya dev DB so
# the series-builder heatmap (issue #258) renders >=2 rows with its outer-hair
# keyline frame + rotated ordering-variable axis label. Closes the visual-
# harness data gap: the dev DB has series but none indexed/ordered.
#
# ---------------------------------------------------------------------------
# WRITABLE COPY ONLY.  This script MUTATES the database it talks to (it creates
# a series + commits a plate). NEVER point it at the original dev DB. Copy the
# dev DB to a temp file first, e.g.:
#     cp ~/projects/himalaya-devdata/himalaya.db /tmp/seed-harness.db
# then `serve` the copy and run this script against that backend.
# ---------------------------------------------------------------------------
#
# API-DRIVEN: every write goes through the live backend (apply_event! pipeline)
# so the schema stays correct. Two calls:
#   1) POST /api/series             -> creates the series (recipe + ordering_variable)
#   2) POST /api/series/{id}/commit -> writes the plate (series_members).
#      Per-member snapshots (incl. confirmed_index) are auto-filled server-side
#      from each exposure's confirmed custom-active index via
#      compute_member_snapshot — the client only supplies {exposure_id, display_order}.
#
# DYNAMIC EXPOSURE DISCOVERY (re-ingest safe): the member exposure ids are
# discovered at runtime by querying the served SQLite DB directly, so the seed
# survives a re-ingest (no hardcoded ids). An exposure qualifies as "confirmed"
# when it has an index_groups row with:
#     kind = 'custom'  AND  active = 1
# joined to an index whose  r_squared >= 0.98  (the CONFIRMED_INDEX_R2_GATE in
# packages/HimalayaUI/src/comparisons.jl). One exposure is taken per distinct
# sample (so the rows are distinct samples), up to N=4. The script FAILS LOUDLY
# if fewer than 2 qualifying exposures exist (a 1-row heatmap is not a harness).
#
# NOT fully idempotent: each run mints a NEW series id (POST). Run against a
# FRESH copy of the dev DB for a clean single-seeded-series result. The script
# prints the new SERIES_ID; use that for the /series/<id> URL.
#
# Usage:
#   HIMALAYA_DB_PATH=/tmp/seed-harness.db ./seed-series.sh [BASE_URL]
#     BASE_URL defaults to http://127.0.0.1:8094
#
#   HIMALAYA_DB_PATH must point at the SAME writable DB file the backend is
#   serving (the script reads it directly with sqlite3 for discovery). If unset,
#   it falls back to ~/projects/himalaya-devdata/himalaya.db — which is the
#   READ-ONLY original, so discovery still works but you MUST be serving a copy.
#
# Requires: backend already serving the writable DB copy with /Volumes/data
# mounted (traces load from *_tot.dat files on that volume); sqlite3, curl,
# python3 on PATH.
# ============================================================================
set -euo pipefail

BASE="${1:-http://127.0.0.1:8094}"
DB="${HIMALAYA_DB_PATH:-$HOME/projects/himalaya-devdata/himalaya.db}"
N=4   # desired number of distinct-sample members (>=2 required)

command -v sqlite3 >/dev/null || { echo "[seed] FATAL: sqlite3 not on PATH" >&2; exit 1; }
[ -f "$DB" ] || { echo "[seed] FATAL: DB not found at HIMALAYA_DB_PATH=$DB" >&2; exit 1; }

# ── 1. Discover confirmed exposures (one per distinct sample) ───────────────
echo "[seed] discovering up to $N confirmed exposures in $DB"
DISCOVERY_SQL="
SELECT e.sample_id, MIN(g.exposure_id) AS exposure_id
FROM index_groups g
JOIN index_group_members m ON m.group_id = g.id
JOIN indices i ON i.id = m.index_id
JOIN exposures e ON e.id = g.exposure_id
WHERE g.kind = 'custom' AND g.active = 1
  AND i.r_squared IS NOT NULL AND i.r_squared >= 0.98
GROUP BY e.sample_id
ORDER BY e.sample_id
LIMIT $N;"

# Read pipe-separated rows: sample_id|exposure_id
SAMPLE_IDS=()
EXPOSURE_IDS=()
while IFS='|' read -r sid eid; do
  [ -z "$sid" ] && continue
  SAMPLE_IDS+=("$sid")
  EXPOSURE_IDS+=("$eid")
done < <(sqlite3 -noheader -separator '|' "$DB" "$DISCOVERY_SQL")

COUNT=${#EXPOSURE_IDS[@]}
echo "[seed] discovered $COUNT qualifying exposure(s):"
for i in "${!EXPOSURE_IDS[@]}"; do
  echo "  sample_id=${SAMPLE_IDS[$i]} -> exposure_id=${EXPOSURE_IDS[$i]}"
done

if [ "$COUNT" -lt 2 ]; then
  echo "[seed] FATAL: only $COUNT qualifying exposure(s) found; need >=2 for a heatmap." >&2
  echo "[seed]        criteria: index_groups kind='custom' active=1, index r_squared>=0.98." >&2
  exit 1
fi

# ── 2. Build the JSON arrays from the discovered ids ────────────────────────
SAMPLES_JSON=""
MEMBERS_JSON=""
for i in "${!EXPOSURE_IDS[@]}"; do
  sep=$([ "$i" -gt 0 ] && echo "," || echo "")
  SAMPLES_JSON+="${sep}{\"sample_id\": ${SAMPLE_IDS[$i]}, \"position\": $i}"
  MEMBERS_JSON+="${sep}{\"exposure_id\": ${EXPOSURE_IDS[$i]}, \"display_order\": $i}"
done

# ── 3. Create the series (recipe + ordering_variable) ───────────────────────
echo "[seed] creating series via POST $BASE/api/series"
CREATE_RESP=$(curl -sS -X POST "$BASE/api/series" \
  -H 'Content-Type: application/json' \
  -d "{
        \"title\": \"Dose ladder (heatmap harness)\",
        \"description\": \"Seeded indexed+ordered series for the #258 series-builder heatmap visual harness. ${COUNT} confirmed-index samples ordered by dose.\",
        \"ordering_variable\": \"dose\",
        \"order_rule\": \"ascending\",
        \"view_grouping_mode\": \"byPhase\",
        \"samples\": [${SAMPLES_JSON}]
      }")
SERIES_ID=$(CREATE_RESP="$CREATE_RESP" python3 -c 'import os,json;print(json.loads(os.environ["CREATE_RESP"])["id"])')
echo "[seed] created series id=$SERIES_ID"

# ── 4. Commit the plate (exposure_id is the only required per-member field) ──
echo "[seed] committing plate via POST $BASE/api/series/$SERIES_ID/commit"
COMMIT_RESP=$(curl -sS -X POST "$BASE/api/series/$SERIES_ID/commit" \
  -H 'Content-Type: application/json' \
  -d "{\"members\": [${MEMBERS_JSON}]}")

# ── 5. Self-verify the committed plate (JSON via env, not stdin) ────────────
COMMIT_RESP="$COMMIT_RESP" SERIES_ID="$SERIES_ID" python3 <<'PY'
import os, json
data = json.loads(os.environ["COMMIT_RESP"])
sid = os.environ["SERIES_ID"]
members = data.get("members", [])
ov = data.get("ordering_variable")
print(f"[seed] series {sid}: {len(members)} members, ordering_variable={ov!r}")
for m in members:
    ci = (m.get("snapshot") or {}).get("confirmed_index") or {}
    print(f"  member exposure_id={m.get('exposure_id')} "
          f"display_order={m.get('display_order')} "
          f"phase={ci.get('phase')!r} r2={ci.get('r_squared')}")
n_indexed = sum(1 for m in members if (m.get("snapshot") or {}).get("confirmed_index"))
assert len(members) >= 2, "FAIL: need >=2 members for heatmap rows"
assert n_indexed >= 2, "FAIL: need >=2 confirmed indices"
assert ov, "FAIL: ordering_variable must be set for the rotated axis label"
print(f"[seed] OK: {len(members)} rows, {n_indexed} indexed, ordering_variable set")
print(f"[seed] series-builder URL path: /series/{sid}")
PY
echo "[seed] done. SERIES_ID=$SERIES_ID"
