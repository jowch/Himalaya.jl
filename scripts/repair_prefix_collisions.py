#!/usr/bin/env python3
"""
One-shot DB repair for the JC_C04 / JC_C04s prefix-collision bug:
exposures whose filename starts with a longer sample's prefix got
duplicated under the shorter sample. The duplicates carry the same
filename as the legit rows under the longer-prefix sample, so we can
safely delete them after the boundary-aware matcher fix lands.

Usage:
  python3 repair_prefix_collisions.py /path/to/himalaya.db          # dry-run (default)
  python3 repair_prefix_collisions.py /path/to/himalaya.db --commit  # actually delete

Strategy: for each filename that exists under >1 sample in the same
experiment, KEEP the row attached to the sample with the fewer total
exposures (the more-specific match: `C4 ⭐️` has 8, `C4` has 16+8=24
→ keep C4⭐️). Delete all other rows + their children (peaks, indices,
index_peaks, index_groups, index_group_members, exposure_tags) in a
single transaction.
"""

import argparse
import sqlite3
import sys
from collections import defaultdict


def find_duplicates(cur):
    """Return dict[filename] -> list of (exposure_id, sample_id, sample_label)
    for filenames that appear under >1 sample in the same experiment."""
    cur.execute("""
        SELECT e.filename, e.id, e.sample_id, s.label, s.experiment_id
        FROM exposures e JOIN samples s ON s.id = e.sample_id
        WHERE e.filename IN (
            SELECT e2.filename FROM exposures e2
            JOIN samples s2 ON s2.id = e2.sample_id
            GROUP BY e2.filename, s2.experiment_id
            HAVING COUNT(*) > 1
        )
        ORDER BY e.filename, e.sample_id
    """)
    out = defaultdict(list)
    for fname, eid, sid, label, exp_id in cur.fetchall():
        out[(exp_id, fname)].append((eid, sid, label))
    return out


def sample_exposure_count(cur, sid):
    cur.execute("SELECT COUNT(*) FROM exposures WHERE sample_id = ?", (sid,))
    return cur.fetchone()[0]


def pick_keeper(cur, candidates):
    """Pick the (eid, sid, label) with the fewer-exposure sample (more specific)."""
    return min(candidates, key=lambda c: sample_exposure_count(cur, c[1]))


def collect_doomed_ids(duplicates, cur):
    """Return list of exposure_ids to delete (one per (sample, file) collision)."""
    doomed = []
    for (exp_id, fname), candidates in duplicates.items():
        keeper = pick_keeper(cur, candidates)
        for c in candidates:
            if c == keeper:
                continue
            doomed.append((c[0], c[1], c[2], fname, keeper))
    return doomed


def child_counts(cur, eid):
    """Per-exposure child row counts for reporting."""
    counts = {}
    counts["peaks"] = cur.execute(
        "SELECT COUNT(*) FROM peaks WHERE exposure_id = ?", (eid,)).fetchone()[0]
    counts["indices"] = cur.execute(
        "SELECT COUNT(*) FROM indices WHERE exposure_id = ?", (eid,)).fetchone()[0]
    counts["index_groups"] = cur.execute(
        "SELECT COUNT(*) FROM index_groups WHERE exposure_id = ?", (eid,)).fetchone()[0]
    counts["exposure_tags"] = cur.execute(
        "SELECT COUNT(*) FROM exposure_tags WHERE exposure_id = ?", (eid,)).fetchone()[0]
    return counts


def cascade_delete(cur, eid):
    """Delete an exposure and all rows that reference it.
    Order matters: deepest children first, since FKs are NO ACTION ON DELETE."""
    # index_peaks references both indices and peaks of this exposure
    cur.execute("""
        DELETE FROM index_peaks WHERE peak_id IN (
            SELECT id FROM peaks WHERE exposure_id = ?
        ) OR index_id IN (
            SELECT id FROM indices WHERE exposure_id = ?
        )
    """, (eid, eid))
    # index_group_members references indices and groups of this exposure
    cur.execute("""
        DELETE FROM index_group_members WHERE index_id IN (
            SELECT id FROM indices WHERE exposure_id = ?
        ) OR group_id IN (
            SELECT id FROM index_groups WHERE exposure_id = ?
        )
    """, (eid, eid))
    cur.execute("DELETE FROM index_groups   WHERE exposure_id = ?", (eid,))
    cur.execute("DELETE FROM indices        WHERE exposure_id = ?", (eid,))
    cur.execute("DELETE FROM peaks          WHERE exposure_id = ?", (eid,))
    cur.execute("DELETE FROM exposure_tags  WHERE exposure_id = ?", (eid,))
    cur.execute("DELETE FROM exposure_sources WHERE source_exposure_id = ? OR averaged_exposure_id = ?",
                (eid, eid))
    cur.execute("DELETE FROM exposures      WHERE id = ?", (eid,))


def main():
    p = argparse.ArgumentParser()
    p.add_argument("db", help="Path to himalaya.db")
    p.add_argument("--commit", action="store_true", help="Actually delete (default: dry-run)")
    args = p.parse_args()

    db = sqlite3.connect(args.db)
    db.execute("PRAGMA foreign_keys = ON")
    cur = db.cursor()

    duplicates = find_duplicates(cur)
    if not duplicates:
        print("No duplicate exposures found — nothing to repair.")
        return

    doomed = collect_doomed_ids(duplicates, cur)

    print(f"Found {len(duplicates)} duplicated filenames across "
          f"{len({d[1] for d in doomed})} sample(s).")
    print(f"Will delete {len(doomed)} duplicate exposure row(s) "
          f"and their cascading children.\n")

    total_children = {"peaks": 0, "indices": 0, "index_groups": 0, "exposure_tags": 0}
    for eid, sid, label, fname, keeper in doomed:
        c = child_counts(cur, eid)
        for k in total_children:
            total_children[k] += c[k]
        keeper_eid, keeper_sid, keeper_label = keeper
        print(f"  delete exposure_id={eid:<5} (sample={sid} {label!r}, fname={fname})")
        print(f"      → keeping  exposure_id={keeper_eid:<5} (sample={keeper_sid} {keeper_label!r})")
        print(f"      → cascade: peaks={c['peaks']}, indices={c['indices']}, "
              f"index_groups={c['index_groups']}, exposure_tags={c['exposure_tags']}")

    print(f"\nCascade totals: {total_children}")

    if not args.commit:
        print("\nDRY RUN — re-run with --commit to apply.")
        return

    print("\nCommitting...")
    cur.execute("BEGIN")
    try:
        for eid, *_ in doomed:
            cascade_delete(cur, eid)
        db.commit()
        print(f"Deleted {len(doomed)} exposure rows and their children.")
    except Exception as e:
        db.rollback()
        print(f"ROLLED BACK: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
