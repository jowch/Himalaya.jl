# Changelog

Covers both the core `Himalaya` package (versioned in the root `Project.toml`)
and the `HimalayaUI` application. Core releases are headed `Himalaya <version>`;
entries under `Unreleased` are application-level unless stated otherwise.

## Unreleased

### Breaking changes

- **`samples.label` removed; replaced by `name` (stable identifier, was column 3) and
  `display_name` (editable label, was column 2). Manifest column meanings swapped.**
  Existing experiments need `himalaya migrate-toml <experiment-dir>` to upgrade
  their `experiment.toml`. Existing DBs auto-migrate on first `open_db` after
  deploy. Issue #88.
- **`/api/export` CSV header changed: `sample_label` → `sample_display_name`.**
  Downstream pipelines parsing this CSV need to update their column names.
- **`PATCH /api/samples/:id` no longer accepts `name`; use `display_name`.**
  `samples.name` is now the stable identifier and is set only at ingest time.
- **`PATCH /api/experiments/:id` no longer accepts any field.** The route is
  retained as defensive surface for future fields. Experiment renames must go
  through `experiment.toml` + reingest.
- **First boot after migration purges `idempotent_responses`.** In-flight
  `client_op_id` retries from before the deploy will get fresh executions
  rather than cached pre-rename response bodies.

## Himalaya 0.6.0 — 2026-07-30

Minor, not patch: `phaseratios` output changes, which renumbers persisted data.
Under the pre-1.0 convention the minor is the breaking slot, and a patch release
would let an existing `Himalaya = "0.5"` bound silently resolve the new physics —
the exact mismatch this release guards against.

### Breaking changes

- **`phaseratios(Hexagonal)` no longer contains `√11`** (issue #304, PR #307).
  `N = 11` is not a Loeschian number (`N = h² + hk + k²` has no integer solution),
  so it is not a permitted 2D hexagonal reflection. The series goes from 14 to 13
  entries and **position 6 is now `√12`**, so every `ratio_position` past 5
  renumbers. Downstream code that indexes into the series by position must be
  re-checked; matching on ratio *values* is stable across this change.
- **All Hexagonal `score` values shift.** `score`'s coverage denominator is
  `sum(1/r for r in 1:length(phaseratios(P)))`, so it moves from `H(14)` to
  `H(13)` — roughly +2.25% for an unchanged found-set, more where positions
  renumber, and *down* for an index that claimed the removed `√11`. Because
  `auto_group` and `remove_subsets` order on score, cross-phase orderings can
  flip: an index that previously lost a subset comparison may now survive.
- **Consumers must resolve `Himalaya >= 0.6`.** `HimalayaUI` and `HimalayaDB`
  now declare `[compat] Himalaya = "0.6"`. `Manifest.toml` is gitignored and
  `HimalayaUI` has no `[sources]`, so a bare `Pkg.instantiate()` previously
  resolved a registry core with different physics under the same version string.
  It now fails to resolve instead — `Pkg.develop` the local core first (see
  `CLAUDE.md`). This is deliberate: a load failure is preferable to silently
  wrong reflections.

### Migrations

- **`hex_sqrt11_removal_v1`** runs automatically on the next `open_db`. It
  deletes Hexagonal `index_peaks` rows (derived state — `persist_analysis!`
  rebuilds them), renumbers Hexagonal `speculative_peak_intents`
  (position 6 dropped, 7+ shifted down one), and clears
  `analysis_inputs_hash` for every exposure holding at least one index so the
  next analysis re-indexes and re-scores. Exposures with no indices are left
  alone — they are provably unaffected.
- The migration **defers without recording its sentinel** if the loaded core
  still lists `√11`, so a stale environment cannot renumber durable data against
  a series it isn't using.

### Action required after upgrading

- **Run `bin/himalaya analyze --all`.** The migration only deletes rows and
  reopens the memoization gate; until a re-analysis runs, Hexagonal indices
  display as claiming no peaks and carry stale scores. Trace `.dat` files must
  be reachable. Take a database backup first.
- **A curated phase call may be withdrawn.** A Hexagonal assignment whose
  evidence included the `√11` can fall below `minpeaks` and stop being
  reproducible; the re-analysis then drops that assignment. `analyze` now
  reports these per-exposure and as a closing summary.
