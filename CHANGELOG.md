# Changelog

Covers core `Himalaya` and the `HimalayaUI` / `HimalayaDB` packages built on it.
Core releases are headed `Himalaya <version>`.

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
  `HimalayaUI` has no `[sources]`, so without the bound a resolve could supply a
  core whose `phaseratios` differs from the one the calling code was written
  against — wrong reflections under the same version string, with no error. The
  bound makes that a resolution failure instead. When working on core itself,
  `Pkg.develop` the local checkout (see `CLAUDE.md`): the bound separates
  *versions*, so an unbumped edit to `phaseratios` still diverges from whatever
  is published.
- **`HimalayaDB.connect` now warns when the database predates this migration.**
  Reading a pre-migration database — an old backup, or a deploy that has not been
  restarted — makes `reconstruct_index` return Hexagonal peaks one reflection too
  high for every position 6 and above, which `Himalaya.fit` turns into a **wrong
  lattice constant with no error**. The package opens `query_only` and cannot
  migrate the database itself, so a warning at connect time is the only available
  signal; it fires only when the sentinel is absent *and* Hexagonal indices exist.
  It warns rather than throws, because reading an old backup is legitimate and
  every other phase is unaffected. Treat Hexagonal results from a warned
  connection as invalid.

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

### Upgrade procedure

Steps 1 and 2 are ordered deliberately — the migration runs on the **first**
`open_db` after the new version is in place, which is whatever command you run
first, not a separate step you opt into.

1. **Back up the database before starting the upgraded build.** By the time you
   can run anything else, `hex_sqrt11_removal_v1` has already deleted Hexagonal
   `index_peaks` rows and renumbered `speculative_peak_intents` — and that
   renumber is **not idempotent**, so it cannot be replayed or reversed in place.
   Restoring the backup is the only way back.
2. **Run `bin/himalaya analyze --all`.** The migration only deletes rows and
   reopens the memoization gate; until a re-analysis runs, Hexagonal indices
   display as claiming no peaks and carry stale scores. Trace files must be
   reachable at the paths each experiment's config resolves.
   - **Rejected exposures are skipped** by `analyze` and therefore never
     rebuilt. A rejected exposure holding a Hexagonal index keeps that index
     with no claimed peaks until someone un-rejects it and re-analyses. Nothing
     unique is lost — auto rows are derived and speculative intents survive —
     but the claimless state persists rather than healing on its own.
3. **Expect a curated phase call to be withdrawn.** A Hexagonal assignment whose
   evidence included the `√11` can fall below `minpeaks` and stop being
   reproducible; the re-analysis then drops that assignment. `analyze` reports
   these per-exposure and as a closing summary — read that summary, because it
   is a human decision disappearing.
