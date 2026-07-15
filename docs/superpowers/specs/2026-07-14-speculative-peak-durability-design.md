# Speculative-index peak durability — design

**Date:** 2026-07-14
**Status:** Draft for review
**Prior investigation:** live repro + prod-DB-copy inspection (this session). Symptoms
reported: focus-page hover does nothing useful for custom (speculative) indices, comb
coloring "off", Miller panel draws no matched peaks.

## Problem

Speculative ("custom") indices lose their peak assignments permanently, and two
independent numeric/rendering bugs degrade what remains.

### 1. Peak-assignment data loss (the reported symptom)

`_reingest_inner!` (`packages/HimalayaUI/src/pipeline.jl`) handles speculative indices
across a reanalysis like this:

1. Snapshot each speculative's `index_peaks` rows **in memory** (q-values via a join to
   the live `auto_peaks` / `peak_curations` rows).
2. Wipe `index_peaks` for all speculatives (required: auto peaks are deleted and
   re-detected with new ids, so old FK refs would dangle).
3. Re-resolve each snapshotted q against the new peaks within `SNAP_TOL`; run a
   basis-driven auto-discovery pass for unfilled ratio positions.
4. If fewer than 2 peaks re-resolve → mark the index `stale` and **commit the wipe**.

Step 4 is destructive: the user's ratio→peak assignments are gone forever. The next
reingest has nothing to re-attach (empty snapshot → stale again). A deleted curation
peak has the same effect on its assignment (the join yields NULL q, the row can never
re-resolve). The `< 2` threshold also strands legitimately-created 1-peak
(anchor-only) speculatives.

Observed in prod (`/opt/Himalaya.jl/data/himalaya.db`, read-only copy): 18 speculative
indices, **11 with zero `index_peaks` rows** and NULL `score`/`r_squared` (5 `stale`,
6 still `candidate`). Empty `peaks` renders exactly the reported symptoms: no Miller
dots, every comb tick permanently in the dim "unmatched" style, hover highlight faint
and useless, card shows "0 peaks / score —".

### 2. Basis-scale corruption for cubic phases (latent, confirmed in source)

Auto indices store the **normalized** basis (q of the first ratio position; core
convention: `predictpeaks = basis × phaseratios(P; normalize=true)`). Speculative
creation (`build_speculative_index`, `speculative.jl`) fits basis against
**un-normalized** ratios, and the pipeline re-attach does the same (`new_basis`,
plus the mixed-scale `basis_for_snap × ratios_normed` discovery predictions).
Every consumer assumes the normalized convention (`predicted_q_for_phase`,
`routes_analysis.jl`). For cubic phases (first raw ratio ≠ 1) a saved speculative's
served `predicted_q` is shrunk by √2 (Pn3m, Im3m), √3 (Fm3m, Fd3m), or √6 (Ia3d) —
and the builder's snap preview (`compute_snap`, already normalized) looks correct
until save. All prod speculatives are Square/Lamellar/Hexagonal (first ratio = 1),
which is why this hasn't been seen; it also means basis-driven healing would search
at the wrong q for cubics, so this must be fixed before the repair pass.

### 3. NaN confidence band in MillerPlot (cosmetic)

`Plot.linearRegressionY` draws a confidence band; with exactly 2 points the band
math is NaN and the browser logs an SVG `<path d="M…,NaN…">` error on every render.
Common for speculatives (2 assigned peaks).

## Goals

- A reanalysis can never permanently destroy a speculative index's peak assignments.
- One stored meaning for `indices.basis` regardless of `kind` (the normalized
  convention auto indices already use).
- Existing peak-less prod indices heal on a reingest under the fixed code.
- A peak-less speculative is legible in the UI, not a silent "0 peaks".

## Non-goals

- Editing a speculative's assignments after creation (no such UI today).
- Predicted-only open circles in the Miller panel for peak-less indices (add later if
  healed data still leaves gaps).
- Any change to auto-index behavior, scoring, or the event log/mutation queue.

## Design

### A. Durable intent: `speculative_peak_intents` table

New table, written at creation, never touched by reingest wipes:

```sql
CREATE TABLE speculative_peak_intents (
    index_id       INTEGER NOT NULL REFERENCES indices(id) ON DELETE CASCADE,
    ratio_position INTEGER NOT NULL,
    q              REAL    NOT NULL,   -- q of the assigned peak at creation time
    PRIMARY KEY (index_id, ratio_position)
);
```

- `insert_speculative_index!` writes one intent row per assignment (anchor +
  additional), using the assigned peak's q.
- `_reingest_inner!` re-attach reads intents (not a join against live peak tables) as
  the snapshot source. `index_peaks` remains the *resolved view*, rebuilt each
  reingest as today.
- Basis-driven auto-discovery results do **not** become intents (machine guesses, not
  user intent); they land only in `index_peaks`, as today.
- Migration backfills intents from existing `index_peaks` rows (best available proxy)
  for speculatives that still have them. The 11 peak-less prod rows have no intent to
  backfill; they heal via the stored-basis discovery path (D).

Rejected alternative — `q_snapshot` column on `index_peaks` with nullable `peak_id`:
keeps one table but changes the API row shape (`peak_id: number | null`) and the
meaning of `peaks` for every consumer. The intent table leaves the API contract and
all frontend data paths untouched.

### B. Non-destructive re-attach

In `_reingest_inner!`'s speculative loop:

- Resolve from intents + discovery as today.
- **0 peaks resolved** → mark `stale`, leave `index_peaks` empty, *intents survive* —
  the next reingest retries from the same intents. Score/r²/lattice/basis keep their
  last known values (they describe the intent-fit, still the best available summary).
- **≥ 1 peak resolved** → write resolved rows, recompute basis/score/r²/lattice_d,
  set `status = 'candidate'` (the existing successful-path UPDATE already resets
  status). This replaces the `< 2 → stale` threshold: 1-peak speculatives are legal
  at creation and stay legal across reingest.

### C. One basis convention (normalized) everywhere

- `build_speculative_index`: fit against `phaseratios(P; normalize=true)`; residual
  ideals likewise. The `Himalaya.Index{P}` it constructs then carries the same basis
  meaning as auto indices (`fit`/`score` are unaffected — `fit` refits d internally
  from peaks; `score` never reads basis).
- Pipeline re-attach: `basis_for_snap` fits and `new_basis`/residuals switch to
  normalized ratios, eliminating the mixed-scale discovery predictions.
- One-shot sentinel migration (existing `schema_migrations` pattern in `db.jl`):
  `UPDATE indices SET basis = basis * first(phaseratios(P))` for speculative rows of
  cubic phases (factor 1 elsewhere — no-op for all current prod rows). Runs on
  `open_db` before any new-code write, so old-scale and new-scale rows can't mix.

### D. Healing existing peak-less rows

No new mechanism needed: with C in place, the existing stored-basis fallback +
auto-discovery in the re-attach loop re-attaches peaks for indices that have a basis
but no rows, and the success path resets `status` to `candidate` and fills
score/r². Ops runbook (prod):

1. Deploy fixed build (rebuild sysimage; prod is currently running an older build —
   schema shows pre-current-main vintage).
2. Back up the prod DB file.
3. `bin/himalaya reingest <experiment-dir>` per experiment (migrations apply on DB
   open). Verify: `SELECT COUNT(*) FROM indices i WHERE kind='speculative' AND NOT
   EXISTS (SELECT 1 FROM index_peaks ip WHERE ip.index_id=i.id)` → only genuinely
   unresolvable indices remain, and those now show the UI chip (E).

### E. UI honesty + NaN band

- `IndexCard` (`PhasePanel.tsx`): when `kind === "speculative" && peaks.length === 0`,
  render a warning chip — `peaks unresolved — reanalyze` — in place of the silent
  "0 peaks" pill styling. No new component; it's a conditional on the existing pill.
- `MillerPlot.tsx`: disable the regression confidence band for groups with fewer than
  3 points (Plot's `ci` option), killing the NaN `<path>` console error.

## Error handling

- Re-attach never throws on unresolvable intents; worst case is `stale` + empty
  resolved view + surviving intents.
- Migration is sentinel-gated and idempotent; a DB missing the intents table gets it
  created via `CREATE TABLE IF NOT EXISTS` in `migrate_schema!`.
- `DELETE /api/indices/:id` needs no change: `ON DELETE CASCADE` clears intents.

## Testing

Julia (`packages/HimalayaUI/test/`):

- Basis convention: create a cubic (Pn3m) speculative; assert
  `predicted_q_for_phase(phase, stored_basis)[anchor_rpos] ≈ anchor q` (fails before
  C, passes after).
- Migration: seed an old-scale cubic speculative row; open DB; assert basis rescaled
  exactly once (sentinel).
- Re-attach round-trip: create speculative → reingest with unchanged data → intents
  and resolved rows survive; with all peaks shifted out of tolerance → `stale`,
  `index_peaks` empty, **intents intact** → shift peaks back, reingest → healed to
  `candidate`.
- Heal-from-basis: seed a peak-less speculative (no intents, basis set) matching the
  prod shape → reingest → rows re-attached, status `candidate`, score non-NULL.
- 1-peak survival: anchor-only speculative stays `candidate` across reingest.

Vitest (`packages/HimalayaUI/frontend/test/`):

- IndexCard: speculative with `peaks: []` renders the unresolved chip; non-empty does
  not.
- MillerPlot: 2-point speculative renders no `<path>` with NaN in its `d` attribute.

Regression floors, not exact counts, where fixtures are real-data-derived (per repo
convention).

## Rollout

Single PR is fine (backend + frontend are small and coupled only through the healed
data). Prod redeploy + reingest per runbook in D. No API shape change; no
frontend-cache or queue implications (`reingest` is CLI-side, not an event-log
mutation).
