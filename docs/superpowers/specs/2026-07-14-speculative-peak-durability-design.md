# Speculative-index peak durability — design

**Date:** 2026-07-14
**Status:** Draft for review (revision 2 — after 4-reviewer fan-out: saxs-physics,
himalaya, frontend, queue)
**Prior investigation:** live repro + prod-DB-copy inspection (this session). Symptoms
reported: focus-page hover does nothing useful for custom (speculative) indices, comb
coloring "off", Miller panel draws no matched peaks.

## Problem

Speculative ("custom") indices lose their peak assignments permanently, and two
independent numeric/rendering bugs degrade what remains.

### 1. Peak-assignment data loss (the reported symptom)

`_persist_analysis_inner!` (`packages/HimalayaUI/src/pipeline.jl:247`) handles
speculative indices across a reanalysis like this:

1. Snapshot each speculative's `index_peaks` rows **in memory** (q-values via a join to
   the live `auto_peaks` / `peak_curations` rows; `pipeline.jl:273-285`).
2. Wipe `index_peaks` for all speculatives (`:305-308`; required: auto peaks are
   deleted and re-detected with new ids, so old FK refs would dangle).
3. Re-resolve each snapshotted q against the new peaks within `SNAP_TOL`; run a
   basis-driven auto-discovery pass for unfilled ratio positions (`:450-476`).
4. If fewer than 2 peaks re-resolve → mark the index `stale` and **commit the wipe**
   (`:478-482`).

This loop is reached from `analyze_exposure!` — i.e. from CLI `analyze`/`init` **and**
from server routes on the SSE hot path (`POST /api/exposures/{id}/analyze`,
`POST /api/experiments/{id}/analyze`, `POST /api/experiments/{id}/reingest`, and every
peak-curation mutation in `routes_peaks.jl`). The wipe+re-add is already atomic per
exposure: `persist_analysis!` wraps the inner function in `SQLite.transaction` under
`_DB_WRITE_LOCK` (`pipeline.jl:229-244`).

Step 4 is destructive: the user's ratio→peak assignments are gone forever. The next
analysis has nothing to re-attach (empty snapshot → empty again). A deleted curation
peak has the same effect on its assignment (the join yields NULL q, the row can never
re-resolve). The `< 2` threshold also strands legitimately-created 1-peak
(anchor-only) speculatives. The `stale` status it writes is a dead enum value: an
unconditional R3.2 normalization flips `stale → candidate` on every `open_db`
(`db.jl:346`) — hash mismatch, not status, is the real staleness signal.

Observed in prod (`/opt/Himalaya.jl/data/himalaya.db`, read-only copy): 18 speculative
indices, **11 with zero `index_peaks` rows** and NULL `score`/`r_squared`. Empty
`peaks` renders exactly the reported symptoms: no Miller dots, every comb tick
permanently in the dim "unmatched" style, hover highlight faint and useless, card
shows "0 peaks / score —".

**Why the broken rows also read as "fresh":** `_persist_analysis_inner!` syncs every
index's `inputs_hash` to the exposure hash at commit (`pipeline.jl:592-597`), and
`analyze_exposure!` skips the whole persist step — including the re-attach loop —
whenever `stored_inputs_hash == new_inputs_hash && indices_count > 0`
(`pipeline.jl:846`, `:916`). For the 11 prod rows the data hasn't changed since the
destructive wipe, so plain reanalysis is a no-op: no heal, no `analyze_run` frame.
Any fix must reopen this gate deliberately (§D).

### 2. Basis-scale corruption for cubic phases (latent, confirmed in source)

Auto indices store the **normalized** basis (q of the first ratio position; core
convention: `predictpeaks = basis × phaseratios(P; normalize=true)`,
`src/index.jl:33`). Speculative creation (`build_speculative_index`,
`speculative.jl:121-136`) fits basis against **un-normalized** ratios, and the
pipeline re-attach does the same (`new_basis` at `pipeline.jl:491`, plus the
mixed-scale `basis_for_snap × ratios_normed` discovery predictions at `:454`).
Every consumer assumes the normalized convention (`predicted_q_for_phase`,
`routes_analysis.jl:64-75`). For cubic phases (first raw ratio ≠ 1) a saved
speculative's served `predicted_q` is shrunk by √2 (Pn3m, Im3m), √3 (Fm3m, Fd3m), or
√6 (Ia3d) — and the builder's snap preview (`compute_snap`, already normalized) looks
correct until save. All prod speculatives are Square/Lamellar/Hexagonal (first ratio
= 1), which is why this hasn't been seen; it also means basis-driven healing would
search at the wrong q for cubics, so this must be fixed before the repair pass.

`Himalaya.fit` and `Himalaya.score` are unaffected by which scale the stored basis
uses: `fit` regresses d internally from peaks and never reads `index.basis`
(`src/index.jl:230-245`); `score` reads only peaks + sharpness (`:265-280`).
Residuals are scale-invariant too — `ratios_unnorm[rpos]·basis_unnorm ≡
ratios_normed[rpos]·basis_norm` (the LSQ slopes differ by exactly `first_ratio`) — so
**`indices.basis` is the only column that needs migrating**.

### 3. NaN confidence band in MillerPlot (cosmetic)

`Plot.linearRegressionY` draws a confidence band; with exactly 2 points the band
math divides by zero residual DOF and the browser logs an SVG `<path d="M…,NaN…">`
error on every render. Common for speculatives (2 assigned peaks). n=1 groups are
already skipped (`MillerPlot.tsx:78`); n≥3 is finite — the failing case is exactly
n=2.

## Goals

- An analysis run can never permanently destroy a speculative index's peak
  assignments.
- One stored meaning for `indices.basis` regardless of `kind` (the normalized
  convention auto indices already use).
- Existing peak-less prod indices heal on one deliberate reanalysis pass under the
  fixed code.
- A peak-less speculative is legible in the UI, not a silent "0 peaks".

## Non-goals

- Editing a speculative's assignments after creation (no such UI today).
- Predicted-only open circles in the Miller panel for peak-less indices (add later if
  healed data still leaves gaps).
- Any change to auto-index behavior, scoring, `auto_group`/`remove_subsets` (verified:
  speculatives never flow through them — `pipeline.jl:267` filters `kind='auto'`;
  speculative group membership is preserved by stable id), or the event-log/mutation-
  queue *contracts* (no new event kinds, no API shape change — but see §Rollout for
  the SSE surface this change is visible on).

## Design

### A. Durable intent: `speculative_peak_intents` table

New table, written at creation, never touched by analysis wipes:

```sql
CREATE TABLE speculative_peak_intents (
    index_id       INTEGER NOT NULL REFERENCES indices(id) ON DELETE CASCADE,
    ratio_position INTEGER NOT NULL,
    q              REAL    NOT NULL,   -- q of the assigned peak at creation time
    PRIMARY KEY (index_id, ratio_position)
);
```

- `insert_speculative_index!` writes one intent row per assignment (anchor +
  additional), using the assigned peak's q. The intent INSERTs sit **after all
  validation** (`build_speculative_index` + `_kind_for` throw before any write today),
  adjacent to the existing `index_peaks` INSERTs, inside the same `with_idempotency`
  outer transaction — an idempotent replay returns the cached response without
  re-executing, so no double-write window exists.
- The re-attach loop reads intents (not a join against live peak tables) as the
  snapshot source. `index_peaks` remains the *resolved view*, rebuilt each analysis
  as today.
- Basis-driven auto-discovery results do **not** become intents (machine guesses, not
  user intent); they land only in `index_peaks`, as today.
- Intent rows are frozen at creation. `DELETE /api/indices/:id` needs no change: the
  route deletes the `indices` row directly inside the event transaction and
  `PRAGMA foreign_keys = ON` (`db.jl:1780`) makes the CASCADE fire with it.
- Migration backfills intents from existing `index_peaks` rows (best available proxy)
  for speculatives that still have them, skipping rows whose join q is NULL (peak
  already deleted — those assignments were orphaned before this fix existed and
  cannot be resurrected). The 11 peak-less prod rows have no intent to backfill; they
  heal via the stored-basis discovery path (§D).

Rejected alternative — `q_snapshot` column on `index_peaks` with nullable `peak_id`:
keeps one table but changes the API row shape (`peak_id: number | null`) and the
meaning of `peaks` for every consumer. The intent table leaves the API contract and
all frontend data paths untouched.

### B. Non-destructive re-attach

In `_persist_analysis_inner!`'s speculative loop:

- Resolve from intents + discovery as today.
- **0 peaks resolved** → leave `index_peaks` empty and *do nothing else*: intents
  survive (next analysis retries from them), and score/r²/lattice/basis keep their
  last known values — the current failure path already does no UPDATE beyond the
  status write, which we drop. No `stale` marking: the status value is dead
  (`db.jl:346` normalizes it away on every open); the UI signal for this state is
  the §E chip, keyed on `peaks.length === 0`.
- **≥ 1 peak resolved** → write resolved rows, recompute basis/score/r²/lattice_d
  via the existing success-path UPDATE (which already sets `status = 'candidate'`).
  This replaces the `< 2 → stale` threshold: 1-peak speculatives are legal at
  creation and stay legal across analysis. Note: a 1-peak fit has undefined R²
  (RSS=0/TSS=0 → NaN; SQLite binds NaN as NULL) — this is pre-existing behavior for
  1-peak creation, and both serializers already NULL-guard `r_squared`
  (`pipeline.jl:645`), so it is benign. `score` is well-defined for 1 peak.

### C. One basis convention (normalized) everywhere

- `build_speculative_index`: fit against `phaseratios(P; normalize=true)`; residual
  ideals likewise (numerically identical, per §Problem 2).
- Pipeline re-attach: `basis_for_snap` fits and `new_basis`/residuals switch to
  normalized ratios, eliminating the mixed-scale discovery predictions.
- One-shot sentinel migration (existing `schema_migrations` pattern):
  `UPDATE indices SET basis = basis * first(phaseratios(P))` for speculative rows of
  cubic phases (factor 1 elsewhere — no-op for all current prod rows).

**Migration placement** (ordering hazards are real; follow the `migrate_compare!` /
`migrate_series!` precedent):

- The intents table is **not** added to the `SCHEMA` const (`create_schema!` runs
  before the PK/autoincrement rebuild, and `indices` is rename-rebuilt by
  `migrate_pk_to_autoincrement!` on legacy DBs, which would rewrite the FK ref).
  It gets its own migration function placed after `migrate_pk_to_autoincrement!` and
  the FK-heal migrations.
- The intents **backfill** runs after `migrate_r2_split_peaks!` (its q values come
  from the `COALESCE(ap.q, pc.q)` join, which needs the post-R2 peak tables).
- The basis-rescale sentinel and the backfill sentinel run after `migrate_series!`
  (which creates `schema_migrations`). `open_db` runs `migrate_schema!` before
  returning and every writer goes through `open_db`, so "migrated before any
  new-code write" holds.
- All migration writes run inside `migrate_schema!`'s existing execution context on
  DB open; creation-time intent writes run inside the POST route's existing
  transaction (§A). No new transaction boundaries anywhere.

### D. Healing existing peak-less rows

The re-attach mechanics need nothing new: with C in place, the existing stored-basis
fallback + auto-discovery (`pipeline.jl:437-476`) re-attaches peaks for indices that
have a basis but no rows, and the success path resets `status`/score/r².

What **does** need a trigger is the gate: reanalysis is memoized shut for exactly
these rows (§Problem 1 — `indexpeaks_skipped`). The migration therefore also runs:

```sql
UPDATE exposures SET analysis_inputs_hash = NULL
 WHERE id IN (SELECT DISTINCT exposure_id FROM indices i
               WHERE i.kind = 'speculative'
                 AND NOT EXISTS (SELECT 1 FROM index_peaks ip
                                  WHERE ip.index_id = i.id));
```

A NULL stored hash can never equal the computed hash, so the next analysis of those
exposures takes the full path (peak re-detection under the current build included —
inherent to healing) and the re-attach loop actually runs. Exposures with healthy
speculatives are untouched. This is a one-shot repair, not a permanent gate bypass:
once intents + §B exist, this loss state can't be recreated, so no standing
`indexpeaks_skipped` carve-out is added.

Ops runbook (prod):

1. Deploy fixed build (rebuild sysimage; prod currently runs an older build —
   schema shows pre-current-main vintage).
2. Back up the prod DB file.
3. `bin/himalaya analyze <experiment-dir>` per experiment (**not** `reingest` —
   reingest only re-syncs manifest metadata and never analyzes). Migrations apply on
   DB open; the hash-NULL from the migration makes analyze actually process the
   affected exposures. The UI analyze-all route works too.
4. Verify: `SELECT COUNT(*) FROM indices i WHERE kind='speculative' AND NOT EXISTS
   (SELECT 1 FROM index_peaks ip WHERE ip.index_id=i.id)` → only genuinely
   unresolvable indices remain, and those now show the UI chip (§E).

### E. UI honesty + NaN band

- `IndexCard` (`PhasePanel.tsx`): when `kind === "speculative" && peaks.length === 0`,
  swap the "{n} peaks" pill for a warning chip — copy: **"peaks unresolved"** (no
  "reanalyze" CTA: after the one-shot heal, reanalysis of an unchanged exposure is
  correctly a no-op, so the CTA would be a dead end). Stable selectors per repo
  convention: `data-unresolved` attribute on the `<li>` (parallel to
  `data-speculative` / `data-low-r2`) and a `data-testid` on the chip. Color from an
  existing `@theme` token (e.g. the `text-error` token already used in this file) —
  no new hex values.
- `MillerPlot.tsx`: pass `ci: 0` to `Plot.linearRegressionY` for regression groups
  with fewer than 3 points. Verified against installed `@observablehq/plot@0.6.17`:
  `ci: 0` hides only the band; the (dashed, speculative-cue) regression line renders
  unconditionally.

## Error handling

- Re-attach never throws on unresolvable intents; worst case is empty resolved view
  + surviving intents.
- Migrations are sentinel-gated and idempotent; `CREATE TABLE IF NOT EXISTS` for the
  intents table.
- A 4xx returned (not thrown) from the speculative POST body commits its writes
  (`with_idempotency` rolls back only on throw) — unchanged from today, and safe
  because all validation precedes the first INSERT; the intents write keeps that
  ordering (§A).

## Testing

Julia (`packages/HimalayaUI/test/`) — drive `analyze_exposure!`, not `reingest!`
(reingest never analyzes); follow the `test_pipeline.jl` curation-lifecycle pattern:

- Basis convention: create a **1-peak (anchor-only)** cubic (Pn3m) speculative;
  assert `predicted_q_for_phase(phase, stored_basis)[anchor_rpos] ≈ anchor q`
  (exact only for 1-peak fixtures — a multi-peak LSQ basis deviates by the fit
  residual; use a residual-scale tolerance if testing multi-peak).
- Migration: seed an old-scale cubic speculative row; open DB; assert basis rescaled
  exactly once (sentinel), and `analysis_inputs_hash` NULLed only for exposures
  owning peak-less speculatives.
- Re-attach round-trip: create speculative → `analyze_exposure!` with unchanged data
  → hash-gate no-op, intents and resolved rows survive; force a full re-analysis with
  all peaks shifted out of tolerance → `index_peaks` empty, **intents intact**,
  status still `candidate` → shift peaks back, analyze → healed with resolved rows
  and score set.
- Heal-from-basis: seed a peak-less speculative (no intents, basis set, NULL exposure
  hash) matching the prod shape → `analyze_exposure!` → rows re-attached, score
  non-NULL.
- 1-peak survival: anchor-only speculative keeps its row and (NULL) r² across
  analysis; serializer emits `r_squared: null`.
- SSE surface (contract layer, per `docs/contract-testing.md`): after a healing
  analyze, the `analyze_run` `post_state` payload
  (`_serialized_indices_for_broadcast`) contains the speculative with non-empty
  `peaks` — the existing plumbing carries healed state to clients; this test pins
  that.

Vitest (`packages/HimalayaUI/frontend/test/`):

- IndexCard: speculative with `peaks: []` renders the unresolved chip (assert via
  `data-unresolved` / chip testid); non-empty does not.
- MillerPlot: assert the mocked `linearRegressionY` is called with `ci: 0` for a
  2-point speculative group and without it (or `ci` unset) for a 3-point group. (The
  existing test file mocks `@observablehq/plot` wholesale, so a "no NaN path in the
  DOM" assertion would pass vacuously — assert on the mark options instead, matching
  the existing `"invokes Plot.dot and Plot.linearRegressionY"` test.)

Regression floors, not exact counts, where fixtures are real-data-derived (per repo
convention).

## Rollout

Single PR (backend + frontend are small and coupled only through the healed data).

**SSE visibility (correction from review):** the re-attach loop runs inside
`_persist_analysis_inner!`, which server routes hit on every peak-curation mutation
and reanalyze — not just the CLI. The behavior change (no more stale-marking, ≥1
threshold, intents) is therefore visible in `analyze_run` `post_state` frames and
flows through the existing `applyRemoteToCache` `analyze_run` branch. No new
frontend branch is needed (post_state already carries all index kinds with status +
peaks); the contract test above pins this. No new event kinds, no `OpKind`, no API
shape change.

Prod redeploy + heal per runbook in §D.
