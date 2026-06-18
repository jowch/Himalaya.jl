# Ingestion redesign — design spec

**Date:** 2026-06-15 (revised 2026-06-18 after a backend/schema/frontend research sweep)
**Status:** Design approved in brainstorm; pending review before an implementation plan.
**Register / design system:** product · "The Print" (see `DESIGN.md`, `PRODUCT.md`).
**Visual reference:** reference-grade mockups in `docs/redesign-mockups/*.html`. These are the artifacts implementation should mirror.

## 1. Goal

Make ingestion frictionless. Today a user hand-authors a `manifest.csv` (sample → exposure-stem grouping) and an `experiment.toml` (geometry + column mappings), then runs CLI `init`/`reingest`. The redesign replaces that with: **point Himalaya at a beamtime directory; it scans, groups the exposures into samples, derives the geometry, and lands you in the corpus.** North star: "dump a directory and it figures it out," with the smoothness of Gradescope's exam-ingestion flow.

The whole motivation is scale and friction: a beamtime is hundreds of exposures across ~a dozen rack-loads, and writing a manifest by hand for that is the cost we are removing.

## 2. Principles

1. **Assume as little as possible.** Infer jointly from timing, filename, and stage position; hard-code the minimum. Frame-count per acquisition is variable; sweep order is not guaranteed; conventions differ per experimenter. Numeric thresholds are *validated defaults*, not constants (§5).
2. **The database is the source of truth.** No hand-maintained manifest file. Ingestion appends rows; the working state lives in the DB.
3. **Additive, no "re-ingest."** Ingestion is one repeatable flow. Re-running it (or rescanning) only adds genuinely-new exposures and never clobbers human edits or curation. There is no reconcile-the-whole-manifest step.
4. **Surface the decision, hide the settled.** Auto-grouping clears the great majority; the human's attention is pulled only to the few ambiguous groupings.
5. **Lean to purpose, per surface.** Each screen does one job (e.g. the Samples screen verifies grouping; it does not show phase).
6. **Reuse the rails we have.** The event log, idempotency, SSE, the per-exposure analysis pipeline, and the background-timer pattern already exist. The redesign extends them; it does not build a parallel stack.

## 3. The flow

```
Point at a directory  →  Auto-process  →  lands in the experiment (corpus)
                          (scan · read PRP/setup · group · analyze)
        ↻ Rescan picks up new files later (cheap change-check; adaptive schedule; never a blind full re-read)
```

- **Point at a directory.** The user creates an experiment by selecting its data directory. An experiment **is** one setup = one root directory.
- **Auto-process** (staged, with live progress over SSE): enumerate exposures (TIFF + PRP + `.dat`); read geometry from PRP + beam center from the analysis-dir setup files (§6); group exposures into samples (§5); run per-exposure analysis (peak-finding + indexing) in the background.
- **Lands in the experiment.** No separate "ingest review" gate; the user lands in the experiment's Samples/corpus view, which populates live as processing proceeds.
- **Rescan.** A rescan never blindly re-reads everything. It first runs a **cheap change-check** (directory listing + max mtime vs. the stored scan signature); if nothing changed it is a no-op. If new files appeared, it ingests *only* those (insert-only, never clobbering human edits) and analyzes them. Rescan runs (a) on an **adaptive auto-schedule** while a beamtime is live and (b) on a manual **Rescan** action. There is no manual "add files." The auto-schedule (§9.4) starts polling when the first samples arrive and **backs off in tiers** (hourly → daily → off) as the run goes quiet, re-arming if new files appear, bounded by the fact that a beamtime lasts at most ~a week.

## 4. Information model

Two classes of information, opposite sources: **grouping + identity** is human knowledge; **everything instrumental** is derivable. The redesign splits them.

| Entity | Fact | Source |
|---|---|---|
| **Experiment** (= one setup = one directory) | name, description | human / auto (dir basename) |
| | geometry: energy, flight-path, detector, pixel pitch, q-units | **PRP** (constant across the experiment) |
| | beam center, calibrated detector distance | **analysis-dir setup files** |
| | data dir, analysis dir, file patterns | config (defaults; editable) |
| | scan signature (last scan time, file-set fingerprint) | computed |
| **Load** (one rack-load / "set"; §5 terminology) | index, session, time range, frame count | computed (time-gap segmentation) |
| **Sample** (a rack slot's specimen) | name, notes | human / auto-generated (§5 naming) |
| | grouping source (auto-position / auto-time / manual / user-merged), name source (auto / human) | computed / human |
| | which exposures (the grouping) | auto (§5) + human correction |
| **Exposure** (one TIFF + PRP + `.dat`) | image / trace / prp paths | computed from disk |
| | accept/reject status, selection | human curation |
| | content fingerprint (dedup) | computed |
| | **timestamp, exposure time** | PRP (the only per-exposure PRP fields kept) |
| | grouping signals: stage H-position, scan id, frame no | PRP / filename |

**No sample "role" field.** References (calibration `agbe`, background `emptycap`/water) are just samples the user adds if they want them; Himalaya does not model them specially. **Most PRP fields are dropped** (motor slits, attenuator, `dispx`, `Phi` — constant or noise to the user).

**Sample identity — id + a single `name` label.** A sample has one numeric identity, `samples.id` (autoincrement PK), used by every FK (`exposures.sample_id`), URL (`/sample/:id`), the mutation queue, and SSE. It has one human label, **`name`** (`HA85 (S01P15)`), auto-populated and user-editable. Today the schema carries *two* text columns — a machine identifier (`name` = the manifest id, e.g. `"JC001"`) and a friendly label (`display_name`); the redesign **drops the machine identifier** and **renames the surviving label `display_name` → `name`**, so a sample ends with exactly one label field. The dropped machine identifier's live uses are all subsumed or dead: display fallback (subsumed — the label is never null now), NavModal search (works on the label), export (emit the label), reingest upsert-by-name (`cli.jl:192`, dead with the manifest), CLI `-s <name>` filter, and the `[A-Za-z0-9._-]` validation (`validate.jl:3`). It was **not** a live permalink slug — URLs are numeric; `/api/resolve?sample=<name>` only back-translated *retired* pre-greenfield `/index/<exp>/<name>` bookmarks, which now degrade to `/samples`. The manifest CSV's own `sample_id` column was never the DB id and also disappears. Across rescans a sample is re-identified by its **(load, slot) coordinate** (§5), not by any text key.

**Provenance.** Every harvested fact is stamped with its source: `human | prp | setup | computed` (the `.dat` is dropped as a provenance source — it is only the Q/I/error integration and carries no identity or geometry). Provenance drives: (a) **dedup** — an exposure's identity is its `(experiment, filename)` + content fingerprint, so re-scanning only adds new files; (b) **never-clobber** — human edits (sample name, status, curation, geometry override) survive rescans while derived facts refresh (this generalizes the `display_name` reingest bug at `cli.jl:203`, see `docs/experiment-config.md`); (c) future LLM-sourced facts would be `proposed` pending confirmation.

**Name provenance is internal, not a Configuration-page badge.** Only *geometry* fields show a visible source tag on the Configuration page (PRP / setup files / you). A sample's `name_source` (auto vs. human) is tracked solely to power never-clobber (a human-edited `name` survives rescans); it is not surfaced as a source chip.

## 5. Grouping algorithm (data-derived)

Derived empirically from the real SSRL dataset (`/Volumes/data/ssrl/2026_04/1p7m`, 682 exposures; full analysis in the brainstorm). The workflow is a **stepping sample rack**: the stage sweeps through fixed slot positions shooting each, then the rack is reloaded and swept again.

**Terminology.** A **load** is one loading of the sample holder (synonyms the user uses: "set", "rack"); the holder has up to ~20 slots here but the count is not fixed in general. We use **load** as the canonical structural term throughout (it matches the mockups' "Load 1, Load 2…"); "set"/"rack" are accepted synonyms. A **slot** is one position within a load; a **sample** is the specimen in one slot; an **acquisition** is one contiguous burst of **frames** (exposures) at a slot.

**Deterministic backbone (no LLM):**
1. **Frames → acquisition.** Consecutive frames at one stage position within a short burst window. A slot boundary is a **position gap large relative to the local within-burst jitter** (gap ≫ the median consecutive-frame position delta), *not* an absolute spacing — so the same logic handles both the HA ~4 mm sweep and the JC well plate (H ≈ 39.5 / 71 / 87 / 107 mm, 16–31 mm spacing) that **coexist in the same directory**. Frame count per acquisition is **variable** (commonly 4 in this data — 161 of 170 slots — but never hard-coded).
2. **Time-gaps → loads.** A gap beyond the load-threshold starts a new load/sweep (clean separation in the data: intra-sweep gaps ≤ 53 s vs. reload gaps ≥ 533 s). Multi-hour gaps mark macro-sessions. The directory yielded ~13 loads, ≤ 24 slots/load, 4 macro-sessions.
3. **Slot identity = stepping-axis position ⨉ filename token, cross-checked.** The stepping axis is the **motor with the dominant inter-burst variance** (Horizontal position in this dataset, where `dispx`/`Phi` are constant — but the axis is beamline-dependent, so it is *detected*, not hard-coded). Position and the filename token agree across this run; the token corroborates and may provide the label. Agreement → near-certain; disagreement or an unparseable convention → flag for the human. **Fallback when no axis varies** (single fixed position, single-exposure samples, or a setup that doesn't record the stepping motor): cluster by time + filename only and flag for review — position is the *preferred* backbone, not a required one.
4. **Each slot within a load = one distinct sample** (no series treatment; the user can group into a series later).

**The rack can stay in place and a subset of positions get reshot in any order** (no reload gap) — a position revisited within a load is another acquisition of that same sample.

**Thresholds are gap-relative and per-load, not global constants.** Slot clustering keys off the *ratio* of a position gap to the local within-burst jitter (slot tolerance = k × median consecutive-frame delta), and the load boundary keys off the bimodal time-gap histogram (intra-sweep ≤ 53 s vs. reload ≥ 533 s) — both derived **per load** from that run's own distribution, since slot spacing is itself bimodal across the directory (4 mm HA sweep vs. 16–31 mm JC plate) and no single global tolerance fits both. The measured figures (within-burst jitter ≈ 0.30 mm median, ±0.5 mm for the tight JC wells, ~4 mm HA spacing, ~300 s reload gap) are sanity bounds and starting points, not magic constants. This honors Principle #1: the numbers are re-validated per run, never baked in.

**Naming.** One human label per sample, `name` (the single label field after the rename in §4):
- **Default:** the filename-derived label with the set/position coordinate as a parenthetical anchor — `<label> (SNNPMM)` (season/episode style: `S` = set/load index, `P` = slot position, zero-padded — e.g. `HA85 (S01P15)`), or a confident per-slot token like `JC C04 (S02P05)`.
- **No character restrictions** — the label is plain TEXT, only whitespace-trimmed (`routes_samples.jl:158`). The old `[A-Za-z0-9._-]` rule was a *manifest-parse* validation (`validate.jl:3`) the redesign retires; samples are addressed by numeric id, so there is no URL/slug constraint.
- The `SNNPMM` coordinate is the deterministic, always-available backbone; the filename label is corroboration. (HA is exactly the case that needs it: the HA filename number is a *global frame counter*, not a slot identity, so the coordinate disambiguates.)
- The user can rename anytime; a human edit sets `name_source = human` and is never clobbered on rescan. Across rescans a sample is re-identified by its **(load, slot) coordinate**, not by its label — so renames are safe.

**Needs the human / label (surfaced, never blind):**
- **Cross-load identity.** Slot number is **not** stable across loads (the same sample can land in a different slot; the same slot holds a different sample after a reload). So cross-load "is this a reshoot of that sample?" comes from the filename label, not position, and is surfaced as a Merge prompt (default: each load's slots are new samples). A counter resetting to 1 signals a fresh rack.
- **Within-load split.** When stage position jumps mid-group, suggest a Split at that boundary.

**Regroup is bookkeeping only — no reanalysis (within one experiment).** Verified against the schema: peaks (`auto_peaks`), curations (`peak_curations`), indices (`indices`/`index_peaks`), and assignments (`assignment_members`) are all keyed to `exposure_id`; `samples` is a pure container with no computed column. Changing which sample an exposure belongs to touches no derived data — **so long as the move stays within one experiment**: `analyze_exposure!` resolves `analysis_dir` + the integration pattern through the exposure's experiment, so a cross-experiment move would change the resolved `.dat` path and require reanalysis. Structural edits are therefore constrained to within-experiment (§9.2). (This retires the earlier brainstorm decision that "editing grouping reanalyzes affected samples" — it predated confirming analysis is exposure-keyed and geometry-independent. A saved Series is preserved by re-pointing membership on merge, never by deleting a sample — §9.3.)

## 6. Geometry derivation

No PRP-parsing code exists today (geometry is hand-entered in `experiment.toml`'s `[beamline]` block and read defensively by `_beamline_from_config` at `routes_experiments.jl:16`). The redesign adds it.

- **From PRP** (constant across the experiment): `Beam energy` → energy_kev (9000 eV → 9.00 keV), `Pipe length` → flight_path_m (1700 mm → 1.70 m), `Detector` → model → pixel pitch via a small **detector→pitch lookup** (a static table in `geometry.jl`, e.g. Pilatus 1M → 172 µm; an unknown detector leaves pitch unset and **flags for manual entry** rather than guessing), q-units.
- **From the analysis-dir setup files:** beam center (`Beam center is at (421.409, 836.946)`) and, if present, the calibrated detector distance (e.g. 1809.5 mm, which is distinct from the PRP `Pipe length`; Himalaya's geometry uses `flight_path_m`, so the calibrated distance is recorded for reference/override but the flight path stays authoritative unless overridden). Multiple setup files → use the latest by mtime.
- **Per-field source tag** (`prp | setup | human`, default `default` before derivation). An **Override** affordance lets the user correct any value; overrides are `human`-sourced and never refreshed on rescan.
- **Multi-setup discrepancy detection.** A directory should hold exactly one setup, but since we read every PRP during the scan anyway, we **verify the constant fields actually are constant** (beam energy, pipe length, detector) and surface any variance as a likely error to correct (a "geometry check found N issues" banner on the Configuration tab), rather than silently averaging or picking the first. This catches data-entry/beamline mistakes for free.
- This removes the hand-edited `[beamline]` block from `experiment.toml`.

## 7. Navigation & information architecture

- **"Experiments" is the nav-default** (and the surface shown when no working experiment is cached). It is presented as a **left sidebar list** of experiments (name + counts), with **+ New**. (Today the frontend is single-corpus with an `activeExperimentId` in Zustand and a `?beamtime` filter chip; this elevates the experiment to a first-class navigable entity.)
- Selecting an experiment **persists the choice to Zustand**; return visits route straight to that experiment's corpus.
- **An experiment is one page with two tabs** at `/experiments/:id` (Configuration | Samples), **not two pages**. A **shared experiment header** sits above the tabs on both: kicker · serif experiment name · rescan status (or ingest progress) · a 2px rule · a stats ledger (`loads · samples · exposures · span · sessions`) · the tab bar. Only the content below the tabs swaps.
- **Soft-lens scope** (settled earlier): the experiment scopes the corpus and is the entry, but the data model is not partitioned — Series can span experiments and samples stay globally addressable.

## 8. UI surfaces (the reference)

All in The Print: warm paper/ink, rationed terracotta, Newsreader for titles, mono for measurements, detector frames as dark windows, flat-except-the-plate, 5px radius. Faithful to `src/print/ui/` primitives (Card, Button, SegmentedControl/StageTabs, Field, EmptyState, ProgressBar, Thumbnail; PhaseStrip is **not** used here). Five states are mocked:

### 8.1 Configuration tab — `experiment-configuration.html`
Below the shared header: editable description + directory; a two-column body of **Geometry** (the auto-derived ledger, source-tagged, with Override, and the multi-setup discrepancy banner when triggered) and **Acquisition** (a timeline plate of the sessions/loads, the one lifted figure); then **Sources** (data/analysis directories + the `{name}.dat/.tiff/.prp` patterns). This is the experiment's identity + setup verification.

### 8.2 Samples tab — `experiment-samples-review.html` (the centerpiece)
Purpose, narrowly: **verify every sample loaded, has its exposures, and is split/grouped correctly.** No phase, no scores. A **three-level fold**:
- **Load** (collapsible): operator · N samples · M exposures · time · status (`✓ grouped cleanly` / `⚠ N to check`).
- **Sample** (collapsible): serif name (e.g. `HA85 (S01P15)`) · "N exposures · time" · status; controls (Rename / Split / for flagged: Merge…/Split… with a dismiss option) in the fold header.
- **Exposure** (leaf, all shown — **no `+N` truncation**): small detector thumbnail · **filename** (mono) · stage H-position · time. Showing every filename + position is what makes the split *verifiable*; a `position jumps X → Y mm` divider marks an auto-suggested split boundary.

Flagged samples are pre-expanded; clean ones fold to one line. A per-exposure **Move…** allows fine regrouping.

### 8.3 First ingest (live unfold) — `experiment-first-ingest.html`
Woven into the Samples tab, not a separate screen. The header shows ingest progress (`Processing · 312 / ~680 exposures`) instead of "up to date"; the stats ledger shows discovered-so-far counts with span/sessions pending; loads appear and fill as they are scanned (done · in-progress with a per-load bar · queued), with a footer "scanning · 5 of ~13 loads." Driven by `ingest_*` SSE frames (§9).

### 8.4 Empty / first-run — `experiments-empty.html`
An EmptyState: "No experiments yet" + the value proposition + a primary **New experiment** action (point at a directory). Sidebar reads "none yet."

### 8.5 Failed scan — `experiment-scan-failed.html`
Honest error (not a masquerading empty): "No exposures found in that directory" + a mono detail box (path · `looked for {name}.tiff + {name}.prp · matched 0 of 0 files`) + recovery (Retry scan / Edit file patterns / Choose another directory). The sidebar flags the experiment.

## 9. Backend & API surface

The redesign is mostly net-new backend work, but it slots onto existing rails.

### 9.1 Shared grouping core (new modules)
A single ingest core, callable from both the HTTP API and the CLI:
- **`prp.jl`** — `parse_prp(path) → NamedTuple` (key=value lines; strip units; defensive: unparseable → `missing`).
- **`geometry.jl`** — `derive_geometry(prp_paths, setup_files) → (geometry, discrepancies)`; samples PRPs, verifies the constant fields, maps detector→pitch, parses the latest setup file for beam center.
- **`grouping.jl`** — the §5 backbone: `scan_directory(dir, patterns)` (filesystem discovery, reusing `resolve_files`/`resolve_file_path` from `config.jl:202`), then `group_into_samples(meta) → loads/samples`. (Named distinctly: `auto_group` already exists in `pipeline.jl:30` for *index*-candidate grouping — an unrelated concept.)
- **`scan_and_group!(db, exp_id, dir; additive)`** — the transactional orchestrator. Reuses `_DB_WRITE_LOCK` (`server.jl:29`), the per-exposure `analyze_exposure!` (`pipeline.jl:913`, idempotent + hash-guarded), and the insert-only upsert discipline already in `_reingest_inner!` (`cli.jl:156`). On rescan it inserts only files absent from the dedup key and analyzes only those.
- **Required edit to `analyze_exposure!`:** today it resolves the experiment (and thus `analysis_dir` + integration pattern) by JOINing `exposures → samples → experiments` (`pipeline.jl:919`), so an exposure with no `sample_id` yet falls out of that JOIN and errors. Rewire it to read the new denormalized `exposures.experiment_id` directly. Ingest order: group into samples first, then analyze (or decouple analysis from `sample_id`) — analysis must not depend on the sample assignment.

### 9.2 REST endpoints (new / changed)
- **`POST /api/experiments`** — create-from-directory: `{ path, name?, patterns? }`. Returns the experiment id immediately and kicks off the first scan asynchronously; progress streams over SSE. (Replaces the pre-create-`experiment.toml`-then-`init` dance.)
- **`POST /api/experiments/{id}/scan`** — rescan: cheap change-check, then additive ingest of new files. Idempotent (a no-change scan does nothing). Supersedes `POST /api/experiments/{id}/reingest`.
- **`GET /api/experiments/{id}`** — extended to include typed geometry + per-field source, ingest status, and the loads/stats roll-up.
- **`GET /api/experiments/{id}/loads`** (or fold into the samples response) — the Load ▸ Sample ▸ Exposures structure for the Samples tab.
- **Structural grouping edits** — `rename` / `move` are thin single-entity event-sourced routes; `merge` / `split` are route-level orchestrations over per-exposure moves + a sample create/retire (§9.3). All are **within-experiment only** — a cross-experiment move would change the resolved `analysis_dir`/integration pattern and invalidate the analysis (§5). None trigger reanalysis for same-experiment moves.
- **Geometry override** — `PATCH /api/experiments/{id}` writes the field + sets its source to `human`; never refreshed on rescan. This **replaces** today's read-only PATCH stub (`routes_experiments.jl:73` returns "experiment metadata is read-only" for every PATCH) and requires widening the frontend `updateExperiment` helper, whose patch type is currently `Record<string, never>` (`api.ts`).

### 9.3 Event kinds

**The event log is single-entity and replayed.** `apply_event!` keys each event to one `entity_id`, and `rebuild_views_from_log!` re-folds events per `(entity_type, entity_id)` from empty. An edit is "free via the dispatcher" only if it mutates exactly one entity and is replay-idempotent. That shapes the classes below.

- **`exposure_moved` and `sample_renamed` are genuine single-entity events** (durable, replay-safe). Moving one exposure to a different sample is keyed to that exposure's `entity_id`; renaming is keyed to the sample. These inherit idempotency + SSE + undo (`undoes_event_id`) and can splice the cache directly.
- **`sample_merged` / `sample_split` are orchestration, not single-entity view-writes.** They touch *two* sample rows plus N exposures, and minting sample rows via `AUTOINCREMENT` is **not** replay-idempotent through `update_view_for_event!` (a replayed "create" mints a new id, breaking FKs). So: (a) implement them in the route as a sequence of per-exposure `exposure_moved` writes plus an explicit sample create/retire, each authored directly — durable `user_actions` rows for audit/SSE, with the route writing the **inverse** for undo (a single `undoes_event_id` stamp cannot reverse a multi-row reassignment); (b) **merge is a re-point, never a delete** — reassign the loser's exposures to the survivor and re-point its `series_samples` / `sample_tags` / `sample_messages` (all key to `sample_id`, and `series_samples` is `ON DELETE CASCADE` at `db.jl:820`, so a naive delete silently drops the sample from saved Series); (c) the **loads/samples grouping view is re-derived** from the `(load_id, slot_index)` coordinate (a derive, per the plotting-redesign "rebuild-not-log-derivable" precedent), not folded from the log.
- **SSE cache handling for merge/split is invalidate-only.** Because they create/destroy sample rows whose ids aren't payload-derivable, their `applyRemoteToCache` handlers **refetch** the `experimentSamples`/`loads` query (the `series_created` precedent), not surgical splices.
- **Transient ingest progress** broadcasts **without** a durable `user_actions` row — but **not** via the `analyze_run` path (that precedent is the inverse: it *writes* the row and *suppresses* the broadcast). `broadcast_event!` currently requires an `event_id` from a `user_actions` INSERT, so this needs a thin **new** helper (`broadcast_progress!`) emitting a frame with a sentinel id. Kinds: `ingest_started {experiment_id, total_estimate}`, `ingest_progress {experiment_id, processed, total, new_load_ids, new_sample_ids}`, `ingest_complete {experiment_id, counts}`, `ingest_failed {experiment_id, error}`. The frontend must route `ingest_*` **before** the `entity_id`-required intake guard (their key is `experiment_id`, not an exposure `entity_id`) into explicit `applyRemoteToCache` cases — never the default branch, which treats `entity_id` as an exposure id and would fire bogus `peaks`/`indices` invalidations — updating the ephemeral `ingestInFlight` store and invalidating the experiment-samples query.

### 9.4 Auto-rescan scheduler
Follows the existing GC-timer pattern (`start_gc_timer!`, `server.jl:144` — a `Timer(0; interval)` held in a module `Ref`):
- A **per-experiment poll Timer** keyed by `experiment_id`, started when the first samples arrive from a scan.
- Each tick runs the **cheap change-check** first; only a changed directory triggers ingest + analysis + broadcast. Unlike the GC timer's fast single-DELETE callback, a scan can outlast its interval, so the tick needs a **reentrancy guard** (skip-if-already-running, or stop/restart the Timer around each scan) and the heavy scan/analyze work must run **off the timer callback** (`@spawn`) so it doesn't stall the libuv timer loop.
- **Tiered backoff** (exponential-retry style): start at the **fast tier** (every 1–6 h) while exposures are actively arriving; after N consecutive empty ticks (the run has gone quiet) back off to **daily**; after a few more empty daily ticks, **stop** the Timer → the experiment drops to **manual-rescan-only**. Bounded by max beamtime ~1 week. A manual Rescan (or auto-tick) that finds new files **re-arms** the schedule at the fast tier.
- Timers are torn down on server shutdown (mirror `stop_gc_timer!`) and on experiment delete — **note no experiment-delete route exists today**, so the scheduler work includes building that delete path with timer teardown. All writes go through `_DB_WRITE_LOCK`.

### 9.5 Fate of the CLI
The CLI stays as a headless/scriptable path and for server boot (`serve` is unchanged). `init`/`reingest` are rewritten to call the **same `scan_and_group!` core** instead of the manifest path; manifest support is retained only as a legacy fallback (or dropped during migration — decided in the plan). The grouping engine is shared, so CLI and API never diverge.

## 10. Schema direction

Additive; the existing schema already carries some bones (source-tagged `sample_tags`/`exposure_tags`, `trace_hash`/`analysis_inputs_hash`, `exposure_sources`). Concrete sketch (exact DDL + a `schema_migrations` sentinel belong to the implementation plan):

- **experiments:** promote geometry to typed columns (`beam_center_x`, `beam_center_y`, `pixel_size_um`, `q_units`, `detector_distance_mm` — today these live only in the TOML blob) alongside the existing `energy_kev`/`flight_path_m`; add a **per-field source** column for each (`*_source ∈ {prp, setup, human, default}`); add `last_scanned_at` + `scan_signature` (max mtime + file count) for the cheap change-check; retire `manifest_path`.
- **loads:** **new table** (`id`, `experiment_id`, `load_index`, `session_id`, `start_time`, `end_time`, `frame_count`, `note`) so Load ▸ Sample ▸ Exposures is first-class and queryable; computed at ingest from the time-gap segmentation.
- **samples:** still a pure container; collapse the two text columns to one. Migration order (to avoid a transient name collision): **drop the manifest machine `name`** and its `UNIQUE(experiment_id, name)` index (`db.jl:1492`), then **`RENAME COLUMN display_name → name`**. Final shape: `(id, experiment_id, name, notes)` + `load_id` (FK), `slot_index`, `grouping_source ∈ {auto_position, auto_time, manual, user_merged}`, `name_source ∈ {auto, human}` (never-clobber: a `human` name survives rescans). Across rescans a sample is re-identified by its **(load_id, slot_index)** coordinate, which replaces the old by-name upsert match. The manifest-era machinery (reingest upsert-by-name `cli.jl:192`, CLI `-s <name>` filter, charset validation `validate.jl:3`) and the by-name `/api/resolve` path are deleted.
- **exposures:** add `experiment_id` (denormalized — the dedup key moves to `(experiment_id, filename)` since `sample_id` is ephemeral during grouping, and `analyze_exposure!` reads it directly to resolve the experiment — §9.1), `prp_path`, `timestamp`, `exposure_time`, `horizontal_position` (or the detected stepping-axis position), `scan_id`, `frame_no`, `load_id`, and a content fingerprint for dedup.
- **dedup key:** `UNIQUE(experiment_id, filename)` (replacing today's `UNIQUE(sample_id, filename)` at `db.jl:1551`); rescan inserts only absent files. **Migration order:** add `exposures.experiment_id` → backfill from the `samples` JOIN → drop the old `UNIQUE(sample_id, filename)` index → create the new index (SQLite cannot enforce the new uniqueness until the column is backfilled and the old index dropped; mirror the dedupe-then-enforce discipline of `migrate_exposures_unique_filename!`).
- **untouched:** `indices`, `auto_peaks`, `peak_curations`, `index_peaks`, `index_groups`, `assignments`, `assignment_members` — all keyed to `exposure_id`; re-grouping reassigns `exposures.sample_id` only and leaves every derived table intact (§5).

## 11. Implementation TODOs (specify before/while building)

- **Net-new core (build spine, phase first):** the modules and surfaces below do not exist yet and gate everything else — `prp.jl` / `geometry.jl` / `grouping.jl` (+ `scan_and_group!`), the `POST /api/experiments` create-from-directory endpoint, the `loads` table + the exposures/samples/experiments migrations (§10), and the `broadcast_progress!` helper (§9.3). Phase the plan core-first: schema + grouping core → scan/rescan API + SSE progress → structural-edit events → frontend `/experiments/:id` page.
- **Guarded structural edits:** `Merge…` / `Split…` / `Move…` open an inline confirm that previews the result and commits with **undo** (the critique's P0 — regrouping is the product, it must be safe). `undoes_event_id` backs single-row undo (rename, one move); merge/split undo is **route-authored inverse writes** (§9.3), since one stamp cannot reverse a multi-row reassignment. Demote primary actions from filled-accent to outline-accent so the dismiss option is co-equal.
- **Auto-split threshold:** the position jump that triggers a split suggestion — gap-relative per §5 (a jump ≫ the local within-burst jitter), derived per load, not a fixed ~4 mm.
- **`name`/`display_name` consolidation sweep:** collapsing to one label column touches the whole stack. (1) Delete the manifest machine `name` and its by-name uses: the by-name `/api/resolve` path + `IndexSlugRedirect` (or repoint `/index/*` straight to `/samples`), reingest upsert-by-name, CLI `-s <name>`, charset validation. (2) Rename `display_name` → `name` everywhere: API responses + `Sample`/`CorpusSample` TS types, the `sampleDisplayName` helper + every `display_name ?? name` fallback chain (now just `name`), NavModal search arrays, comparison/scoping labels, the PATCH route, tests, **and the queue layer** — the `update_sample` branch in `applyRemoteToCache.ts` (the SSE payload key flips to `name`) and `lib/queue/mutators/trivial.ts` (`UpdateSampleInput`, `patchOf`, the `onSuccess` that reads both fields). **Delete** (don't rename) NavModal's secondary `name`-vs-`display_name` line — it's an artifact of the two-label model. Sequence so no reader references a column mid-migration; land it as a focused refactor commit ahead of the feature work.
- **Accessibility / efficiency:** `aria-expanded` on the folds, a keyboard path through Load/Sample/Exposure, and a bulk-select for renaming/regrouping many samples at once (none shown in the mockups; required for 170-sample scale).
- **Rescan cost model:** cheap change-check (listing + mtime vs. `scan_signature`) before any hashing; hash/parse only new or changed files. The adaptive schedule (§9.4) bounds background work.
- **Frontend scaffolding:** new `/experiments/:id` layout route with Configuration/Samples child tabs (`AppRoutes.tsx`); ephemeral `ingestInFlight` store field — added to `AppState` + a named action but **NOT** to `partialize` and **NOT** a persist-version bump (it's ephemeral; bumping `version` without a key-preserving migrate risks wiping the persisted blob); `queryKeys.experimentConfig(id)` / `experimentSamples(id)` (confirm the latter is the Load▸Sample▸Exposure roll-up — a distinct shape/key from the existing `samples(experimentId)` key, not a collision); explicit `ingest_*` cases in `applyRemoteToCache` routed ahead of the `entity_id` guard (§9.3); widen `updateExperiment`'s patch type from `Record<string, never>` to a typed geometry partial; a live-unfold skeleton pattern (append new rows, don't reorder).

## 12. Out of scope / deferred

- **LLM-assisted extraction** (PromptingTools + a small local model reading the experimenter's notebook/spreadsheet to draft names/groupings). The deterministic backbone made this non-essential for v1; it remains a future enhancement that produces the same form, always human-reviewed.
- **Phase on the Samples screen** (removed by design — phase is not this screen's job and may not be computed yet).
- **Series treatment of slots** (each slot is one sample; series grouping is a separate, later user action).
- **Multi-setup auto-splitting** (we *detect and flag* geometry variance within a directory; we do not auto-split one directory into multiple experiments).

## 13. Open questions

- **Cross-load merge defaults per convention.** We default to "each load's slots are new samples; merge on matching label." Confirm this holds across experimenter conventions, or whether the counter-reset signal should auto-drive it.
- **Manifest legacy path.** Keep `manifest.csv` as a read-only legacy fallback for existing experiments, or migrate them and drop it entirely?
- **Geometry override granularity** (per-field vs. whole-block) and whether `flight_path_m` should ever defer to the setup file's calibrated detector distance.
- **Auto-schedule defaults.** Confirm the backoff tiers (1–6 h → daily → off): the fast-tier interval, the N empty ticks before dropping to daily, and the count of empty daily ticks before stopping. Expose as a per-experiment setting or keep implicit?
- **Corpus vs. per-experiment Samples IA.** The new `/experiments/:id` Samples tab overlaps the existing corpus contact sheet (`/samples` + the `?beamtime` filter chip). Does the redesign **retire** the flat-corpus `/samples` root and the chip in favor of experiment-scoped navigation, or do **both** coexist (a global corpus view *and* a per-experiment Samples tab)? §7's soft-lens implies coexistence, but that means two Samples surfaces with overlapping purpose — needs a deliberate call so the live-unfold lands on the right one.

## References

- Design system: `DESIGN.md`, `PRODUCT.md`.
- Today's ingestion contract: `docs/experiment-config.md`, `packages/HimalayaUI/src/{cli.jl,config.jl,manifest.jl,pipeline.jl}`.
- Backend rails the redesign reuses: `server.jl` (lifecycle, `_DB_WRITE_LOCK`, SSE, `start_gc_timer!`), `events.jl` (`apply_event!`, `update_view_for_event!`), `idempotency.jl` (`with_idempotency`), `pipeline.jl` (`analyze_exposure!`).
- Mockups + captures: `docs/redesign-mockups/*.html`.
- Memory: `project_ingestion_redesign` (decision log + data-derived parameters).
- Research sweep (2026-06-18): five-explorer workflow mapping API/CLI, pipeline, schema, PRP/geometry, and frontend — findings folded into §§4–11.
- Review sweep (2026-06-18): five-reviewer validation (himalaya/queue/saxs-physics/frontend reviewers + consistency) — corrections to the event architecture (§9.3), `analyze_exposure!` resolution (§9.1/§5), gap-relative grouping (§5), and the rename/scheduler/migration details folded in.
