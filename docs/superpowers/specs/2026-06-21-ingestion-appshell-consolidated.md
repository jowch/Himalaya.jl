# Ingestion + App-Shell — Consolidated implementation spec

**Date:** 2026-06-21 · branch `ingestion-redesign` · Status: **active working reference**
**Supersedes (for implementation purposes):** `2026-06-15-ingestion-redesign-design.md` and
`2026-06-20-app-shell-unification-design.md`. Those remain the decision-history of record; this
doc is the single *current* target. These specs are disposable scaffolding — kept only until the
cutover lands.

**Canonical visual reference:** `docs/redesign-mockups/app-shell/` (rescued from `/tmp`, committed
`19407ea9`). The Jun-18 `docs/redesign-mockups/*.html` set is **stale** — do not use it.
Key renders: `p1-gallery.png` · `p2-corpus.png` · `bars.html`/`b1–b4.png` (dock per surface) ·
`p1-new.png` · `p1-scanning.png` · `p1-grouping.png` · `p1-config.png` · `p1-failed.png` · `flow.html`.

---

## 1. Goal

Bring the live UI/UX **and supporting backend** up to the app-shell-unification mockups,
**pixel-faithful** (improving on the mockups only where it is a genuine upgrade). The hard backend
(grouping core, events, scheduler, migration, structural mutators, directory-picker endpoints) and
the unified TopNav already exist. The remaining work is overwhelmingly **presentational + wiring +
a few backend stat additions**.

## 2. Where the branch actually stands (live sweep, 2026-06-21)

Served a copy of the real regrouped DB (3 experiments · 46 loads · 505 samples · 1985 exposures)
from source and walked every surface against `docs/redesign-mockups/app-shell/`.

**Landed & correct:** unified `TopNav` (wordmark · Experiments | Series); single shell; the `Dock`
*component* (grammar-correct on Focus / Loupe / Series); grouping core, geometry, events, scheduler,
migration, structural mutators, `routes_fs` directory-picker, `get_loads_rollup` + flags.

**An earlier iteration, NOT at the mockups** (the gap this spec closes): the experiment surfaces
were minimally wired, never refined.

## 3. Settled IA decisions (this session)

1. **No Corpus | Configuration tab bar.** Configuration moves behind the **⚙ gear** on the experiment
   header; **Corpus is the experiment home**. Grouping review stays a banner-reached surface.
2. **Grouping review = live "review as loads land"** with a **Confirm-groups gate** (mockup
   `p1-grouping` / `p1-scanning`), not the static post-scan list. Loads appear and unfold during the
   scan; settled loads are immediately reviewable; Confirm unlocks when the scan finishes.
3. **Sequencing: backend stats first** — land the list-endpoint stats + sessions/span before the
   frontend surfaces (unblocks the gallery cards + stat ledger).
4. Everything else: match the mockups.

## 4. Target state per surface (the deltas)

Severity: **P0** core daily-loop · **P1** spec'd surface wrong/missing · **P2** degraded · **P3** nit.

### 4.1 Experiments gallery — `/experiments` · ref `p1-gallery.png` — **P0, large**
- **Now:** sparse list (status dot + name + dir).
- **Target:** year-grouped **card gallery** + a slim left **timeline rail** (year · N sessions);
  each card = serif name (`SSRL April 2026 · 1p7m`, middot not hyphen) + **state chip**
  (`up to date` / `N to review` / `scanning…`) + mono `data_dir` + counts line
  (`170 samples · 682 exp · 13 loads`); a `scanning…` card shows `indexing N files…`; title
  `All experiments`; `+ New experiment` (terracotta). ⚙ on the top bar far-right.
- **Backend prerequisite (§5.1):** the list endpoint must carry per-card stats.

### 4.2 Experiment Corpus — `/experiments/:id/corpus` · ref `p2-corpus.png` — **P0, large**
- **Now:** Corpus/Config tabs; columns `FRAMES KEPT · TAGS · STATUS(Index→)`; "Needs review" on
  every row; **no Dock**; empty EXPOSURES column; terse banner; 4 stats; `SESSIONS 0`.
- **Target:** **no tabs**; experiment header (kicker · serif name · dir · `● scanned 2 min ago ·
  up to date` · **Rescan** · **⚙**); **5-stat ledger** (`loads · samples · exposures · span · sessions`);
  explanatory amber banner (`3 samples need a grouping check. Two stage-position jumps and one
  counter reset were ambiguous.` → `Review grouping →`); sheet columns `SAMPLE · EXPOSURES · KEPT ·
  PHASE`; rows = checkbox + serif name + id + `slot N` chip + **detector thumbnails** (current frame
  outlined, `+N`) + `KEPT n/m` + **phase chip** (`Pn3m`/`Im3m`/`not indexed`) + `Open ›`; current row
  terracotta-tinted; **the bottom Dock** (§4.3). Drop the per-row "Needs review" pill (noise).

### 4.3 The Dock (two-tier shell, bottom) — ref `bars.html` / `b1–b4.png` — **P0, medium**
- **Now:** renders on Focus/Loupe/Series; missing count readouts, `Sample`/`Frame` labels, key-chips.
- **Target:** persistent on **every** surface incl. **Corpus** (currently absent). Grammar
  `‹ up-link | cursor (steppers w/ N/total readout) | cull (where frames exist) … destinations`,
  with key-chips (`L`, `X`, `K`, `↵`). Per-surface table (spec §3.3 of the app-shell design):
  | Surface | up-link | cursor | cull | destinations |
  |---|---|---|---|---|
  | Corpus | `‹ Experiments` | Sample ↑↓ · Frame ‹› | Drop·Keep·Restore | Loupe · **Focus** |
  | Focus | `‹ Corpus` | Sample ↑↓ | — | Loupe |
  | Loupe | `‹ Corpus` | Sample ↑↓ · Frame ‹› | Drop·Keep·Set-rep·Restore | **Focus** |
  | Series | `‹ All series` | Sample ↑↓ | — | **Focus** |

### 4.4 Grouping review — `/experiments/:id/grouping` · ref `p1-grouping.png` — **P1, medium**
- **Now:** static post-scan review (good bones: 3-level fold, split divider, Rename/Split, filter).
- **Target (decision §3.2):** live-unfold during scan — progress line (`↻ Parsing exposures… 418 /
  ~682 · 8 loads · 102 samples ▓▓▓ · 3 flags to review`); compact `load NN ▸ slot N` rows with inline
  flag reason (`stage jump → maybe same sample as slot 4`) + inline `Keep separate` / `Merge ↑` /
  `Keep` / `Split…`; `✓ clean` loads collapse to one line with `+ expand`; still-loading loads show
  `unfolding…`; a footer **Confirm-groups** bar (disabled while scanning). Reuses the existing
  structural-edit events; "Confirm groups" is the terminal commit. (Keep the detailed exposure-leaf
  view available on expand.)

### 4.5 New-experiment funnel — `/experiments/new` · ref `p1-new.png` — **P1, small–medium**
- **Now:** bare input + inline Cancel/Review.
- **Target:** carded **Directory** picker with live autocomplete dropdown + validation line
  (`✓ directory exists  ✓ not already an experiment`); a footer action bar (`Nothing is created —
  the next step indexes & lets you review.` · Cancel · `Review →`). Then Configuration approval gate
  → create + scan (`p1-config` → `p1-scanning`). Scan-failed recovery per `p1-failed.png`.

### 4.6 Configuration — behind ⚙ · ref `p1-config.png` — **P2**
- Components exist (`GeometryLedger` / `AcquisitionTimeline` / `SourcesCard`). Move off the tab onto
  the ⚙ route; fidelity pass against `p1-config.png`. `AcquisitionTimeline` needs real sessions (§5.2).

## 5. Supporting-backend changes

### 5.1 Widen `GET /api/experiments` (list) with per-card stats — **gates 4.1** — P0
Add to each list row: `load_count`, `sample_count`, `exposure_count`, `session_count`,
`span_hours` (max−min exposure timestamp), and `review_count` (flagged samples = the gallery's
`N to review`). Today only the *detail* endpoint computes stats; the list endpoint returns bare rows.
Reuse the `_experiment_stats` helper (`routes_experiments.jl:74`); fold it into the list builder.

### 5.2 Compute & expose **sessions** and **span** — gates the 5-stat ledger + gallery — P1
`SESSIONS` reads `0` everywhere and `span` is absent. Macro-session segmentation (multi-hour gaps,
spec §5 step 2) is conceptually in the grouping core but not persisted/surfaced. Decide: persist a
`session_id`/session table at ingest, or derive span+session-count at read time from exposure
timestamps. Read-time derivation is the lazier path if no surface needs session *entities* yet.

### 5.3 Everything else is already serving correctly
Loads rollup + flags, structural edits, geometry override, scheduler, `broadcast_progress!`,
directory-picker endpoints. No backend change needed for the Dock / grouping-live-unfold beyond
the existing `ingest_*` SSE frames (already wired).

## 6. Milestone plan (backend-first)

- **M0 — harness.** Mini-experiment: a local 2–3-load slice of real data
  (`~/projects/himalaya-devdata/mini-1p7m/`, exp1 loads 1–3 = 142 exposures incl. the agbe
  calibration) + a fresh dev DB, so the live scan/unfold/Confirm-groups UX is fast to iterate. *(in
  progress)*
- **M1 — backend stats** (§5.1, §5.2). List-endpoint stats + sessions/span. TDD against the route
  buckets. Unblocks M2.
- **M2 — gallery** (§4.1). Year-grouped card gallery + timeline rail + state chips + counts.
- **M3 — corpus + dock** (§4.2, §4.3). Remove tabs (⚙ for config); reskin sheet (KEPT/PHASE +
  thumbnails); wire the Dock onto Corpus; explanatory banner; 5-stat ledger.
- **M4 — funnel** (§4.5). Carded picker + validation + footer bar; Configuration approval →
  scan → scan-failed recovery.
- **M5 — grouping live-unfold** (§4.4). Progress line + load unfold + compact slot rows + Confirm
  gate. (Uses M0's mini-experiment to debug.)
- **M6 — configuration fidelity** (§4.6) + dock polish (key-chips, count readouts) sweep.

Each milestone: TDD where logic exists, pixel-faithful to the named mockup, **verified live in the
browser before moving on** (per the standing render-verify mandate). Branch stays unmerged.

## 7. Open questions
- **5.2 sessions:** persist a session entity vs derive at read time? (Lean derive unless a surface
  needs session entities.)
- **Gallery card meta line:** the mockup shows `SSRL · April 2026 · lipid mesophase series` — richer
  than the data model carries (`name` + null `description`). Source the meta from name parsing, or
  populate `description` during ingest, or drop the third clause?
- **Legacy `/samples`:** the flat global corpus (`SamplesPage`, still dock-wired) is slated for
  retirement (app-shell §2). Confirm removal vs keep as a hidden power-user route.
