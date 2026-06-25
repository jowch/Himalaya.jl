# Ingestion-redesign regression sweep

Live walk of every UI surface against a **copy of real production data**, to catch
regressions the mocked test suites cannot — the mocks always feed *ingested* data
(loads + grouping flags), but a real upgraded deployment has **0 loads** until it
rescans.

- **Date:** 2026-06-20
- **Data:** copy of `~/projects/himalaya-devdata/himalaya.db` (3 experiments, 139
  samples, 678 exposures; old schema, still carrying `display_name`).
- **Method:** served julia-from-source on a scratch port, drove a real browser
  (Playwright), checked every surface for render correctness + console errors.
- **Branch HEAD at sweep:** `01aedd30` → fixes landed at `019602fc`.

**Reproduce:** copy the dev DB to a scratch path, point `HIMALAYA_DB_PATH` at the copy,
`serve`, and open each URL below. (Finding #1 blocks the server from booting on the raw
copy — dedup the 30 colliding rows first, see below.)

Severity: **P0** blocks upgrade / data loss · **P1** broken surface · **P2** degraded ·
**P3** nit.

## Findings at a glance

| # | Sev | Surface | Status |
|---|-----|---------|--------|
| 1 | P0 | Migration cannot upgrade existing DBs | **Open — needs decision** (root cause resolved, fix recommended) |
| 2 | P1 | Experiment **Corpus tab is a blank unwired stub** | **Open — design call** |
| 3 | P2 | Experiments **home page has no nav chrome** | **Open — design call** (ties to single-shell rework) |
| 4 | P3 | Home experiment cards sparse vs mockup | Open — documented |
| 5 | P2 | Grouping empty-state shows misleading `No samples match ""` | **Fixed `019602fc`** |
| 6 | P3 | `--` placeholder; blank Acquisition card | **Fixed `019602fc`** |

**Legacy surfaces: ZERO regressions.** The `display_name → name` collapse (`7ee296c6`)
broke nothing. All of `/samples`, `/sample/:id`, `/series`, `/series/:id` (folio /
reading / builder / picker), and `/samples/loupe/:id` render real data with **0 console
errors**.

| Surface | URL | Result |
|---------|-----|--------|
| Contact sheet | `/samples` | 140 rows, real names + detector thumbs + phase chips ✅ |
| Focus | `/sample/2` | trace + peaks, Im3m a=220Å (0.41), candidates, detector rings, comb, filmstrip ✅ |
| Series folio | `/series` | "Bundle A", d3 overlay of 2 samples, filters/sort ✅ |
| Series reading | `/series/1` | 2-trace overlay (1.20× offset), Compose, Edit ✅ |
| Series builder | `/series/1` → Edit | drag-reorder real members, Add/Save/Cancel ✅ |
| Sample picker | builder → Add sample | grouped by beamtime, real samples + slugs ✅ |
| Loupe | `/samples/loupe/2` | detector frame + rings, Keep/Drop, rep-frame, filmstrip, kbd ✅ |

The new experiment surfaces render (StatBar wires correctly: `0 loads / 31 samples /
192 exposures / 0 sessions`), but the **experiment-scoped IA is load-centric**, so a
migrated-but-un-rescanned experiment shows real counts in the StatBar while every panel
below is empty — and nothing tells the user to rescan. That theme drives #2, #5, #6.
(The happy path — loads + flags — works; verified separately with a seeded DB.)

---

## #1 — P0 — Migration hard-fails on real data

`migrate_exposures_experiment_id!` (`packages/HimalayaUI/src/db.jl:573`) swaps the
exposures unique key from `(sample_id, filename)` → `(experiment_id, filename)` and
**`error()`s** when the same filename appears under two samples in one experiment. The
dev DB has **30 such groups** (exp 3 + 4). `open_db` re-throws on every boot, so the
server never starts — **no existing deployment can upgrade.**

This is a design-premise mismatch, not a typo: under the old model a file legitimately
belonged to more than one sample, and the migration anticipates the collision but
refuses to guess a merge.

### Resolved: are the duplicates the same underlying file?

**Yes.** Of the 30 collision groups, **all 30 share one `image_path`, 0 differ.**
Example: exposures `881` (sample 270, `D30-2`) and `905` (sample 244, `D30`) both point
to `/Volumes/data/ssrl/2026_01/0p7/data/JC117_S153_0_001.tif`. One row carries
`trace_hash` + `analysis_inputs_hash`, the other is empty (across all collisions: 32
analyzed rows, 28 not — about one analyzed + one un-analyzed pointer per pair).

So the duplicates are **redundant pointers to one physical file** (an old per-sample
re-grouping artifact, e.g. a re-shoot variant), **not genuinely-distinct exposures.**
Dedup-merge is therefore **safe, not lossy.**

### Recommended migration (supersedes the fail-fast)

For each `(experiment_id, filename)` collision where all rows share `image_path`: keep
the best row (prefer a non-empty `trace_hash` / `analysis_inputs_hash`; tiebreak lowest
id), drop the redundant rows, **then** create the unique index. Keep the hard `error()`
**only** for a collision whose rows have *differing* `image_path`s — a genuinely
ambiguous distinct-file case, which does not occur in this data.

Status: **awaiting go-ahead to implement.** Alternatives considered and rejected as
worse: relax the key (ingestion dedup relies on it); ship a separate manual cleanup tool
(strands operators); reconsider the key entirely (larger blast radius).

---

## #2 — P1 — Experiment Corpus tab is a blank unwired stub

`/experiments/:id/corpus` → `ExperimentCorpusPage` renders an empty `corpus-sheet-slot`.
The comment says "E2 mounts the scoped SheetTable here … wired in E2", but the E2 plan
scoped only the grouping-review surface (`/grouping`) — the corpus contact sheet was
never wired by any phase. The **primary view of the new IA** (an experiment's samples)
shows nothing, for *all* data, not just migrated data. The grouping-review banner only
appears when flagged loads exist, so a migrated experiment gets a fully empty tab above
a StatBar that reads "31 SAMPLES".

Verified DOM: `<div data-testid="corpus-sheet-slot" class="text-sm text-ink-soft"></div>`
(empty).

---

## #3 — P2 — Experiments home page has no nav chrome

`/experiments` (`ExperimentsHomePage`) renders a bare `PageFrame` — no wordmark, no
Experiments | Series nav. It is a dead-end: the only ways out are clicking an experiment
card or browser-back. `ExperimentShell` (the `/experiments/:id` layout) *does* render its
own chrome; the home page and `/experiments/new` do not.

This is one symptom of the larger nav architecture: the branch **added** a parallel
`ExperimentShell` + `ExperimentTopNav` rather than evolving `CorpusShell` /
`CorpusTopbar`, yielding three different top-chromes:

| Route | Shell | Top nav |
|-------|-------|---------|
| `/experiments`, `/experiments/new` | none (bare `PageFrame`) | **none** |
| `/experiments/:id/*` | `ExperimentShell` | `ExperimentTopNav` + own header + Corpus\|Config tabs |
| `/samples`, `/series`, `/sample/:id` | `CorpusShell` | `CorpusTopbar` (with an added Experiments tab) |

The spec (§7/§9.6) justified placing the experiments tree outside `CorpusShell` so the
header would not stack on `CorpusTopbar`, but the fragmentation is real. **Design
direction: lean toward collapsing back to one shell** (single nav system, chrome-ful
home page). Deferred as a deliberate design decision.

---

## #4 — P3 — Home experiment cards are sparse vs mockup

Cards show only a status dot + name + `data_dir`. The mockup specified a date range,
counts, and an "N need grouping review" hint. The list endpoint (`useExperiments`)
doesn't carry stats — only the detail endpoint does — so counts can't render on the
gallery without widening the list response. Degraded, not broken.

---

## #5 — P2 — Grouping empty-state was misleading copy — **fixed**

`GroupingReviewPage` collapsed every empty case to `No samples match "${q}"`, even when
the filter was empty, so the common all-clear case rendered `No samples match ""` (empty
quotes) and read as an error. Fixed (`019602fc`) by splitting into three honest states:

- no loads at all → "No loads yet" + rescan hint
- "Needs review" filter, nothing flagged → "Nothing needs review" (the all-clear)
- active filter, no match → the literal `No samples match "<q>"`

Pinned by `packages/HimalayaUI/frontend/test/GroupingReviewPage.empty.test.tsx`.

## #6 — P3 — Copy nits — **fixed**

- The grouping filter placeholder used a double-hyphen (`Filter samples -- name…`) →
  reworded to `Filter samples by name or glob…`.
- The Configuration **Acquisition** card rendered a blank body when there is no exposure
  timeline → now shows an empty-state ("No acquisition timeline yet. Rescan…").

Both landed in `019602fc`.
