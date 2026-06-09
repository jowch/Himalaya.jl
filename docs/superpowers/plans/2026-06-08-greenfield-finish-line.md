# Greenfield UI Rebuild — Finish-Line Tracker

> **For agentic workers / the `/loop`:** this is the living checklist that drives the rebuild to completion. Each loop iteration: read this file + `git log --oneline -8`, pick the next unchecked unit, execute it subagent-driven, gate it, check it off, commit. **Definition of "finish line" = the branch is READY TO MERGE (all milestones ✅, all gates green, live-validated) — NOT merged.** The branch lands as one piece on the user's explicit say-so; never merge or run `finishing-a-development-branch` autonomously.

**Branch:** `worktree-greenfield-ui-rebuild` (worktree `.claude/worktrees/greenfield-ui-rebuild`). Stays UNMERGED.

**User profile / law:** the orchestrator's behavior is governed by the core memory `user_jonathan_core.md`. The non-negotiables below are extracted from it.

---

## Standing constraints (every unit, every subagent)
- **Cadence:** subagent-driven — fresh implementer per task → **spec-compliance review THEN code-quality review** (`frontend-reviewer`/`himalaya-reviewer`/`queue-reviewer`/`saxs-physics-reviewer` as fits) → fix loop. **Dispatch implementers SEQUENTIALLY** (git index lock). The **orchestrator (root session) runs the slow suites**; subagents run only targeted/foreground checks. TDD: failing test → minimal impl → verify → commit, each step its own commit.
- **No-legacy-port:** `src/print/**` imports NEW (`src/print/**`) + CARRIED (`queries`/`api`/`state`/`lib`/`hooks`) ONLY — never OLD (`src/components/**`,`src/pages/**`). Reimplement from print components + mockup; reuse logic not presentation. Never reintroduce an affordance "we already agreed wouldn't be in the new frontend."
- **Design:** closed-look/open-placement — appearance in `ui/`/composites, placement-only `className`. Controls-don't-lie. Honesty over fake states. Accent (terracotta hue 38) is interaction-only — never let it leak onto auto-peaks/borders, never let dark-era tokens re-enter. Colour=meaning (measured line neutral ink). Extensible bases, not hardcoded variants.
- **Git:** NEVER `git add -A` — stage only named files. Every commit ends `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Never stage `registry.ts`/`*.bones.json` except a sanctioned boneyard capture, nor stray `.js`.
- **Gate (must be green before a unit is ✅):** `npx tsc -p tsconfig.build.json --noEmit` (0) · `npm run lint:design` (0) · targeted `npm test` (orchestrator runs FULL `npm test`) · `npm run e2e` (affected specs) · `npm run build`. Clean stray `.js` before trusting Vite/Vitest ([[feedback_greenfield_stray_js_shadow]]).
- **Verify by rendering, not reading** for any visual change — screenshot-compare mockup ⇄ live ⇄ Storybook (live-audit harness recipe: serve a writable dev-DB COPY from Julia source + Vite + Playwright MCP — see `feedback_live_audit_harness`). Production (:8080) untouched.
- **When a genuine PRODUCT decision appears** (not an execution detail), STOP and ask the user (`AskUserQuestion`) — don't guess product direction. Mid-unit discoveries that are out of scope → note them here as a new unit, don't scope-creep.

---

## Milestones (work top-down; check off in the same commit that lands the unit)

### M-A — New-series creation flow (the missing primary flow) ✅ (2026-06-09)
**DONE + live-verified end-to-end:** corpus checkbox picker (CheckCircle) → ComposeBar "+ New series" → `location.state` seed → scoping (warm propose / cold assign) → "Confirm & build" = two sequential queue ops (tag write `scopeSeries` → single-write `createSeries`) → navigate to `/series/:id` on `createSeries.data.id`. Backend `POST /api/series` now resolves the plate from the recipe on create (`_resolve_series_plate!`), so the new builder shows traces immediately. Keystone saga (3 queue bugs: idempotency op-id collision → compound-op deferred race → SSE-synthesized-partial nav guard) — see [[feedback_queue_sse_wins_partial]]. Commits: 22ce155·1448bba·f194373·01c26a8·4b34440·22cea09·dfcbd92·9c47272·c0f9915·0ded204·73c01124·17b1c0b·7e68103 + plan ca78609.
**Why:** there is currently NO way to create a series from the UI (`POST /api/series` exists + `api.saveSeries(no id)` wraps it, but nothing calls it; the CullBar offers only Drop/Restore/Clear; scoping's "Confirm & build" only writes ordering tags). Highest-value gap. NOTE: the machine proposal (`proposeOrdering`) only fires when an ordering variable already exists in sample tags/names — for an unstructured corpus it proposes nothing, so a **sample picker on the corpus page is the necessary path, not just an accelerator**.
**Design (LOCKED 2026-06-08 with user):** B-shaped flow with a corpus-page **sample picker**:
1. **Corpus page sample picker** — add a NEW **left-most column with a checkbox (checkmark)** to the SheetTable, a **sample-grain** selection DISTINCT from the existing frame-level cull (thumbnails keep their ⇧click→cull + dark CullBar). Selecting sample rows raises a **compose bar**: "N samples · + New series · Clear". (Closed-look: the column + compose bar are print primitives; honest — compose bar only shows when ≥1 sample checked.)
2. **Carry selection → scoping** (Zustand or route state) on "+ New series" → navigate `/series/new` seeded with the selected sample ids (scoping reads the seed, not the whole corpus).
3. **Scoping** proposes the ordering variable from the seeded samples; **cold/unstructured case** (no proposable variable) gets an honest "name the ordering variable / assign values" path instead of the dead-end empty state.
4. **"Confirm & build" CREATES the series** — writes the ordering tags AND `POST /api/series` with the recipe (selected samples + ordering) → navigate to the new `/series/:id` builder.
- Plan-first: write the M-A implementation plan, review it, THEN execute subagent-driven.
- Acceptance: from `/samples`, check sample rows → "+ New series" → (scoping confirm/assign) → land on a real new `/series/:id` builder with those members; e2e covers the full path; live-validated in the browser. Works for BOTH a name-encoded corpus (machine proposes) AND an unstructured one (user assigns).

### M-B — Builder rail full-height + builder polish ✅ (2026-06-09)
**DONE + live-verified.** Commits `ffa50e3` (items 1–4) + `d38e6ae` (item 5). Reviewer (`frontend-reviewer`) clean — no Critical/Important. Gates: tsc 0 · lint:design 0 · targeted vitest (BuilderRail/SeriesBuilderPage/SeriesCard/folio + route/autogroup/adapters) green.
- (1) **Rail full-height** ✅ — `SeriesBuilderPage` replaced `<PageFrame width="builder">` + `grid-cols-[1fr_336px] items-start` with `<div data-testid="builder-workspace" className="grid min-h-full grid-cols-1 xl:grid-cols-[minmax(0,1fr)_336px]">`; the `BuilderRail` is now a **direct grid child**, so default `align-items:stretch` stretches it to the work column. **Live-measured: rail 703px == work column 703px, byte-identical span (top 56 → bottom 759), `railIsDirectChild:true`.** Structurally equivalent to FocusPage's `focus-workspace` (also a bare auto-height grid in a `<Skeleton>` with no height context — so both share the same short-page gap; matching it IS the goal). The old `items-start` was what pinned the rail to 625px content height. Plate keeps 1180px via inner `mx-auto max-w-[1180px]`. (`min-h-full` is an inert no-op given the auto-height ancestor chain — harmless, left in place; does NOT diverge from Focus.)
- (2) **Adjust hidden in draft** ✅ — page passes `onAdjust` only when no draft live (`{...(liveDraft ? {} : { onAdjust: ensureDraft })}`, exactOptionalPropertyTypes-clean conditional spread); BuilderRail renders Adjust only when `onAdjust` provided. **Live-clicked: read state = Confirm(disabled)+Adjust; draft state = Confirm(enabled)+Cancel, Adjust GONE.**
- (3) **Single q-scale toggle** ✅ — kept the plate's SegmentedControl, removed the rail DISPLAY one; `scale`/`setScale` still owned by `SeriesBuilderPage` and threaded to the plate (not orphaned). **Live-confirmed: one `log q / linear q` group on the plate toolbar, rail DISPLAY has only the Trace-offset slider.**
- (4) **Annotation toggles gate on peaks** ✅ — `seriesHasPeaks = rows.some(r => r.anchors.length > 0)`; `AnnotationToggleRow` `enabled` → buttons `disabled` + never `armed` when false. **Live-confirmed: peakless 2-sample series shows both Peak-ticks/Peak-labels `[disabled]`.**
- (5) **Dangling "· by"** ✅ — `SeriesCard` gates `· by {variable}` on `variable.trim() !== ""`. **Live-confirmed on folio: "Bundle A (server-edited)" renders "2 samples" with no dangling "by"; variable-bearing cards show "· by temperature".**

### M-C — Drag-to-reorder for builder recipe rows ✅ (2026-06-09)
**DONE + live-verified.** Commit `f800712`. Reviewer (`frontend-reviewer`) clean — no issues. Wired `useDragReorder(onReorder)` into `RecipeEditor` (`SeriesBuilderPage.tsx`), mirroring the scoping-page pattern byte-for-byte: `dragItemProps(i)` spread + `group relative cursor-grab` + `data-dragging`→`opacity-50` + edge-gated `bg-accent` drop-line + `<GripHandle/>` (aria-hidden) as first child. The ▲▼ Move-up/Move-down IconButtons stay UNCHANGED as the keyboard path (GripHandle's aria-hidden contract requires a keyboard reorder path). Same `onReorder` → `reorderSeriesSample` (one Zustand action, two entry points — drag is additive). **Fixed a controls-don't-lie defect: the rail label "Traces — drag to reorder" was previously lying (drag unimplemented).** Gates: tsc 0 · lint:design 0 · SeriesBuilderPage 18/18 · useDragReorder 10/10. **Live drag verified in-browser: dragging row 1 onto row 0 reordered `[HEPES, LL37]` → `[LL37, HEPES]`; each row `draggable=true` + `cursor:grab` + carries grip + ▲▼.** (Note: faithful live DnD verification needs a frame-yield between dragstart/dragover/drop so React commits `draggingIndex` between ticks — a synchronous dispatch no-ops and falsely reads as broken.)

### M-D — Plan 7: dead-code sweep ✅ (2026-06-09)
**Full gate GREEN:** lint:design 0 · full `npm test` **2055/2055 passed (260 files)** · `npm run build` 0 · tsc 0. (Test count ~2186→2055 = 19 deleted test files + trims, expected.)
**Frontend tail DONE.** Verified-orphan sweep (every candidate grepped for live importers, incl. `import type`, before `git rm`; tsc the arbiter after each batch). Commits `67ad992` (18 components) · `7094f5e` (19 orphan-only tests) · `87c2535` (3 shared-test trims) · `ab51a02` (lib-tail). Net ~−5800 lines.
- **Deleted 18 components:** AutogroupCard, FigureExportControls, GroupingModeToggle, MillerPlot, NoUsableExposureNotice, OffsetSlider, Pn3mIcon, SamplePickerRow, ScaleToggle, SeriesPhaseStrip, AnnotationToggles, MemberTraceLayer; clusters SeriesTrackingOverlay+PlotSurface, BandResizeDivider+ActiveBandContext, MemberMetaRow+RowActionZone. + lib-tail `lib/series/trackGeometry.ts`+`anchors.ts` (only consumer was the deleted SeriesTrackingOverlay).
- **KEEP (verified live, NOT deleted):** `RepresentationToggle` (figure-export type-imports `Representation` — the Explore "no-JSX-import" pass falsely flagged it; caught pre-delete), `MemberHeatmapLayer`+`CrossTraceTrackingLayer` (figure-export `buildMemberHeatmapMarks`/`buildCrossTraceTrackingMarks`), DetectorImage/ThumbnailGallery/AssignmentCart/ResolvingFallback (live in print pages/routing).
- **Tests:** 19 orphan-only deleted; 3 surgically trimmed to keep live coverage (`SegmentedControl.test` −GroupingModeToggle; `labelDodge.test` kept live `layoutPeakLabels`, −MemberTraceLayer cases; `MemberHeatmapLayer.test` kept, comment-only fix). `src/pages/` already gone (tracker note was stale).
- **Backend dead-code = OUT OF SCOPE (correctness, not avoidance):** the `comparison*` tables + `apply_event!` dispatcher branches are RETAINED BY DESIGN for historical event replay (per AppRoutes comment) — no safe backend deletion exists.
- **PENDING:** orchestrator full `npm test` + `npm run build` + `lint:design` (running, bg `brwf2lj8h`) → flip to ✅ when green. tsc already 0; targeted suites (SegmentedControl/MemberHeatmapLayer/labelDodge 38/38) green.

### M-E — Plan 8: shell reskin ⬜
- The carried `CorpusShell` / app-shell reskin to the greenfield look (the last legacy presentation under the per-route bodies). Plan-first; mockup-driven.

### M-F — Deferred polish ⬜
- Boneyard skeleton captures for the greenfield pages (sanctioned capture recipe; `registry.ts` IS committable for this) + a per-page live visual-fidelity pass vs the mockups.
- Decide with the user: enforce the e2e suite. **Actual pre-existing red count (measured 2026-06-09 M-A run): 10 failing tests across 4 specs** — `permalinks` (1: legacy /index deep-link), `figure-export` (2: PNG/SVG download), `qlink` (3: ring↔row hover), `smoke` (4: new-user tutorial, curate, reanalyze ×2). Assignment-model + slug rot from #280; NOT M-A regressions (new-series-creation 4/4 green). Fix-or-quarantine + CI-gate, or leave for post-merge.
- NOTE (2026-06-09): the 2 long-"known" `bonnet_lattice` Julia dev-link errors are GONE in this worktree once core is `Pkg.develop`-linked (done this session); the full HimalayaUI suite is fully green.

### FINAL — whole-branch review, then STOP ⬜
- Run the full gate one last time (tsc · lint:design · full `npm test` · `npm run e2e` · `npm run build`) + a final `frontend-reviewer` over the whole branch + a live click-through of every page.
- Update `user_jonathan_core.md` / `project_greenfield_composite_layer.md` to mark the rebuild ready.
- **STOP and report "FINISH LINE: branch ready to merge".** Do NOT merge, do NOT run `finishing-a-development-branch` — wait for the user.

---

## Open decision — RESOLVED 2026-06-08
**How is a new series born?** LOCKED: B-shaped with a corpus-page sample picker (a new left-most checkbox column at the sample grain) → "+ New series" carries the selection → scoping proposes-or-assigns the ordering variable → "Confirm & build" creates the series (`POST /api/series`) → builder. See M-A above. (Rationale: pure machine-proposal dead-ends on unstructured corpora; the sample-grain picker is the necessary feeder. Frame-level cull stays unchanged on the thumbnails.)

---

## The loop prompt (dynamic `/loop`)
Run `/loop` (no interval → dynamic self-pacing) with the prompt below verbatim.

```
Drive the greenfield UI rebuild to the finish line, one unit per iteration, governed by
docs/superpowers/plans/2026-06-08-greenfield-finish-line.md (the tracker) and the core
memory user_jonathan_core.md. Work on branch worktree-greenfield-ui-rebuild; it STAYS UNMERGED.

EACH ITERATION:
1. Orient: read git log --oneline -8 and the tracker's Milestones checklist. Clean any stray
   src/**/*.js first. Figure out what the last iteration finished so you never redo it.
2. Pick the next UNCHECKED unit (top-down: M-A → M-B → … → FINAL). Honor the tracker's
   "Open decision" — if M-A's create-flow is not yet locked, STOP and AskUserQuestion; do not
   guess product direction.
3. If the unit needs a plan (M-A, M-E, or any non-trivial unit), FIRST write/refine a short
   plan and have a project reviewer sanity-check it; then execute. Trivial layout/wiring units
   (M-B, M-C) can go straight to implementation.
4. Execute subagent-driven: fresh implementer (sequential dispatch — git index lock) → spec-
   compliance review → code-quality review (frontend-reviewer / himalaya-reviewer / etc.) →
   fix loop. TDD. Honor every Standing Constraint in the tracker (no-legacy-port, closed-look,
   never git add -A, commit footer, etc.).
5. Gate (orchestrator runs the slow ones): tsc -p tsconfig.build.json --noEmit · npm run
   lint:design · the relevant npm test · npm run e2e (affected) · npm run build. For any VISUAL
   change, verify by RENDERING via the live-audit harness (writable dev-DB copy, Playwright MCP),
   not by reading code.
6. Check the unit off in the tracker and commit (named files only, co-author footer). Update
   the memory tracker (project_greenfield_composite_layer.md) when a milestone completes.
7. Re-assess: if a genuine product decision surfaced, AskUserQuestion and pause the loop (don't
   keep looping on an unanswered decision). If the same failure recurs 3×, AskUserQuestion.

STOP CONDITIONS (do NOT schedule another wakeup):
- FINAL reached → run the whole gate + a final whole-branch review + a live click-through →
  report "FINISH LINE: branch ready to merge" and STOP. Do NOT merge / do NOT run
  finishing-a-development-branch — that's the user's call.
- A product decision or a 3×-stuck unit needs the user.

PACING: dynamic. When you've dispatched harness-tracked work (subagents, a backgrounded
npm test) you'll be re-invoked automatically on completion — schedule a LONG fallback
(1200s+) rather than polling. When genuinely idle between units, a short re-check is fine.
Keep iterations focused so the prompt cache stays warm.
```

**To stop early:** press `Esc` while the loop is waiting, or send a message. The loop also self-exits at FINAL or when it needs a decision.
