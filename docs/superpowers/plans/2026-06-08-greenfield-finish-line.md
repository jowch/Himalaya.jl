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

### M-A — New-series creation flow (the missing primary flow) ⬜
**Why:** there is currently NO way to create a series from the UI (`POST /api/series` exists + `api.saveSeries(no id)` wraps it, but nothing calls it; the CullBar offers only Drop/Restore/Clear; scoping's "Confirm & build" only writes ordering tags). Highest-value gap. NOTE: the machine proposal (`proposeOrdering`) only fires when an ordering variable already exists in sample tags/names — for an unstructured corpus it proposes nothing, so a **sample picker on the corpus page is the necessary path, not just an accelerator**.
**Design (LOCKED 2026-06-08 with user):** B-shaped flow with a corpus-page **sample picker**:
1. **Corpus page sample picker** — add a NEW **left-most column with a checkbox (checkmark)** to the SheetTable, a **sample-grain** selection DISTINCT from the existing frame-level cull (thumbnails keep their ⇧click→cull + dark CullBar). Selecting sample rows raises a **compose bar**: "N samples · + New series · Clear". (Closed-look: the column + compose bar are print primitives; honest — compose bar only shows when ≥1 sample checked.)
2. **Carry selection → scoping** (Zustand or route state) on "+ New series" → navigate `/series/new` seeded with the selected sample ids (scoping reads the seed, not the whole corpus).
3. **Scoping** proposes the ordering variable from the seeded samples; **cold/unstructured case** (no proposable variable) gets an honest "name the ordering variable / assign values" path instead of the dead-end empty state.
4. **"Confirm & build" CREATES the series** — writes the ordering tags AND `POST /api/series` with the recipe (selected samples + ordering) → navigate to the new `/series/:id` builder.
- Plan-first: write the M-A implementation plan, review it, THEN execute subagent-driven.
- Acceptance: from `/samples`, check sample rows → "+ New series" → (scoping confirm/assign) → land on a real new `/series/:id` builder with those members; e2e covers the full path; live-validated in the browser. Works for BOTH a name-encoded corpus (machine proposes) AND an unstructured one (user assigns).

### M-B — Builder rail full-height + builder polish ⬜
- Builder `BuilderRail` should stretch to full page height flush under the header (match the Focus page's grid), not sit as a content-height grey box (measured live: builder 625px vs Focus 1481px). Layout-only fix in `SeriesBuilderPage`/`PageFrame`.
- Minor polish surfaced in the 2026-06-08 live audit: (a) "Adjust" should not linger enabled in draft state alongside Cancel (BuilderRail bakes the Confirm+Adjust pair — withhold/hide Adjust when a draft is live); (b) two redundant log/linear-q toggles on the builder (plate + rail DISPLAY) — keep one; (c) annotation toggles default both-armed on a peakless form-factor series — default off or gate on peak availability; (d) folio card meta "2 samples · by" dangling author.

### M-C — Drag-to-reorder for builder recipe rows ⬜
- `useDragReorder` already exists and is wired on the scoping page; the builder hand-rolled ▲▼ buttons. Wire `useDragReorder` into the builder's recipe rows (mouse drag) while KEEPING the ▲▼/aria path for keyboard a11y (the `GripHandle` is `aria-hidden` by design — drag-only would regress keyboard users). Reuse the existing hook, don't reinvent.

### M-D — Plan 7: dead-code sweep ⬜
- **Frontend tail** (left by the Series-builder cutover): legacy plot helper-layers (`MemberTraceLayer`/`MemberHeatmapLayer`/`PlotSurface`/`CrossTraceTrackingLayer`/`SeriesTrackingOverlay`/`BandResizeDivider`/etc.) + now-orphaned carried controls (`GroupingModeToggle`/`AnnotationToggles`/`FigureExportControls`/`SeriesPhaseStrip`) — grep each for ANY surviving consumer (incl. tests/stories) before `git rm`.
- **Backend dead-code** track per the original Phase-4 cutover strategy.
- **Orphaned by M-A:** `useScopeSeries` (queries.ts) lost its last caller when M-A Task 7 swapped scoping's build to `useScopeAndCreateSeries` — delete it, and check whether `scopeSeriesMutator` is now dead too (it may still be referenced by `mutatorRegistry` for foreign `add_tag` replay — grep before removing).
- Acceptance: full `npm test` + `tsc` green after each deletion (the suite is the proof a deletion was safe).

### M-E — Plan 8: shell reskin ⬜
- The carried `CorpusShell` / app-shell reskin to the greenfield look (the last legacy presentation under the per-route bodies). Plan-first; mockup-driven.

### M-F — Deferred polish ⬜
- Boneyard skeleton captures for the greenfield pages (sanctioned capture recipe; `registry.ts` IS committable for this) + a per-page live visual-fidelity pass vs the mockups.
- Decide with the user: enforce the e2e suite (the 4 pre-existing red smoke specs — assignment-model rot from #280) — fix-or-quarantine + CI-gate, or leave for post-merge.

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
