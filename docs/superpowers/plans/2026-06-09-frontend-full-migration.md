# Full frontend migration — eliminate the old `src/components/` tree

> Goal: the ENTIRE frontend lives under `src/print/`; `src/components/` is deleted. Remove the vestigial old `src/components/ui/` design-system barrel + the old duplicate components, then relocate the live shell into `src/print/`. User-chosen scope (2026-06-09): full relocation.

**Branch:** `worktree-greenfield-ui-rebuild` (stays UNMERGED). Work from `packages/HimalayaUI/frontend/`.

**Arbiters after every stage:** `npx tsc -p tsconfig.build.json --noEmit` (0) · `npm run lint:design` (0) · targeted `npm test` · (orchestrator runs full `npm test` + `npm run build` at the end). Stage the work; each stage its own commit(s). NEVER `git add -A`; commit footer `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`.

## Load-bearing gotchas (verified 2026-06-09)
- **`scripts/check-design.mjs` is path-coupled.** `COLOR_AUTHORING_ALLOW` hard-codes `components/MemberHeatmapLayer.tsx`, `components/DetectorImage.tsx`, `components/FocusDetectorPanel.tsx`, `main.tsx`, `phases.ts`, `lib/comparison/coloring.ts`; renderer-exempt PREFIXES are `components/ui/`, `print/ui/`, `print/plot/`, `print/detector/`, `print/comb/`, `print/export/`; and the **import-boundary guard** flags any `src/print/**` file that relatively imports from `components/` or `pages/`. When a file moves, its allowlist entry must move with it (or it lands under a `print/*` exempt prefix), AND the exclude comment at the top should be updated. After EVERY move, run `npm run lint:design`.
- **IconButton `tone="accent"`** is used at `NavModal.tsx:324`. Confirm `print/ui/IconButton` `IconButtonTone` includes `"accent"` before swapping; if not, add it to the print primitive (don't silently downgrade the tone). tsc will flag a missing union member.
- **tsc is the relocation arbiter:** a missed import-path update → TS2307. Run tsc after each move batch; fix every break before committing.
- Old + greenfield DUPLICATES exist for AssignmentCart (`print/components/`), ThumbnailGallery (`print/components/`), DetectorImage (`print/detector/`). The live pages already import the greenfield ones; only TESTS import the old ones. Repoint those tests to the greenfield path (or delete the old test if a greenfield test already covers it — verify before deleting).

---

## Stage A — repoint old-barrel imports → print/ui (pure swaps)
1. `src/App.tsx:7` — `ToastContainer` from `./components/ui` → `./print/ui`.
2. `src/lib/figure-export/marks/traceExportMarks.ts:8` — `peakMark` from `../../../components/ui/peakMark` → `../../../print/ui/peakMark`.
3. `src/lib/figure-export/marks/multiTraceExportMarks.ts:15` — `peakMark` same swap.
Gate: tsc 0 · lint:design 0 · `npm test -- Toast figure-export`. Commit.

## Stage B — repoint the live shell off the old barrel
4. `src/components/NavModal.tsx:7` — `{ IconButton, ModalShell }` from `./ui` → `../print/ui`. VERIFY `tone="accent"` (line 324) compiles against print IconButton; if `"accent"` isn't in the union, ADD it to `print/ui/IconButton` (match the old tone's class) rather than changing NavModal's tone.
5. `src/components/OnboardingFlow.tsx:5` — `{ Button, ModalShell, Kicker }` from `./ui` → `../print/ui`. Verify Button/Kicker APIs match (variant/tone props).
Gate: tsc 0 · lint:design 0 · `npm test -- NavModal OnboardingFlow`. Commit.

## Stage C — delete the old duplicate components + vestigial RepresentationToggle
6. Repoint any TEST importing old `components/AssignmentCart`/`ThumbnailGallery`/`DetectorImage` to the greenfield path (`print/components/AssignmentCart`, `print/components/ThumbnailGallery`, `print/detector/DetectorImage`). If the old test duplicates a greenfield test, delete the old test (verify the greenfield test exists + covers it).
7. `git rm src/components/AssignmentCart.tsx src/components/ThumbnailGallery.tsx src/components/DetectorImage.tsx` (delete ThumbnailGallery before DetectorImage — DetectorImage's only non-test importer is old ThumbnailGallery).
8. `RepresentationToggle.tsx`: its `Representation` type is the only thing imported live (`multiTraceExportMarks.ts:17`, type-only) and it has a dead `./ui/SegmentedControl` import. Inline `export type Representation = "waterfall" | "heatmap"` (confirm the exact union) into `multiTraceExportMarks.ts` (and `multiTraceAdapter.ts` if it imports it too), then `git rm src/components/RepresentationToggle.tsx`.
9. Update `scripts/check-design.mjs`: remove the now-dead `components/DetectorImage.tsx` (and `components/FocusDetectorPanel.tsx` if already absent) from `COLOR_AUTHORING_ALLOW`.
Gate: tsc 0 · lint:design 0 · `npm test` (the affected suites). Commit per logical group.

## Stage D — delete the orphaned old barrel
10. Confirm zero importers: `grep -rn "components/ui" src test` returns nothing (ignore `.stories`/comments). Then `git rm -r src/components/ui/` (15 files).
11. Update `check-design.mjs`: remove `components/ui/` from the renderer-exempt prefix list + the top exclude comment (the dir no longer exists).
Gate: tsc 0 · lint:design 0 · `npm test`. Commit.

## Stage E — relocate the live shell into src/print/, delete src/components/
Remaining src/components files after A–D: AppRoutes, CorpusShell, CorpusTopbar, IndexSlugRedirect, InfrastructureBanner, NavModal, OnboardingFlow, ResolvingFallback, StaleUrlPage + the two export utilities MemberHeatmapLayer, CrossTraceTrackingLayer.
12. Create `src/print/shell/` and `git mv` the 9 shell files there. Move the 2 export utilities to `src/print/export/` (already a design-guard-exempt prefix — keeps their pixel-painting appearance allowed) OR `src/print/lib/`; if `print/export/` chosen, their COLOR_AUTHORING needs are covered by the prefix (verify MemberHeatmapLayer authors color → `print/export/` prefix exemption covers it; update `COLOR_AUTHORING_ALLOW` to drop the old `components/MemberHeatmapLayer.tsx` entry).
13. Update EVERY import path: `src/App.tsx` (imports AppRoutes/OnboardingFlow/NavModal/InfrastructureBanner/ToastContainer), the shell files' imports of EACH OTHER, `src/lib/figure-export/marks/multiTraceExportMarks.ts` (MemberHeatmapLayer/CrossTraceTrackingLayer paths), and ALL test/story files referencing the moved files. tsc is the arbiter — repeat until 0.
14. **Import-boundary guard:** the relocated shell now lives under `src/print/` and the guard forbids `src/print/**` importing `components/`/`pages/`. Since src/components is being deleted this resolves naturally, BUT verify no relocated file still references an old path. Also confirm CorpusTopbar's existing `../print/ui` imports become `../ui` (now that it's inside print/) — adjust relative paths.
15. `git rm` the now-empty `src/components/` (should be nothing left). Update `check-design.mjs` import-boundary + exclude comments (`src/components/**` no longer exists; the guard's `components/` checks become dead — simplify or leave inert, but update the stale path references in `COLOR_AUTHORING_ALLOW`/comments).
Gate: tsc 0 · lint:design 0 · `npm test` · orchestrator runs full `npm test` + `npm run e2e` + `npm run build` + a live smoke render.

## Acceptance
- `src/components/` no longer exists; the entire frontend is under `src/print/` (+ `src/App.tsx`, `src/lib/`, `src/state.ts`, `src/queries.ts`, `src/api.ts`, `src/hooks/`).
- tsc 0 · lint:design 0 · full `npm test` green · `npm run e2e` 35/3/0 · build 0.
- Live smoke: every page still renders (shell, samples, focus, folio, builder, loupe) with 0 console errors.
- `check-design.mjs` has no stale path references; `lint:design` still guards appearance + the import boundary (now within print/).
