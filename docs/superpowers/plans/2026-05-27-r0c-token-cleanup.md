# R0c — Residual per-component token cleanup Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:executing-plans. Steps use checkbox (`- [ ]`).

**Goal:** Sweep the components that carry semantically-wrong tokens the R0a value-remap could not auto-fix, and remove dead theme scaffolding the keystone left behind, so no dark-theme/ice-blue chrome remains on any Print surface.

**Architecture:** Pure token/class swaps to the canonical Print tokens (DESIGN.md `colors:`), one drift-correction to `--color-success`, and deletion of dead `MutationObserver`/comment scaffolding. No new UI/affordances.

**Tech Stack:** React + TypeScript + Tailwind v4 (`@theme` CSS vars), Vitest.

---

## Bucket A — semantically-wrong tokens

### Task 1: Kept-dot sage drift (T-4)
DESIGN.md status block: `success: oklch(0.520 0.120 162)` (sage) and explicitly notes the kept verdict dot is sage (T-4). styles.css `--color-success` drifted to `oklch(0.72 0.08 155)` (mint). Correct the token to its DESIGN.md value; the kept-dot keeps `bg-success` (now correctly sage).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/styles.css` (`--color-success`)
- Test: `packages/HimalayaUI/frontend/test/loupeKeptDot.test.tsx` (new) + a styles assertion

- [ ] Step 1: Add a render test asserting LoupeSidebar kept-dot uses `bg-success` (and NOT `bg-accent`) when not rejected, and the styles `--color-success` value equals DESIGN.md sage.
- [ ] Step 2: Run, expect FAIL (success token value mismatch).
- [ ] Step 3: Set `--color-success: oklch(0.520 0.120 162);` in styles.css.
- [ ] Step 4: Run, expect PASS.
- [ ] Step 5: Commit.

### Task 2: ThumbnailGallery selection chrome (T-7)
Mockup `.thumb`: bg `frame-edge`; hover ring = `hair-strong`; selected ring = `accent` (terracotta); indexing badge = accent bg + light text. Shipped uses `ring-accent`/`bg-accent/85`/`text-bg` (the `text-bg`=paper utility is fragile + the hover ring uses `fg-muted`). Swap to explicit Print tokens.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ThumbnailGallery.tsx`

- [ ] Step 1: Selected ring → `ring-2 ring-print-accent shadow-sm shadow-print-accent/30`; hover ring → `ring-1 ring-hair/60 group-hover:ring-2 group-hover:ring-hair-strong`; indexing badge → `bg-print-accent/90 text-paper`.
- [ ] Step 2: `npm run build`. Commit.

### Task 3: LoupeFrame + FocusDetectorPanel detector frame + caption (T-8)
The detector frame is a dark window (`frame-edge`), not paper. After R0a `bg-bg` resolves to light paper → dark caption invisible. Mockup: `.big-frame` bg `frame-edge`; `.frame-tag` light mono (`oklch(0.82 0.01 80)`); dropped-tag near-white on accent.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/LoupeFrame.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/FocusDetectorPanel.tsx`
- Modify: `packages/HimalayaUI/frontend/src/styles.css` (add `--color-frame-tag` light token for captions on the dark frame)

- [ ] Step 1: Add `--color-frame-tag: oklch(0.82 0.01 80);` to styles.css `@theme`.
- [ ] Step 2: LoupeFrame frame wrapper: `border-frame-edge bg-frame-edge`; caption → `text-frame-tag`; dropped badge → `bg-print-accent text-paper`.
- [ ] Step 3: FocusDetectorPanel inner frame: `border-frame-edge bg-frame-edge`.
- [ ] Step 4: `npm run build`. Commit.

### Task 4: MultiTracePlot / MemberMetaGutter chrome (B-C)
Render-core tooltip + selection + insertion-line. After R0a these auto-remap to Print values, but switch to explicit Print token names matching DESIGN.md (tooltip plate/hair/ink; selection + insertion accent → print-accent).

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/MultiTracePlot.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/MemberMetaGutter.tsx`

- [ ] Step 1: MultiTracePlot tooltip `border-border bg-bg-elevated text-fg`/`text-fg-muted` → `border-hair bg-plate text-ink`/`text-ink-soft`; band overlay `ring-accent`/`bg-accent/*` → `ring-print-accent`/`bg-print-accent/*`.
- [ ] Step 2: MemberMetaGutter hover ring `ring-accent/50` → `ring-print-accent/50`; drop-indicator `bg-accent` → `bg-print-accent`.
- [ ] Step 3: `npm run build`. Commit.

## Bucket B — dead theme scaffolding

### Task 5: Remove dead DetectorImage MutationObserver
The `<html>` theme-class toggle producer was removed by R0a; the observer + its comment are dead.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/DetectorImage.tsx`

- [ ] Step 1: Delete the `useEffect` MutationObserver (and the now-false comment).
- [ ] Step 2: `npm run build`. Commit.

### Task 6: Fix stale CorpusTopbar comment
- [ ] Step 1: Rewrite the "T theme shortcut still toggles" sentence to reflect R0a retirement.
- [ ] Step 2: Commit.

### Task 7: Tidy inert `theme:"dark"` seeds + Hexagonal comment
KEEP migration-exercising seeds (persistMigrate/mergePersistedState tests; useGlobalShortcuts assertion; e2e `version:3` blobs). Tidy inert unit seeds + non-migration e2e seeds. Align phases.ts Hexagonal comment to "Rose".

**Files:**
- Modify (unit, remove `theme:"dark"`): FocusDetectorPanel.test, useSyncActiveSampleFromRoute.test, seriesFolioRoute.test, seriesBuilderRoute.test, FocusWorkspacePage.routing.test, FocusWorkspacePage.layout.test, smoke.test
- e2e: LEFT UNTOUCHED. On inspection EVERY e2e theme seed lives in a
  `version: 3` localStorage blob, i.e. each deliberately exercises the
  v3→v5 migrate-strip. Per the kit ("KEEP those that exercise the migrate-
  strip"), all are retained: smoke/permalinks/figure-export/multiplayer +
  corpus-culling/series-scoping/qlink/series-builder.
- KEEP: persistMigrate.test, mergePersistedState.test, useGlobalShortcuts.test
- Modify: phases.ts Hexagonal comment `magenta-rose` → `rose`

- [ ] Step 1: Remove the inert seeds; align the comment.
- [ ] Step 2: `npm test` + `npm run build`. Commit.

## Final
- [ ] Grep sweep clean: no `bg-success` misuse, no ice-blue `ring-accent`/`bg-accent`/`text-bg` selection chrome on swept components, no dead scaffolding.
- [ ] `npm run build` + `npm test` green.
