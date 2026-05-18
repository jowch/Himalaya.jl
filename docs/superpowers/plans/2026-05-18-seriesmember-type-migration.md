# `SeriesMember` Type Migration Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Introduce a standalone `SeriesMember` TypeScript type and retype the multi-trace render pipeline onto it, with no change to render behaviour.

**Architecture:** `SeriesMember` is a full, self-standing interface in `api.ts` mirroring the I2.2 series route response (`fetch_series_with_plate`) — field-for-field identical to `ComparisonMember` except `comparison_id` → `series_id`. The render pipeline never reads the member's foreign key, so the nine source modules migrate by a pure type-name rename. `Compare.tsx` and `draftFactories.ts` keep feeding comparison data into the retyped pipeline via a localized bridge that I3.6 deletes wholesale.

**Tech Stack:** TypeScript (`exactOptionalPropertyTypes` on), React, Vitest. Frontend root: `packages/HimalayaUI/frontend/`.

---

## Conventions for every task

- **All commands run from `packages/HimalayaUI/frontend/`.**
- **Typecheck:** `npx tsc --noEmit` (uses `tsconfig.json`, which includes `test/`).
- **Unit tests:** `npm test -- <pattern>` runs Vitest one-shot (the repo's settings hook supplies `--run`).
- **`tsc --noEmit` is RED from Task 2 through Task 9 — this is expected.** A type migration cannot keep `tsc` green between every step: the moment a component is retyped, its not-yet-migrated consumers fail to typecheck. Per-task safety comes from **Vitest**, which strips types without checking them — a green Vitest run after each task proves render *behaviour* is unchanged. `tsc --noEmit` is the gate in **Task 10**, after every file is consistent.
- The nine source modules contain **zero references to `comparison_id`** (verified by grep). In those modules, replacing the identifier `ComparisonMember` with `SeriesMember` is the entire change.

## File map

Created:
- (none — `SeriesMember` is added to the existing `api.ts`)

Modified — source (`src/`):
- `api.ts` — add `SeriesMember`
- `lib/figure-export/adapters/multiTraceAdapter.ts`
- `lib/figure-export/marks/multiTraceExportMarks.ts`
- `lib/comparison/coloring.ts`
- `lib/comparison/labels.ts`
- `components/MemberTraceLayer.tsx`
- `components/MultiTracePlot.tsx`
- `components/MemberMetaRow.tsx`
- `components/MemberMetaGutter.tsx`
- `lib/comparison/draftFactories.ts` — `memberFromSaved` migrates; `fromComparison`/`fromComparisonAsFork` bridge
- `pages/Compare.tsx` — bridge (deleted at I3.6)

Modified — tests (`test/`), fixtures retyped to `SeriesMember`:
- `test/figure-export/multiTraceAdapter.test.ts`
- `test/coloring.test.ts`
- `test/MemberTraceLayer.test.tsx`, `test/AnnotationToggles.test.tsx`
- `test/MultiTracePlot.test.tsx`, `test/peakTooltip.test.tsx`
- `test/MemberMetaRow.test.tsx`, `test/MemberMetaRow.collapse.test.tsx`, `test/MemberMetaRow.drag.test.tsx`, `test/hoverPhaseColoring.test.tsx`
- `test/MemberMetaGutter.reorder.test.tsx`, `test/MemberMetaGutter.resize.test.tsx`, `test/MemberReorder.test.tsx`

Untouched (verified): `lib/comparison/yBands.ts`; `test/draftPersistence.test.ts`, `test/MentionChip.test.tsx`, `test/queue/saveComparison.test.tsx` (these reference `ComparisonMember` but do not import a migrated module — they stay on `ComparisonMember`).

---

## Task 1: Define `SeriesMember` in `api.ts`

**Files:**
- Modify: `src/api.ts` (insert after the `ComparisonMember` interface, which ends at the `}` on line 441; `Comparison` begins on line 443)

- [ ] **Step 1: Add the interface**

Insert this block immediately after the closing `}` of `export interface ComparisonMember` and before `export interface Comparison`:

```ts

/**
 * Per-member shape returned by the series GET / POST endpoints.
 * Mirrors `fetch_series_with_plate` (packages/HimalayaUI/src/series.jl).
 * Field-for-field identical to `ComparisonMember` except `comparison_id` →
 * `series_id`. This is the render pipeline's going-forward input type.
 */
export interface SeriesMember {
  id: number;
  series_id: number;
  exposure_id: number | null;
  display_order: number;
  band_height: number;
  y_offset: number;
  normalization: string;
  color_override: string | null;
  label_override: string | null;
  q_window_min: number | null;
  q_window_max: number | null;
  peak_display: unknown;
  snapshot: MemberSnapshot | null;
  is_stale: boolean;
  created_by: number | null;
  created_at: string | null;
}
```

- [ ] **Step 2: Verify the typecheck is still green**

Run: `npx tsc --noEmit`
Expected: PASS (the addition is purely additive — nothing consumes `SeriesMember` yet).

- [ ] **Step 3: Commit**

```bash
git add src/api.ts
git commit -m "Add the SeriesMember type (I3.2)"
```

---

## Task 2: Migrate the figure-export modules

**Files:**
- Modify: `src/lib/figure-export/adapters/multiTraceAdapter.ts` (`ComparisonMember` at lines 2, 13, 21, 121, 122, 124)
- Modify: `src/lib/figure-export/marks/multiTraceExportMarks.ts` (`ComparisonMember` at lines 5, 16)
- Test: `test/figure-export/multiTraceAdapter.test.ts`

- [ ] **Step 1: Rename in `multiTraceAdapter.ts`**

Replace every occurrence of the identifier `ComparisonMember` with `SeriesMember`. There are 6: the `import type` on line 2, and the type annotations on lines 13, 21, 121, 122, 124. After Step 1 line 2 reads:

```ts
import type { SeriesMember, Trace } from "../../../api";
```

- [ ] **Step 2: Rename in `multiTraceExportMarks.ts`**

Replace every occurrence of `ComparisonMember` with `SeriesMember` (2: import on line 5, annotation on line 16). Line 5 becomes:

```ts
import type { SeriesMember } from "../../../api";
```

- [ ] **Step 3: Retype the test fixture**

In `test/figure-export/multiTraceAdapter.test.ts`, replace every occurrence of `ComparisonMember` with `SeriesMember`, and in the member-fixture object(s) replace the `comparison_id` field with `series_id`.

- [ ] **Step 4: Run the affected Vitest file**

Run: `npm test -- multiTraceAdapter`
Expected: PASS — same tests, same assertions; the rename is type-only and Vitest erases types.

- [ ] **Step 5: Commit**

```bash
git add src/lib/figure-export/adapters/multiTraceAdapter.ts src/lib/figure-export/marks/multiTraceExportMarks.ts test/figure-export/multiTraceAdapter.test.ts
git commit -m "Migrate the figure-export multi-trace modules to SeriesMember (I3.2)"
```

---

## Task 3: Migrate `coloring.ts` and `labels.ts`

**Files:**
- Modify: `src/lib/comparison/coloring.ts` (`ComparisonMember` at lines 29, 73, 79, 89, 114, 151, 158)
- Modify: `src/lib/comparison/labels.ts` (`ComparisonMember` at lines 17, 20)
- Test: `test/coloring.test.ts`

- [ ] **Step 1: Rename in `coloring.ts`**

Replace every occurrence of the identifier `ComparisonMember` with `SeriesMember` (7 occurrences). Line 29 becomes:

```ts
import type { SeriesMember } from "../../api";
```

- [ ] **Step 2: Rename in `labels.ts`**

Replace every occurrence of `ComparisonMember` with `SeriesMember` (2 occurrences). Line 17 becomes:

```ts
import type { SeriesMember, Exposure, Sample } from "../../api";
```

- [ ] **Step 3: Retype the test fixture**

In `test/coloring.test.ts`, replace every occurrence of `ComparisonMember` with `SeriesMember`, and replace the `comparison_id` field with `series_id` in the member fixture(s).

- [ ] **Step 4: Run the affected Vitest file**

Run: `npm test -- coloring`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add src/lib/comparison/coloring.ts src/lib/comparison/labels.ts test/coloring.test.ts
git commit -m "Migrate comparison coloring + labels to SeriesMember (I3.2)"
```

---

## Task 4: Migrate `MemberTraceLayer.tsx`

**Files:**
- Modify: `src/components/MemberTraceLayer.tsx` (`ComparisonMember` at lines 36, 67, 118, 125)
- Test: `test/MemberTraceLayer.test.tsx`, `test/AnnotationToggles.test.tsx`

- [ ] **Step 1: Rename in `MemberTraceLayer.tsx`**

Replace every occurrence of the identifier `ComparisonMember` with `SeriesMember` (4 occurrences). Line 36 becomes:

```ts
import type { Trace, SeriesMember, MemberSnapshotPeak } from "../api";
```

(`Trace` and `MemberSnapshotPeak` are unchanged — only `ComparisonMember` is renamed.)

- [ ] **Step 2: Retype the test fixtures**

In `test/MemberTraceLayer.test.tsx` and `test/AnnotationToggles.test.tsx`, replace every occurrence of `ComparisonMember` with `SeriesMember`, and replace `comparison_id` with `series_id` in each member fixture.

- [ ] **Step 3: Run the affected Vitest files**

Run: `npm test -- MemberTraceLayer AnnotationToggles`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/components/MemberTraceLayer.tsx test/MemberTraceLayer.test.tsx test/AnnotationToggles.test.tsx
git commit -m "Migrate MemberTraceLayer to SeriesMember (I3.2)"
```

---

## Task 5: Migrate `MultiTracePlot.tsx`

**Files:**
- Modify: `src/components/MultiTracePlot.tsx` (`ComparisonMember` at lines 31, 128, 177)
- Test: `test/MultiTracePlot.test.tsx`, `test/peakTooltip.test.tsx`

- [ ] **Step 1: Rename in `MultiTracePlot.tsx`**

Replace every occurrence of the identifier `ComparisonMember` with `SeriesMember` (3 occurrences). Line 31 becomes:

```ts
import type { Trace, SeriesMember } from "../api";
```

- [ ] **Step 2: Retype the test fixtures**

In `test/MultiTracePlot.test.tsx` and `test/peakTooltip.test.tsx`, replace every occurrence of `ComparisonMember` with `SeriesMember`. In each `makeMember`-style fixture builder, replace the `comparison_id` field with `series_id` (e.g. `comparison_id: 100` → `series_id: 100`).

- [ ] **Step 3: Run the affected Vitest files**

Run: `npm test -- MultiTracePlot peakTooltip`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/components/MultiTracePlot.tsx test/MultiTracePlot.test.tsx test/peakTooltip.test.tsx
git commit -m "Migrate MultiTracePlot to SeriesMember (I3.2)"
```

---

## Task 6: Migrate `MemberMetaRow.tsx`

**Files:**
- Modify: `src/components/MemberMetaRow.tsx` (`ComparisonMember` at lines 22, 31)
- Test: `test/MemberMetaRow.test.tsx`, `test/MemberMetaRow.collapse.test.tsx`, `test/MemberMetaRow.drag.test.tsx`, `test/hoverPhaseColoring.test.tsx`

- [ ] **Step 1: Rename in `MemberMetaRow.tsx`**

Replace every occurrence of the identifier `ComparisonMember` with `SeriesMember` (2 occurrences). Line 22 becomes:

```ts
import type { SeriesMember } from "../api";
```

- [ ] **Step 2: Retype the test fixtures**

In `test/MemberMetaRow.test.tsx`, `test/MemberMetaRow.collapse.test.tsx`, `test/MemberMetaRow.drag.test.tsx`, and `test/hoverPhaseColoring.test.tsx`, replace every occurrence of `ComparisonMember` with `SeriesMember`, and replace `comparison_id` with `series_id` in each member fixture.

- [ ] **Step 3: Run the affected Vitest files**

Run: `npm test -- MemberMetaRow hoverPhaseColoring`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/components/MemberMetaRow.tsx test/MemberMetaRow.test.tsx test/MemberMetaRow.collapse.test.tsx test/MemberMetaRow.drag.test.tsx test/hoverPhaseColoring.test.tsx
git commit -m "Migrate MemberMetaRow to SeriesMember (I3.2)"
```

---

## Task 7: Migrate `MemberMetaGutter.tsx`

**Files:**
- Modify: `src/components/MemberMetaGutter.tsx` (`ComparisonMember` at lines 22, 31)
- Test: `test/MemberMetaGutter.reorder.test.tsx`, `test/MemberMetaGutter.resize.test.tsx`, `test/MemberReorder.test.tsx`

- [ ] **Step 1: Rename in `MemberMetaGutter.tsx`**

Replace every occurrence of the identifier `ComparisonMember` with `SeriesMember` (2 occurrences). Line 22 becomes:

```ts
import type { SeriesMember } from "../api";
```

- [ ] **Step 2: Retype the test fixtures**

In `test/MemberMetaGutter.reorder.test.tsx`, `test/MemberMetaGutter.resize.test.tsx`, and `test/MemberReorder.test.tsx`, replace every occurrence of `ComparisonMember` with `SeriesMember`, and replace `comparison_id` with `series_id` in each member fixture.

- [ ] **Step 3: Run the affected Vitest files**

Run: `npm test -- MemberMetaGutter MemberReorder`
Expected: PASS.

- [ ] **Step 4: Commit**

```bash
git add src/components/MemberMetaGutter.tsx test/MemberMetaGutter.reorder.test.tsx test/MemberMetaGutter.resize.test.tsx test/MemberReorder.test.tsx
git commit -m "Migrate MemberMetaGutter to SeriesMember (I3.2)"
```

---

## Task 8: Migrate `draftFactories.ts` (with the comparison-factory bridge)

**Files:**
- Modify: `src/lib/comparison/draftFactories.ts` (import line 16; `memberFromSaved` signature line 21; `memberFromSaved` call sites at line 77 and line 107)
- Test: `test/lib/comparison/draft.test.ts` (exercises these factories)

`memberFromSaved` is the reusable, render-adjacent helper — it migrates to `SeriesMember`. `fromComparison` / `fromComparisonAsFork` are comparison-only factories (parameter type `Comparison`) that I3.6 deletes; they keep taking `Comparison` and bridge the now-`SeriesMember`-typed `memberFromSaved` call. `memberFromNewExposure` has no member-type parameter — leave it alone.

- [ ] **Step 1: Update the import**

Change line 16 from:

```ts
import type { Comparison, ComparisonMember, MemberSnapshot } from "../../api";
```

to:

```ts
import type { Comparison, MemberSnapshot, SeriesMember } from "../../api";
```

- [ ] **Step 2: Migrate the `memberFromSaved` signature**

Change line 21 from `  saved: ComparisonMember,` to:

```ts
  saved: SeriesMember,
```

- [ ] **Step 3: Bridge the `fromComparison` call site**

`Comparison.members` is `ComparisonMember[]`, which lacks `series_id`. Add the field via a spread so the object satisfies `SeriesMember` (a spread-introduced `comparison_id` is exempt from excess-property checking; the value of `series_id` is irrelevant — `memberFromSaved` never reads it). Change line 77 from:

```ts
    members: c.members.map((m) => memberFromSaved(m, qc)),
```

to:

```ts
    // I3.2 bridge — comparison members lack series_id; the render pipeline
    // never reads it. Deleted with this factory at I3.6.
    members: c.members.map((m) => memberFromSaved({ ...m, series_id: 0 }, qc)),
```

- [ ] **Step 4: Bridge the `fromComparisonAsFork` call site**

Change line 107 from `      const dm = memberFromSaved(m, qc);` to:

```ts
      const dm = memberFromSaved({ ...m, series_id: 0 }, qc); // I3.2 bridge — see fromComparison
```

- [ ] **Step 5: Run the affected Vitest file**

Run: `npm test -- draft`
Expected: PASS — `draftFactories` behaviour (snapshot recompute, id handling) is unchanged.

- [ ] **Step 6: Commit**

```bash
git add src/lib/comparison/draftFactories.ts
git commit -m "Migrate memberFromSaved to SeriesMember with the I3.6-scoped bridge (I3.2)"
```

---

## Task 9: Bridge `Compare.tsx`

**Files:**
- Modify: `src/pages/Compare.tsx` (import lines 75-77; `ComparisonMember` at lines 399, 571, 978, 1050; `comparison_id` at line 574)

`Compare.tsx` is deleted at I3.6 (#177). Until then, it feeds the retyped render pipeline. The only render-bound member array is `plotMembers`, built fresh by the local `draftToMember` converter — `draftToMember` already fabricates a placeholder foreign key (`comparison_id: 0`), so this is a pure rename, no cast needed.

- [ ] **Step 1: Update the api import block**

Replace lines 75-77:

```ts
import type {
  Comparison, ComparisonMember, ComparisonMemberInput, SaveComparisonBody,
} from "../api";
```

with:

```ts
import type {
  Comparison, ComparisonMemberInput, SaveComparisonBody, SeriesMember,
} from "../api";
```

- [ ] **Step 2: Retype the `draftToMember` converter**

On line 571, change the return type `function draftToMember(d: DraftMember): ComparisonMember {` to:

```ts
function draftToMember(d: DraftMember): SeriesMember {
```

On line 574, change the placeholder field `    comparison_id: 0,` to:

```ts
    series_id: 0,
```

- [ ] **Step 3: Retype the two `sampleIdFor` callbacks and the `plotMembers` memo**

- Line 399: change `(m: ComparisonMember): number | null =>` to `(m: SeriesMember): number | null =>`.
- Line 1050: change `(m: ComparisonMember): number | null =>` to `(m: SeriesMember): number | null =>`.
- Line 978: change `useMemo<ComparisonMember[]>` to `useMemo<SeriesMember[]>`.

- [ ] **Step 4: Run the Compare Vitest suite**

Run: `npm test -- Compare`
Expected: PASS — all `Compare.*` specs behave identically.

- [ ] **Step 5: Commit**

```bash
git add src/pages/Compare.tsx
git commit -m "Bridge Compare.tsx onto the SeriesMember render pipeline (I3.2)"
```

---

## Task 10: Full typecheck + test gate

**Files:**
- (verification only)

- [ ] **Step 1: Run the full typecheck**

Run: `npx tsc --noEmit`
Expected: PASS — every consumer is now consistent on `SeriesMember`.

If any error remains, it points at a file with a stray `ComparisonMember` fixture passed to a migrated component (or a `comparison_id` field still in such a fixture). Apply the same rename — `ComparisonMember` → `SeriesMember`, `comparison_id` → `series_id` — to that file and re-run. Do **not** touch `test/draftPersistence.test.ts`, `test/MentionChip.test.tsx`, or `test/queue/saveComparison.test.tsx`: those legitimately keep `ComparisonMember` (they do not feed a migrated component) and should not appear in the error list.

- [ ] **Step 2: Run the full unit-test suite**

Run: `npm test`
Expected: PASS — no behavioural change anywhere.

- [ ] **Step 3: Run the production build**

Run: `npm run build`
Expected: PASS (`tsc --noEmit` + `vite build` both green).

- [ ] **Step 4: Confirm `yBands.ts` is untouched**

Run: `git diff --name-only main -- src/lib/comparison/yBands.ts`
Expected: empty output — `yBands.ts` is a pure numeric module excluded from this migration.

- [ ] **Step 5: Commit (if Step 1 required a residual fix)**

```bash
git add -A
git commit -m "Resolve residual SeriesMember typecheck errors (I3.2)"
```

If Step 1 passed with no residual fix, skip this commit.

---

## Acceptance criteria (from the spec)

- [ ] `SeriesMember` is defined as a standalone interface in `api.ts`, mirroring the I2.2 route response; optional fields modelled `T | null`.
- [ ] The nine listed modules are retyped onto `SeriesMember`.
- [ ] `lib/comparison/yBands.ts` is unchanged.
- [ ] `tsc --noEmit` passes.
- [ ] The render pipeline's existing Vitest passes unchanged.

## Coordination note

If the series event-kind cluster (#166/#167/#168) has merged a `Series` parent type into `api.ts` before this lands, reconcile so there is a single `SeriesMember` (this issue owns it) and `Series.members` references it. The additions are textually separate interface blocks — low conflict — but check at rebase. See the spec's "Coordination notes" section.
