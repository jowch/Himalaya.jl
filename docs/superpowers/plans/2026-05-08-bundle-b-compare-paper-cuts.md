# Bundle B — Compare paper cuts (#73, #74, #75, #77, #78) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Land five small post-sprint Compare-page bug fixes — drift-seam cleanup, conditional hint, popover dismissal affordances, conflict-modal cold-cache prefetch, and URL/Zustand sync gap.

**Architecture:** All five issues are independent surface-level fixes localized to one or two files each. Tasks are TDD: failing test → minimal implementation → green → commit. They can be worked sequentially or dispatched in parallel via subagents (no shared state between tasks).

**Tech Stack:** React 18 + TypeScript strict, Vitest + React Testing Library, Zustand (client state), TanStack Query (server state), react-router-dom v6.

---

## File map

Per-task touched files. Each task is its own commit.

| Task | Issue | Implementation file(s) | Test file |
|------|-------|------------------------|-----------|
| 1 | #73 | `src/components/MemberMetaRow.tsx` | `test/MemberMetaRow.test.tsx` |
| 2 | #78 | `src/pages/ComparePageEdit.tsx` | `test/ComparePageEdit.test.tsx` |
| 3 | #75 | `src/components/ForksPopover.tsx` | `test/ForksPopover.test.tsx` |
| 4 | #74 | `src/components/ConflictModal.tsx` | `test/ConflictModal.test.tsx` |
| 5 | #77 | `src/components/AppShell.tsx` | `test/AppShell.test.tsx` (new) |

All paths are relative to `packages/HimalayaUI/frontend/`.

---

## Task 1: #73 — Drop `MemberMetaRow.defaultLabel`

**Files:**
- Modify: `src/components/MemberMetaRow.tsx` (lines 51, 68–72, 235–237)
- Test: `test/MemberMetaRow.test.tsx`

**Why:** `MemberMetaRow.tsx` defines a local `defaultLabel(member)` (lines 68–72) used as a fallback when `displayLabel` is not passed. The canonical chain lives in `lib/comparison/labels.ts::resolveDisplayLabels` (extended in PR #70 to `label_override → "${sample.name} · ${exposure.filename}" → exposure.filename → sample.name → "Exposure #N" → "(deleted exposure)"`). The local fallback is **truncated** (skips three rungs) and currently unreachable because both `ComparePage.tsx` and `ComparePageEdit.tsx` always pass `displayLabel`. Eliminating the dead code closes the drift seam — same class as #62/#69. Make `displayLabel` required so future callers can't silently regress.

- [ ] **Step 1: Write a failing test that asserts `displayLabel` is the only label source**

Add this test to `test/MemberMetaRow.test.tsx`. Place it inside the existing `describe("MemberMetaRow — review mode", ...)` block (after the existing tests, before the closing brace).

```tsx
it("renders the displayLabel prop verbatim (no internal fallback chain)", () => {
  render(
    <MemberMetaRow
      member={makeMember()}
      top={0}
      height={50}
      mode="review"
      displayLabel="JC068P · run-007.dat"
    />,
  );
  expect(screen.getByTestId("member-meta-label").textContent).toBe(
    "JC068P · run-007.dat",
  );
});
```

- [ ] **Step 2: Run test to verify it passes today (it should — current code uses `displayLabel ?? defaultLabel(member)`)**

Run from `packages/HimalayaUI/frontend/`:
```
node_modules/.bin/vitest run test/MemberMetaRow.test.tsx -t "displayLabel prop verbatim"
```
Expected: PASS. (This test is a regression pin for after we delete `defaultLabel`.)

- [ ] **Step 3: Make `displayLabel` required and delete `defaultLabel`**

In `src/components/MemberMetaRow.tsx`:

1. Change the `displayLabel?: string` prop declaration (line 51) to required:
   ```ts
     /** Pre-resolved display label from `resolveDisplayLabels` (lib/comparison/labels.ts). */
     displayLabel: string;
   ```
   Drop the lengthy inline doc; the one-liner pointing to the canonical resolver is enough.

2. Delete the `defaultLabel` function (lines 68–72 — the entire `function defaultLabel(member: ComparisonMember): string { ... }` block).

3. At the two `displayLabel ?? defaultLabel(member)` call sites (lines 235 and 237), drop the fallback so they read just `displayLabel`:
   ```tsx
   <span
     data-testid="member-meta-label"
     className="truncate min-w-0 max-w-[16ch] text-fg"
     title={displayLabel}
   >
     {displayLabel}
   </span>
   ```

- [ ] **Step 4: Update existing tests in `test/MemberMetaRow.test.tsx` that don't pass `displayLabel`**

Every `<MemberMetaRow ...>` render call in the test file currently omits `displayLabel`. Add `displayLabel="…"` to every render call. The simplest mechanical fix is to pass a stable string like `"row-label"` on tests that don't care about the label content; tests that already pass `displayLabel` are unchanged.

Do this with a single search across the file: any line with `<MemberMetaRow` followed by `member={…}` lacking `displayLabel=` gets `displayLabel="row-label"` added before the closing `/>`.

- [ ] **Step 5: Run the full test file to verify everything passes**

```
node_modules/.bin/vitest run test/MemberMetaRow.test.tsx
```
Expected: PASS for all tests.

- [ ] **Step 6: Run `tsc --noEmit` to verify no callers broke**

```
npm run build
```
Expected: clean build. If `ComparePage.tsx`, `ComparePageEdit.tsx`, or any other file fails because they don't pass `displayLabel`, fix the call site (they already construct `displayLabel` via `resolveDisplayLabels` — the issue body confirms this).

- [ ] **Step 7: Commit**

```
git add packages/HimalayaUI/frontend/src/components/MemberMetaRow.tsx \
        packages/HimalayaUI/frontend/test/MemberMetaRow.test.tsx
git commit -m "refactor(compare): drop MemberMetaRow.defaultLabel drift seam (#73)

Make displayLabel required; delete the truncated local fallback.
The canonical chain lives in lib/comparison/labels.ts::resolveDisplayLabels.
Closes the seam from #62/#69."
```

---

## Task 2: #78 — Edit hint conditional

**Files:**
- Modify: `src/pages/ComparePageEdit.tsx` (lines 605–612)
- Test: `test/ComparePageEdit.test.tsx`

**Why:** The right-slot HintText "Use '+ Add traces' to populate the comparison." is rendered unconditionally in edit mode. It's correct for empty drafts (`/compare/.../new`) but misleading when editing a populated comparison. Gate on `draft.members.length === 0`.

- [ ] **Step 1: Inspect the existing test file structure**

Open `test/ComparePageEdit.test.tsx` and look at how it sets up Zustand `activeDraft`. Find the existing render helper or describe block that boots a populated draft. The test you'll add follows the same pattern.

- [ ] **Step 2: Write a failing test that asserts the hint disappears when members.length > 0**

Add this test to `test/ComparePageEdit.test.tsx` inside the most appropriate `describe(...)` block (likely the one covering edit-mode rendering, near where other right-slot or empty-state tests live):

```tsx
it("right-slot hint hides when the draft has members (#78)", async () => {
  // Arrange: render the edit page with a draft that has at least one member.
  // Use whatever helper the file already has for booting edit mode with a
  // populated draft (e.g. the same setup the "Save button is enabled" test uses).
  await renderEditWithMembers({ memberCount: 2 });

  // The hint testid is "compare-edit-right-hint". When draft has members,
  // the hint must not render at all.
  expect(screen.queryByTestId("compare-edit-right-hint")).toBeNull();
});

it("right-slot hint shows when the draft is empty (#78)", async () => {
  await renderEditWithMembers({ memberCount: 0 });
  expect(screen.getByTestId("compare-edit-right-hint")).toBeInTheDocument();
});
```

If `renderEditWithMembers` doesn't exist in the file, use the existing render helper (whatever it's called — `renderEdit`, `renderPage`, etc.) and seed the draft via `useAppState.setState({ activeDraft: { ..., members: [...] } })` before the render call. Look at neighbouring tests for the exact shape; do not invent a new helper.

- [ ] **Step 3: Run the new tests to verify the empty-case passes and the populated-case fails**

```
node_modules/.bin/vitest run test/ComparePageEdit.test.tsx -t "right-slot hint"
```
Expected: the empty-case test passes, the populated-case test FAILS (hint is rendered unconditionally today).

- [ ] **Step 4: Make the right slot conditional on members.length**

In `src/pages/ComparePageEdit.tsx` lines 605–612, change the `right={...}` prop on the `WorkspaceGrid` to gate on the draft's member count. The component already has `draft` available in scope as the result of `useAppState((s) => s.activeDraft)` — confirm by re-reading the surrounding code if uncertain.

Replace:
```tsx
right={
  <div
    data-testid="compare-edit-right-hint"
    className="h-full flex items-center justify-center p-4 text-center"
  >
    <HintText>Use “+ Add traces” to populate the comparison.</HintText>
  </div>
}
```

With:
```tsx
right={
  (draft?.members.length ?? 0) === 0 ? (
    <div
      data-testid="compare-edit-right-hint"
      className="h-full flex items-center justify-center p-4 text-center"
    >
      <HintText>Use “+ Add traces” to populate the comparison.</HintText>
    </div>
  ) : null
}
```

- [ ] **Step 5: Run the tests to verify both cases pass**

```
node_modules/.bin/vitest run test/ComparePageEdit.test.tsx -t "right-slot hint"
```
Expected: PASS for both.

- [ ] **Step 6: Run `tsc --noEmit` to confirm no type regressions**

```
npm run build
```
Expected: clean.

- [ ] **Step 7: Commit**

```
git add packages/HimalayaUI/frontend/src/pages/ComparePageEdit.tsx \
        packages/HimalayaUI/frontend/test/ComparePageEdit.test.tsx
git commit -m "fix(compare): hide edit-mode right-slot hint when draft is populated (#78)

The 'Use + Add traces to populate' copy was rendered unconditionally,
showing on populated comparisons in edit mode and confusing users.
Gate on draft.members.length === 0."
```

---

## Task 3: #75 — `ForksPopover` Esc / outside-click / aria

**Files:**
- Modify: `src/components/ForksPopover.tsx`
- Test: `test/ForksPopover.test.tsx`

**Why:** `ForksPopover.tsx` currently has only a click-to-toggle trigger. No Esc, no outside-click dismiss, no `aria-expanded`. On fork-heavy comparisons the popover panel partially occludes the trigger so click-again-to-close is invisible. Add all three affordances.

- [ ] **Step 1: Write three failing tests in `test/ForksPopover.test.tsx`**

Add inside the existing `describe("ForksPopover", ...)` block:

```tsx
it("Esc closes the popover (#75)", () => {
  renderPopover({ parentId: 100, experimentId: 7, forks: [] });
  fireEvent.click(screen.getByTestId("comparison-forks-trigger"));
  expect(screen.getByTestId("comparison-forks-popover")).toBeInTheDocument();

  fireEvent.keyDown(document, { key: "Escape" });
  expect(screen.queryByTestId("comparison-forks-popover")).toBeNull();
});

it("outside click closes the popover (#75)", () => {
  const { container } = renderPopover({ parentId: 100, experimentId: 7, forks: [] });
  fireEvent.click(screen.getByTestId("comparison-forks-trigger"));
  expect(screen.getByTestId("comparison-forks-popover")).toBeInTheDocument();

  // Click an element outside the popover and outside the trigger.
  fireEvent.mouseDown(container);
  expect(screen.queryByTestId("comparison-forks-popover")).toBeNull();
});

it("trigger reflects open state via aria-expanded (#75)", () => {
  renderPopover({ parentId: 100, experimentId: 7, forks: [] });
  const trigger = screen.getByTestId("comparison-forks-trigger");
  expect(trigger).toHaveAttribute("aria-expanded", "false");
  expect(trigger).toHaveAttribute("aria-haspopup", "menu");

  fireEvent.click(trigger);
  expect(trigger).toHaveAttribute("aria-expanded", "true");
});
```

- [ ] **Step 2: Run tests to verify all three FAIL**

```
node_modules/.bin/vitest run test/ForksPopover.test.tsx
```
Expected: 3 FAILs (Esc has no handler; outside click does nothing; aria attributes missing).

- [ ] **Step 3: Implement Esc + outside-click + aria in `src/components/ForksPopover.tsx`**

Add `useEffect, useRef` to the existing `useState` import:
```ts
import { useEffect, useRef, useState } from "react";
```

Inside the `ForksPopover` function, add refs for the trigger and panel:
```ts
const triggerRef = useRef<HTMLButtonElement>(null);
const panelRef = useRef<HTMLDivElement>(null);
```

Add the two effects right after the `forks` const declaration (line 43):
```ts
useEffect(() => {
  if (!open) return;
  const onKey = (e: KeyboardEvent): void => {
    if (e.key === "Escape") setOpen(false);
  };
  document.addEventListener("keydown", onKey);
  return () => document.removeEventListener("keydown", onKey);
}, [open]);

useEffect(() => {
  if (!open) return;
  const onDown = (e: MouseEvent): void => {
    const target = e.target as Node | null;
    if (target === null) return;
    if (panelRef.current?.contains(target)) return;
    if (triggerRef.current?.contains(target)) return;
    setOpen(false);
  };
  document.addEventListener("mousedown", onDown);
  return () => document.removeEventListener("mousedown", onDown);
}, [open]);
```

Wire the refs and aria attributes into the markup. Update the `<button>` (lines 47–55):
```tsx
<button
  ref={triggerRef}
  type="button"
  data-testid="comparison-forks-trigger"
  aria-expanded={open}
  aria-haspopup="menu"
  onClick={() => setOpen((v) => !v)}
  className="px-2 py-0.5 rounded border border-border text-fg-muted text-xs
             hover:bg-bg-elevated"
>
  Forks ({forks.length}) →
</button>
```

Update the popover panel `<div>` (line 57) to attach the ref:
```tsx
<div
  ref={panelRef}
  data-testid="comparison-forks-popover"
  ...
>
```

- [ ] **Step 4: Run all `ForksPopover` tests**

```
node_modules/.bin/vitest run test/ForksPopover.test.tsx
```
Expected: all PASS, including the new three.

- [ ] **Step 5: Run `tsc --noEmit`**

```
npm run build
```
Expected: clean.

- [ ] **Step 6: Commit**

```
git add packages/HimalayaUI/frontend/src/components/ForksPopover.tsx \
        packages/HimalayaUI/frontend/test/ForksPopover.test.tsx
git commit -m "fix(compare): ForksPopover Esc + outside-click + aria-expanded (#75)

Closing affordances were missing — Esc did nothing, clicking outside the
panel did nothing, and trigger had no aria-expanded. On fork-heavy
comparisons the panel occluded the trigger so the click-again-to-close
target was invisible. Mirror the NavModal dismissal pattern."
```

---

## Task 4: #74 — `ConflictModal` cold-cache prefetch

**Files:**
- Modify: `src/components/ConflictModal.tsx` (lines 36–46 imports, 75–114 `buildOverwritePayload`, 192–197 `handleOverwrite`)
- Test: `test/ConflictModal.test.tsx`

**Why:** `buildOverwritePayload` (lines 75–113) is symmetric to `ComparePageEdit::handleSave` (lines 189–262), but it reads `analysis_inputs_hash` from the TanStack cache without prefetching cold members. `handleSave` was extended in #49 to prefetch four cache keys (`exposure`, `peaks`, `indices`, `groups`) when any is missing; `buildOverwritePayload` did not get the same fix. Today the bug is masked because Overwrite is reached via Save → 409 (so the prefetch already ran), but the regression mode is real (long-idle conflict modal, future codepaths, test fixtures).

- [ ] **Step 1: Write a failing test in `test/ConflictModal.test.tsx`**

Add a new test that:
1. Mounts the modal with `pendingConflict` set (so the modal is open).
2. Loads a draft whose member references an exposure id whose `queryKeys.exposure(id)` / `peaks(id)` / `indices(id)` / `groups(id)` are NOT pre-seeded in the QueryClient (cold cache).
3. Spies on the `getExposure` / `listPeaks` / `listIndices` / `listGroups` API helpers (or asserts via fetch-mock) that they are called when the user clicks Overwrite.
4. Assert: clicking `[data-testid="conflict-overwrite"]` calls all four prefetch helpers for the cold member's exposure id BEFORE `save.mutate` is called.

The simplest test shape (mirroring how `ComparePageEdit.test.tsx` does this — go read it first if you're writing this fresh):

```tsx
it("Overwrite prefetches cold member caches before computing snapshots (#74)", async () => {
  const fetchSpy = vi.spyOn(global, "fetch");
  // Helpers in the test file should already exist for seeding draft + opening
  // the modal. If not, follow the pattern in this file's existing
  // "Overwrite re-submits" test and adapt: seed a draft with a member whose
  // exposure_id is, say, 999, and DO NOT call qc.setQueryData(...) for any of
  // the four queryKeys for id=999.

  await renderModalWithConflict({
    draftMemberExposureIds: [999],   // cold
    serverHash: "sha256:newer",
  });

  const overwriteBtn = screen.getByTestId("conflict-overwrite");
  await userEvent.click(overwriteBtn);

  // Assert all four endpoints were hit for the cold exposure.
  await waitFor(() => {
    const urls = fetchSpy.mock.calls.map((c) =>
      typeof c[0] === "string" ? c[0] : String(c[0]),
    );
    expect(urls.some((u) => u.includes("/api/exposures/999"))).toBe(true);
    expect(urls.some((u) => u.includes("/api/exposures/999/peaks"))).toBe(true);
    expect(urls.some((u) => u.includes("/api/exposures/999/indices"))).toBe(true);
    expect(urls.some((u) => u.includes("/api/exposures/999/groups"))).toBe(true);
  });
});
```

Naming details (helpers, exact api fetch URL shapes) may differ — read existing tests in this file to use the same fixtures and helpers.

- [ ] **Step 2: Run the test to verify it FAILS today**

```
node_modules/.bin/vitest run test/ConflictModal.test.tsx -t "prefetches cold member caches"
```
Expected: FAIL — without prefetch, `fetch` is not called for the cold exposure id at all when Overwrite fires (the snapshot is computed against an empty cache, then the mutation runs).

- [ ] **Step 3: Implement the prefetch in `buildOverwritePayload`**

In `src/components/ConflictModal.tsx`:

1. Extend the imports to pull in the four api helpers + `queryKeys`:
   ```ts
   import { computeMemberSnapshot } from "../lib/comparison/snapshot";
   import { comparePath, type CompareScope } from "../lib/comparison/routes";
   import { getExposure, listPeaks, listIndices, listGroups } from "../api";
   import { queryKeys } from "../queries";
   ```
   (Verify the exact path of `queryKeys` by re-reading `ComparePageEdit.tsx` line 44 — the import must match.)

2. Convert `buildOverwritePayload` from sync to `async`:

   Replace the existing function body (lines 75–114) with the version below. The change set is: prefetch any cold exposures (mirroring `ComparePageEdit.tsx::handleSave` lines 189–221) before computing snapshots; everything else is identical to today's body.

   ```ts
   async function buildOverwritePayload(
     draft: ActiveDraft,
     serverHash: string,
     qc: ReturnType<typeof useQueryClient>,
   ): Promise<SaveComparisonBody & { id?: number }> {
     // Mirror handleSave's cold-cache prefetch (#49) — without it, snapshots
     // for never-visited members land with analysis_inputs_hash = "" and the
     // server marks them stale on the next view fold. See issue #74.
     const coldExposureIds = draft.members
       .map((m) => m.exposure_id)
       .filter((id): id is number => id !== null)
       .filter((id) =>
         qc.getQueryData(queryKeys.exposure(id)) === undefined
         || qc.getQueryData(queryKeys.peaks(id)) === undefined
         || qc.getQueryData(queryKeys.indices(id)) === undefined
         || qc.getQueryData(queryKeys.groups(id)) === undefined
       );
     if (coldExposureIds.length > 0) {
       await Promise.all(
         coldExposureIds.flatMap((id) => [
           qc.fetchQuery({
             queryKey: queryKeys.exposure(id),
             queryFn: () => getExposure(id),
           }),
           qc.fetchQuery({
             queryKey: queryKeys.peaks(id),
             queryFn: () => listPeaks(id),
           }),
           qc.fetchQuery({
             queryKey: queryKeys.indices(id),
             queryFn: () => listIndices(id),
           }),
           qc.fetchQuery({
             queryKey: queryKeys.groups(id),
             queryFn: () => listGroups(id),
           }),
         ]),
       );
     }

     const members: ComparisonMemberInput[] = draft.members.map((m) => {
       const snapshot = m.exposure_id !== null
         ? computeMemberSnapshot(m.exposure_id, qc)
         : (m.snapshot ?? {
             effective_peaks: [],
             confirmed_index: null,
             analysis_inputs_hash: "",
           });
       const out: ComparisonMemberInput = {
         exposure_id: m.exposure_id,
         display_order: m.display_order,
         band_height: m.band_height,
         y_offset: m.y_offset,
         normalization: m.normalization,
         snapshot,
       };
       if (m.id !== undefined) out.id = m.id;
       if (m.color_override !== undefined) out.color_override = m.color_override;
       if (m.label_override !== undefined) out.label_override = m.label_override;
       if (m.q_window_min !== undefined) out.q_window_min = m.q_window_min;
       if (m.q_window_max !== undefined) out.q_window_max = m.q_window_max;
       if (m.peak_display !== undefined) out.peak_display = m.peak_display;
       return out;
     });
     const payload: SaveComparisonBody & { id?: number } = {
       title: draft.title,
       members,
       expected_content_hash: serverHash,
     };
     if (draft.id !== undefined) payload.id = draft.id;
     if (draft.description !== "") payload.description = draft.description;
     if (draft.forkedFromId !== undefined) payload.forked_from_id = draft.forkedFromId;
     if (draft.forkedAtHash !== undefined) payload.forked_at_hash = draft.forkedAtHash;
     return payload;
   }
   ```

3. Update `handleOverwrite` to await the now-async builder. Replace lines 192–197:

   ```ts
   const handleOverwrite = useCallback(async () => {
     if (draft === null || serverHash === null) return;
     if (save.isPending || overwriteInFlightRef.current) return;
     overwriteInFlightRef.current = true;
     try {
       const payload = await buildOverwritePayload(draft, serverHash, qc);
       save.mutate(payload);
     } catch {
       // Prefetch failed — release the in-flight guard so the user can retry.
       overwriteInFlightRef.current = false;
     }
   }, [draft, serverHash, save, qc]);
   ```

   The `try/catch` mirrors `handleSave`'s guard-release on prefetch failure (`ComparePageEdit.tsx` lines 263–266).

- [ ] **Step 4: Run the new test**

```
node_modules/.bin/vitest run test/ConflictModal.test.tsx -t "prefetches cold member caches"
```
Expected: PASS.

- [ ] **Step 5: Run the full ConflictModal test suite to confirm no regressions**

```
node_modules/.bin/vitest run test/ConflictModal.test.tsx
```
Expected: all PASS.

- [ ] **Step 6: Run `tsc --noEmit`**

```
npm run build
```
Expected: clean.

- [ ] **Step 7: Commit**

```
git add packages/HimalayaUI/frontend/src/components/ConflictModal.tsx \
        packages/HimalayaUI/frontend/test/ConflictModal.test.tsx
git commit -m "fix(compare): ConflictModal Overwrite prefetches cold caches (#74)

Mirror handleSave's #49 fix — without prefetch, members whose exposure
caches are cold land analysis_inputs_hash='' in the overwrite payload,
silently regressing #49. Today this is masked because Overwrite is reached
via Save→409 (caches already warm); the fix removes the regression mode."
```

---

## Task 5: #77 — Compare cold-load empty body

**Files:**
- Modify: `src/components/AppShell.tsx` (lines 51–58 — extend the URL/Zustand sync)
- Test: `test/AppShell.test.tsx` (new file)

**Why:** When the user reloads at `/` with persisted `activePage: "compare"` in localStorage, the rocker highlights "Compare" but `ZustandShellPage` (lines 27–35) renders nothing because it only branches on `index` and `inspect`. There's a URL→Zustand sync (`onComparePath` → `setActivePage("compare")`) but no symmetric Zustand→URL sync. Add the converse: when `activePage === "compare"` AND we're not on a compare path, navigate to `/experiments/:eid/compare` (if `experimentId` is set) or `/compare/all`.

- [ ] **Step 1: Create `test/AppShell.test.tsx` with a failing test**

```tsx
/**
 * AppShell — URL/Zustand sync regression tests (#77).
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppShell } from "../src/components/AppShell";

function makeQc() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function renderShell(initialPath: string) {
  const qc = makeQc();
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[initialPath]}>
        <AppShell />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppShell — Zustand → URL sync (#77)", () => {
  beforeEach(() => {
    // Reset Zustand back to defaults between cases.
    useAppState.setState({
      activePage: "index",
      activeExperimentId: undefined,
    });
  });

  it("activePage='compare' + URL '/' navigates to /compare/all", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: undefined });
    renderShell("/");

    // Compare page should mount; the empty rocker-only render is the bug.
    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });

  it("activePage='compare' + URL '/' navigates to /experiments/:eid/compare when experiment is set", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: 7 });
    renderShell("/");

    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
    // Optional: assert the URL itself, e.g. via a useLocation probe component.
    // For the regression test, "ComparePage rendered" is the load-bearing
    // assertion — the empty body is the bug.
  });
});
```

Note: `data-testid="compare-page"` should be checked against the actual ComparePage root. If the testid is different (likely something like `compare-page-list` or similar), update the assertion. Run a quick `grep -rn "data-testid=\"compare-page" packages/HimalayaUI/frontend/src/pages/ComparePage.tsx` first to find the canonical testid.

- [ ] **Step 2: Run tests to verify they FAIL**

```
node_modules/.bin/vitest run test/AppShell.test.tsx
```
Expected: both FAIL — empty body is rendered today.

- [ ] **Step 3: Add the symmetric sync effect in `src/components/AppShell.tsx`**

Right after the existing `useEffect` at lines 56–58 (the URL→Zustand sync), add the converse Zustand→URL sync. Get `activePage` from Zustand (the existing code only sets it; we need to read it):

Insert this block right after line 58:

```ts
const activePage = useAppState((s) => s.activePage);

// Symmetric: when activePage is "compare" but URL isn't on a compare path,
// navigate so the URL-routed Compare pages mount. Without this, a reload
// at "/" with persisted activePage='compare' renders the rocker but no
// page body (issue #77).
useEffect(() => {
  if (activePage !== "compare") return;
  if (onComparePath) return;
  const url = experimentId !== undefined
    ? `/experiments/${experimentId}/compare`
    : "/compare/all";
  navigate(url, { replace: true });
}, [activePage, onComparePath, experimentId, navigate]);
```

- [ ] **Step 4: Run the new tests**

```
node_modules/.bin/vitest run test/AppShell.test.tsx
```
Expected: both PASS.

- [ ] **Step 5: Run any neighbouring AppShell-touching tests to confirm no regression**

```
node_modules/.bin/vitest run test/TabRocker.test.tsx test/state.test.ts
```
Expected: PASS.

- [ ] **Step 6: Run `tsc --noEmit`**

```
npm run build
```
Expected: clean.

- [ ] **Step 7: Commit**

```
git add packages/HimalayaUI/frontend/src/components/AppShell.tsx \
        packages/HimalayaUI/frontend/test/AppShell.test.tsx
git commit -m "fix(compare): symmetric activePage→URL sync (#77)

When persisted activePage='compare' but URL is '/', the existing URL→Zustand
sync did nothing (URL was never on a compare path) and ZustandShellPage
rendered empty body (it only branches on index/inspect). Add the converse
Zustand→URL effect so the Compare routes mount on cold load."
```

---

## Final verification

After all five tasks are committed, run a full check:

- [ ] **Run the full Vitest suite**

```
cd packages/HimalayaUI/frontend
npm test > /tmp/vitest.out 2>&1
grep -E "FAIL|fail|✗" /tmp/vitest.out | head -20
tail -20 /tmp/vitest.out
```
Expected: 0 fails. Per CLAUDE.md, capture once and grep the file rather than re-running.

- [ ] **Run `npm run build`**

```
npm run build
```
Expected: `tsc --noEmit` clean + Vite build succeeds.

- [ ] **Optional: open the dev server and smoke-test each fix**

```
npm run dev
```
Then in the browser:
- Open a populated comparison → click Edit → confirm right slot is empty (no hint copy).
- Open a comparison with forks → click Forks → press Esc → popover closes; reopen → click outside → closes; check accessibility tree for `aria-expanded`.
- Reload at `/` with localStorage `himalaya-ui:state` having `activePage: "compare"` → page body mounts (not just rocker).

(#73 and #74 are not directly user-visible without specific scenarios — covered by the unit tests.)

---

## Self-review checklist

- **Spec coverage:** All five issues map to one task each (#73→T1, #78→T2, #75→T3, #74→T4, #77→T5). ✓
- **No placeholders:** Each step has exact code or exact commands. ✓
- **Type consistency:** `displayLabel` made required in T1 → callers in `ComparePage.tsx` and `ComparePageEdit.tsx` already pass it (per issue #73 body, the prop is always passed today). T4 makes `buildOverwritePayload` async, and the only call site (`handleOverwrite`) is updated to `await`. T5's `useEffect` is added in `AppShell` where `experimentId` and `navigate` are already in scope (lines 39, 42). ✓
- **TDD discipline:** Each task starts with a failing test, then minimal implementation, then green. ✓
- **Frequent commits:** One commit per task; no batching. ✓
