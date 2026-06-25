# Ingestion Redesign — Phase E2: Grouping-Review Surface + Configuration Components Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the interactive grouping-review surface (the centerpiece — Load ▸ Sample ▸ Exposures fold tree with rename / split / merge / move, fuzzy/glob filter with selection that persists across filtering, bulk-merge, and free undo) and the Configuration-tab internals (geometry ledger with per-field provenance + override + revert, an editable Sources list, and an acquisition plot built on the existing trace-plot infrastructure). Each structural edit is a queue mutation: it mints a `client_op_id`, writes an optimistic cache effect, and reconciles against the Phase-D backend event kinds over SSE.

> **Backend event-kind reality (spec §9.3, resolved 2026-06-19) — load-bearing for every cache arm and mutator below:** there is **no `sample_merged` event**. A *merge* fans out as a burst of per-exposure `exposure_moved` frames (+ a `sample_renamed`/`update_sample` if the survivor is renamed) and retires the loser via `merged_into_id`; the frontend reconciles a foreign merge from that `exposure_moved` burst — it has **no `sample_merged` arm**. A *split* emits `sample_created` (entity = the new sample) plus `sample_split` (entity = the source sample). The full kind set the frontend handles is therefore: `exposure_moved`, `sample_renamed`, `sample_created`, `sample_split`, and `grouping_flag_dismissed`. Canonical payloads (all carry `experiment_id`; the frontend reads payload fields, never `remote.entity_id`, for scope): `exposure_moved {sample_id, from_sample_id, experiment_id}`; `sample_renamed {name, experiment_id}`; `sample_created {experiment_id}`; `sample_split {new_sample_id, exposure_ids, experiment_id}`; `grouping_flag_dismissed {flag_kind, merge_with_sample_id?, experiment_id}`.

This phase owns the optimistic **send** path (the five mutators + their `onSuccess` `loads(experimentId)` own-op reconcile, the structural `OpKind`s, `resolveMutatorForEvent` wiring, and the mutation hooks). The SSE **receive** path — the `applyRemoteToCache` structural arms (all five kinds, invalidate-only on `payload.experiment_id`), the `ingest_*` progress arms (broadcast-only), and the `ingestInFlight` store — is **delivered by Phase E1** (E1 is the contract owner and lands first; the arms are pure invalidate-only refetches with no dependency on E2's mutators). E2's Task 7 only VERIFIES E1's structural arms are present and correct.

**Architecture:**

- **E2 IMPORTS E1's contract; it never redefines it.** The roll-up types (`Load`, `LoadSample`, `LoadExposure`, `GroupingFlag`), the `listLoads`/`renameSample`/`moveExposure`/`mergeSamples`/`splitSample`/`dismissGroupingFlag` fetchers, the read hooks (`useLoads`, `useExperiment`, `useUpdateExperiment`), and `queryKeys.loads(id)` are **pinned once in E1's `api.ts`/`queries.ts` and imported here** (spec §8.8: "pinned once in `api.ts`/`queries.ts` by Phase E1 and imported (never redefined) by Phase E2"). The canonical `queryKeys.loads(id)` is **`["experiment", id ?? "none", "loads"]`** (mirrors the live `queryKeys.samples` family — an experiment-prefix invalidation also refreshes loads), NOT a bare `["loads", id]`. The canonical `Load` is the **nested §8.8 shape** (`Load.samples[].exposures[]`, keyed `load_id`; exposures keyed `.id`). E2 adds only the *send-path* artifacts (mutators, structural `OpKind`s, the structural-edit mutation hooks, undo extraction) and the surface components — the `applyRemoteToCache` structural RECEIVE arms are **E1-owned** (E2 only verifies them, Task 7). **If a test or fixture in this plan hard-codes `["loads", id]`, it is wrong — use E1's `queryKeys.loads(id)`.**
- **Reuse, don't rebuild.** The fold tree composes the shipped `Card`, `Checkbox`, `Input variant="title"`, `Thumbnail`, `FlagButton`, `IconButton`, `Menu`, `SearchInput`, `ModalShell`, `EmptyState` primitives. New appearance-carrying leaves (`LoadFold`, `SampleFold`, `ExposureLeaf`, `GeometryLedger`, `SourcesCard`) live under `print/components/` (which is **NOT** design-guard-exempt — only `print/ui/`, `print/plot/`, `print/detector/`, `print/comb/`, `print/export/` are) and carry only placement/layout `className`; any literal appearance they need goes in a `print/ui/` primitive or an existing token utility, so `scripts/check-design.mjs` (the `lint:design` build step) stays green.
- **Undo is a shared hook.** `useUndoStack<T>()` is extracted from `SeriesScopingPage`'s inline `HistoryEntry[]` so grouping-review and scoping share one undo implementation. Single-row edits (rename, one move) carry a session-local re-issue-the-inverse undo via the toast action; multi-row edits (merge, split, bulk-merge) are session-local undo only (and v1 shows the toast WITHOUT an undo action — one stamp cannot reverse a multi-row reassignment server-side, spec §9.3). Undo surfaces through the bottom-center toast — an `action` slot is added to the **imperative `showToast` API** (the `Toast`/`lib/toast.ts` singleton is NO-PROP; there is no `toasts` prop to extend).
- **Mutators wire the optimistic send path.** Five new mutators (`moveExposure`, `renameSample`, `mergeSamples`, `splitSample`, `dismissGroupingFlag`) register in `mutatorRegistry`. `moveExposure`/`renameSample`/`dismissGroupingFlag` splice the `loads(id)` cache surgically in `onMutate`; `mergeSamples`/`splitSample` do an optimistic tree edit but reconcile invalidate-only (ids/rows aren't payload-derivable — the `series_save` invalidate-only precedent in `saveSeries.ts`). **All five MUST invalidate `loads(experimentId)` in their `onSuccess`** (the `saveSeries.ts` precedent): the replay coordinator's own-op path resolves the deferred + aborts the HTTP request and **never calls `applyRemoteToCache`**, so an SSE-arm-only reconcile would leave split's negative-id placeholder unreconciled (a 404 on the next op). The complementary `applyRemoteToCache` structural arms (E1-owned) cover FOREIGN tabs. All five call through `useQueueMutation`, minting the `client_op_id` inside `mutationFn`.
- **The filter keeps a persistent selection.** `GroupingReviewPage` owns selection as an **ordered list** (so "first-selected = survivor" is deterministic for bulk-merge) and the filter/search state; changing the filter never clears selection (the cross-load merge gesture). A floating bulk bar appears when ≥1 is selected; Merge requires ≥2.
- **Configuration internals are presentational + lifted state.** `GeometryLedger` and `SourcesCard` are presentational; `ConfigurationBody` owns the geometry override/undo state and feeds them, and is mounted into E1's `ExperimentConfigurationPage` shell. The acquisition plot reuses the focus/series d3 trace-plot scaffolding (`print/plot/`'s `PlotFrame`/`Axis`/`makeProjection`), NOT a new chart engine.

**Tech Stack:** React 18 + Vite + TypeScript strict (`exactOptionalPropertyTypes`) + TailwindCSS 4. TanStack Query (server state) + Zustand (client state). Vitest + Testing Library (`@testing-library/react`, JSDOM) for unit tests; Playwright for the end-to-end grouping flow. All paths under `packages/HimalayaUI/frontend/`.

**Spec:** `docs/superpowers/specs/2026-06-15-ingestion-redesign-design.md` — §8.2 (grouping-review surface), §8.1 (Configuration tab), §8.6 (component reuse map), §8.7 (new components), §8.8 (the pinned `Load`/`LoadSample`/`LoadExposure`/`GroupingFlag` shapes + `queryKeys.loads(id)` = `["experiment", id ?? "none", "loads"]`), §9.3 (event kinds — **no `sample_merged`; merge fans out as `exposure_moved`; split emits `sample_created` + `sample_split`; `grouping_flag_dismissed` is durable**; merge/split SSE reconcile is invalidate-only; undo is session-local for multi-row), §9.6 (frontend wiring — the hook names `useMoveExposure`/`useRenameSample`/`useMergeSamples`/`useSplitSample`/`useDismissGroupingFlag`, `display_name → name` sweep). Read §8.2, §8.8, and §9.3/§9.6 in full before starting.

**Depends on (HARD — E1 must land first):**
- **Phase E1** (`docs/superpowers/plans/2026-06-18-ingestion-phase-e1-shell-routing.md`): provides `AppRoutes` with `/experiments/:id/{corpus,grouping,config}`, `ExperimentShell`, `ExperimentCorpusPage`, `ExperimentConfigurationPage` (shell only — this plan fills its body), the `ingestInFlight` store + `ingest_*` SSE receive arms, `PageFrame` `home`/`experiment` width keys, **and the pinned `api.ts` additions: the nested `Load`/`LoadSample`/`LoadExposure`/`GroupingFlag` types (§8.8), `listLoads`, `createExperiment`, `triggerScan`, the widened `updateExperiment(id, patch: ExperimentPatch)` (the canonical frontend patch type, renamed from `ExperimentGeometryPatch` — DECISION), the per-field geometry `*_source` columns + `ingest_status`/`scan_signature`/`last_scanned_at` on `Experiment`, `description TEXT` additive column + `image_pattern`/`metadata_pattern`/`integration_pattern` typed columns on `Experiment` (DECISION: additive migration in E1), editable `name`/`description` in `ExperimentShell`, editable pattern rows in the Sources surface, the `Sample.display_name → name` collapse, and the `queryKeys.loads(id)` key + `useLoads`/`useExperiment` read hooks.** E2 imports all of these — it does not redefine any of them.

  > ⚠ **Contract reconciliation note for the implementer (read before Task 2):** E1's plan currently shows a *flat* `Load` in its Task 1 (`{id, experiment_id, load_index, session_id, …}` with **no** `samples`). That contradicts spec §8.8, which pins the **nested** shape (`Load.samples[].exposures[]`, keyed `load_id`). The nested shape is the contract (it mirrors `get_loads_rollup`'s JSON and is what every E2 mutator/component consumes). **E1 must ship the nested §8.8 shape**; if E1 lands the flat one, that is an E1 bug to fix in E1 — do not work around it by redefining `Load` in E2. Coordinate with the E1 owner.

- **Phase D** (`docs/superpowers/plans/2026-06-18-ingestion-phase-d-structural-events.md`): provides the backend routes (`PATCH /api/samples/{id}/name`, `POST /api/exposures/{id}/move`, `POST /api/samples/{id}/merge`, `POST /api/samples/{id}/split`, `POST /api/samples/{id}/dismiss-flag` or equivalent), the `GET /api/experiments/{id}/loads` endpoint, and the event kinds (`exposure_moved`/`sample_renamed`/`sample_created`/`sample_split`/`grouping_flag_dismissed`; **no `sample_merged`**). This plan calls those routes; it does not change the backend. Confirm the exact route paths + the dismiss-flag route shape against Phase D before writing the fetchers (Task 3).
- **The `display_name → name` label sweep** (spec §11): OWNED BY E1 (a single sequenced task lands the `api.ts` `Sample` collapse + `SCHEMA_VERSION` 4→5 bump + the trivial-mutator/`applyRemoteToCache` `update_sample` half). E2's Task 0 is **verification-only**.

**Source of truth for current code:** anchors were line-verified 2026-06-18 but drift — confirm each with a quick `grep`/Read before editing.

---

## File Structure

| File | Responsibility | This plan |
|---|---|---|
| `src/api.ts` | `renameSample`/`moveExposure`/`mergeSamples`/`splitSample`/`dismissGroupingFlag` fetchers (the types `Load`/`LoadSample`/`LoadExposure`/`GroupingFlag` + `listLoads` are E1's — IMPORTED) | MODIFY (Task 3) |
| `src/queries.ts` | structural-edit hooks `useMoveExposure`/`useRenameSample`/`useMergeSamples`/`useSplitSample`/`useDismissGroupingFlag` + `useUpdateExperiment` (E1 ships read hooks `useLoads`/`useExperiment` + `queryKeys.loads` — IMPORTED) | MODIFY (Task 6, 9) |
| `src/lib/queue/types.ts` | add the five structural `OpKind`s | MODIFY (Task 3) |
| `src/lib/queue/mutators/grouping.ts` | `moveExposure`/`renameSample`/`mergeSamples`/`splitSample`/`dismissGroupingFlag` mutators | CREATE (Tasks 4, 5) |
| `src/lib/queue/mutatorRegistry.ts` | register the five mutators (outbound + event) | MODIFY (Tasks 4, 5) |
| `src/lib/queue/applyRemoteToCache.ts` | structural RECEIVE arms (`exposure_moved`/`sample_renamed`/`sample_created`/`sample_split`/`grouping_flag_dismissed`; NO `sample_merged`) — **E1-OWNED**; E2 only VERIFIES (Task 7) | (E1) — E2 verify only |
| `src/hooks/useUndoStack.ts` | shared undo hook (extracted from scoping) | CREATE (Task 1) |
| `src/print/pages/SeriesScopingPage.tsx` | adopt `useUndoStack` | MODIFY (Task 1) |
| `src/lib/toast.ts` + `src/print/ui/Toast.tsx` | optional `action` slot on the imperative `showToast` API + render it | MODIFY (Task 8) |
| `src/print/components/ExposureLeaf.tsx` | one exposure row: thumb · filename · H · time · Move… | CREATE (Task 10) |
| `src/print/components/SampleFold.tsx` | sample accordion: name (edit-in-place) · meta · flag · Rename/Split · merge-prompt · split-divider · leaves | CREATE (Task 11) |
| `src/print/components/LoadFold.tsx` | load accordion: title · meta · status · sample folds | CREATE (Task 12) |
| `src/print/components/GroupingReviewPage.tsx` | the centerpiece: back/summary/filter+search/folds/bulk bar; owns ordered selection + undo | CREATE (Tasks 13, 14, 15) |
| `src/print/components/GroupingBulkBar.tsx` | dark floating bulk-action bar for the persistent selection (local, not a CullBar) | CREATE (Task 13) |
| `src/print/components/GeometryLedger.tsx` | per-field value + provenance chip + Override/Revert + discrepancy banner | CREATE (Task 16) |
| `src/print/components/SourcesCard.tsx` | editable `{role, path}` rows | CREATE (Task 17) |
| `src/print/components/AcquisitionTimeline.tsx` | placement wrapper for the acquisition bars | CREATE (Task 18) |
| `src/print/plot/AcquisitionChart.tsx` | the SVG bar chart (design-guard-exempt render layer) | CREATE (Task 18) |
| `src/print/components/ConfigurationBody.tsx` | composes GeometryLedger + AcquisitionTimeline + SourcesCard; owns override/undo state | CREATE (Task 19) |
| `src/styles.css` | `--color-accent-wash` token (if used by GeometryLedger) | MODIFY (Task 16, conditional) |
| `test/*.test.tsx` / `e2e/*.spec.ts` | per-task tests | CREATE (all tasks) |

**Out of scope (other plans):** E1's routing/shell/`ingestInFlight`/`ingest_*` receive arms/**the `applyRemoteToCache` structural RECEIVE arms (all five structural kinds — E1-owned; E2 only verifies, Task 7)**/`createExperiment` flow/`DirectoryPickerField`/`LiveIngestUnfold`/`ExperimentsHomePage`/`FailedScanPage`/the `StatBar` primitive/the `Dropdown` primitive/the `Load` & geometry-source types/`listLoads`/`queryKeys.loads`/`useLoads`/`useExperiment`/the `display_name → name` collapse + `SCHEMA_VERSION` bump; Phase D's backend; the auto-rescan scheduler.

---

## Conventions for every task

- **Run a single Vitest file** from `packages/HimalayaUI/frontend/`:
  `npm test -- <relative test path>`
  (the repo's hook forces `--run`; Vitest is one-shot). The full suite is slow — run only the touched file per task.
- **Type-check before commit** when a task touches a `.ts(x)` public type:
  `npx tsc --noEmit -p tsconfig.json` (NEVER `tsconfig.build.json` — it emits stray `.js` shadows; see the greenfield stray-`.js` gotcha).
- **Run the design guard** when a task adds/edits a component outside `print/ui/**`:
  `node scripts/check-design.mjs` — must exit 0. Appearance goes in a primitive; `className` is placement-only.
- **Commit after each task** once its test passes. Exact `git add` lists are per task. Never `git add -A`.
- **Don't assert on Tailwind class strings.** Use `data-testid` / `data-*` attributes (frontend AGENTS.md anti-pattern).
- **Mint `client_op_id` inside `mutationFn`**, never at hook construction.

---

## Task 0: VERIFY the `display_name → name` label collapse landed (E1-owned; verification-only)

Spec §11 / §8.8: the `display_name → name` collapse (the `api.ts` `Sample` type, the `lib/queue/mutators/trivial.ts` `UpdateSampleInput`/`patchOf`/`onSuccess`, the `applyRemoteToCache.ts:314` `update_sample` arm, the `updateSample` patch key) **and** the `persistence.ts` `SCHEMA_VERSION` 4→5 bump are landed by **exactly one sequenced E1 task** (spec §8.8: "owned by exactly one sequenced task in E1; E2's collapse step is verification-only"). This task does NOT re-apply the collapse — it only **verifies** E1 landed it before E2's cache writes (which all read/write `s.name`) proceed. If E1 has not landed it, STOP and escalate to the E1 owner; do not race E1 by collapsing the columns here.

**Files:** none modified — verification + (optionally) a guard test that documents the invariant.

- [ ] **Step 1: Verify the collapse is in place**

Run from `packages/HimalayaUI/frontend/`:
```bash
# All four must return ZERO hits if E1 landed the collapse:
grep -n "display_name" src/api.ts
grep -n "display_name" src/lib/queue/mutators/trivial.ts
grep -n "display_name" src/lib/queue/applyRemoteToCache.ts
grep -n "^export const SCHEMA_VERSION" src/lib/queue/persistence.ts   # expect = 5
```
- If `display_name` hits remain, or `SCHEMA_VERSION` is still `4`: **E1 has not landed its label-collapse task. STOP** — coordinate with the E1 owner. Do not proceed; every downstream cache write in this plan assumes `Sample.name: string` (non-null) and no `display_name`.
- If all four are clean: proceed.

> **Live state at plan-authoring time (2026-06-18, E1 UNMERGED):** `src/api.ts:31-38` still declares `Sample.display_name: string | null`; `trivial.ts:41` `UpdateSampleInput = { display_name?; notes? }`; `trivial.ts:86` `onSuccess` still copies `response.display_name`; `applyRemoteToCache.ts:314-318` `update_sample` arm does `qc.setQueryData(queryKeys.sample(id), old => ({ ...old, ...(payload ?? {}) }))` (an untyped spread — the cross-version `display_name` hazard the spec §9.6 flags; the fix is E1's); `persistence.ts:14` `SCHEMA_VERSION = 4`. **All of these are E1's to change**, not E2's.

This task produces no commit (verification-only).

---

## Task 1: Extract `useUndoStack<T>()` and adopt it in SeriesScopingPage

Spec §8.6: the inline `HistoryEntry[]` undo in `SeriesScopingPage.tsx` is extracted to a shared hook so grouping-review and scoping share one implementation. The hook is generic over the entry type, pushes pure entries (no live-state reads — preserves StrictMode-safe single-apply), and exposes `push`, `undo`, `canUndo`, `clear`, and the top entry.

**Files:**
- Create: `src/hooks/useUndoStack.ts`
- Modify: `src/print/pages/SeriesScopingPage.tsx` (replace inline state with the hook)
- Test: `test/useUndoStack.test.ts`

- [ ] **Step 1: Write the failing test**

Create `test/useUndoStack.test.ts`:
```ts
import { describe, it, expect } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { useUndoStack } from "../src/hooks/useUndoStack";

type Entry = { label: string; undo: () => void };

describe("useUndoStack", () => {
  it("starts empty, canUndo=false", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    expect(result.current.canUndo).toBe(false);
    expect(result.current.top).toBeUndefined();
  });

  it("push then pop returns the entry and removes it", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "move", undo: () => {} }));
    expect(result.current.canUndo).toBe(true);
    expect(result.current.top?.label).toBe("move");
    let popped: Entry | undefined;
    act(() => { popped = result.current.pop(); });
    expect(popped?.label).toBe("move");
    expect(result.current.canUndo).toBe(false);
  });

  it("pop on empty returns undefined (no throw)", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    let popped: Entry | undefined = { label: "x", undo: () => {} };
    act(() => { popped = result.current.pop(); });
    expect(popped).toBeUndefined();
  });

  it("clear empties the stack", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "a", undo: () => {} }));
    act(() => result.current.clear());
    expect(result.current.canUndo).toBe(false);
  });

  it("push uses a functional updater (StrictMode double-invoke pushes once per act)", () => {
    const { result } = renderHook(() => useUndoStack<Entry>());
    act(() => result.current.push({ label: "x", undo: () => {} }));
    expect(result.current.depth).toBe(1);
  });
});
```

- [ ] **Step 2: Run test → FAIL**

`npm test -- test/useUndoStack.test.ts`
Expected: FAIL — `useUndoStack` does not exist (module resolution error).

- [ ] **Step 3: Implement the hook**

Create `src/hooks/useUndoStack.ts`:
```ts
import { useCallback, useRef, useState } from "react";

/**
 * A generic LIFO undo stack. Entries are opaque to the hook — the caller
 * decides each entry's shape and provides an `apply` callback at undo time.
 *
 * Every state transition uses a FUNCTIONAL updater so React StrictMode's
 * double-invoke in dev does NOT double-push or double-pop (the SeriesScoping
 * lesson: an impure updater that reads a captured `stack` pushed twice).
 */
export interface UndoStack<T> {
  /** Push one entry. Functional-updater safe (StrictMode double-invoke pushes
   *  once per `act` because the updater is pure). */
  push: (entry: T) => void;
  /**
   * Remove and RETURN the most recent entry (or undefined when empty). The
   * caller runs the inverse effect OUTSIDE any setState — `pop` performs NO
   * side effect of its own, so a StrictMode double-invoke of the internal
   * updater cannot double-apply the caller's effect (the SeriesScoping lesson:
   * never run a side effect inside a setState updater).
   */
  pop: () => T | undefined;
  /** Empty the stack (e.g. on route leave or successful build). */
  clear: () => void;
  canUndo: boolean;
  /** The entry that `pop` would return next, or undefined when empty. */
  top: T | undefined;
  depth: number;
}

export function useUndoStack<T>(): UndoStack<T> {
  const [stack, setStack] = useState<T[]>([]);
  // A ref mirror lets `pop` return the popped entry synchronously without
  // running a side effect inside the setState updater.
  const ref = useRef<T[]>(stack);
  ref.current = stack;

  const push = useCallback((entry: T) => {
    setStack((prev) => [...prev, entry]);
  }, []);

  const pop = useCallback((): T | undefined => {
    const cur = ref.current;
    if (cur.length === 0) return undefined;
    const entry = cur[cur.length - 1] as T;
    setStack((prev) => prev.slice(0, -1)); // pure updater — no side effect
    return entry;
  }, []);

  const clear = useCallback(() => setStack([]), []);

  return {
    push,
    pop,
    clear,
    canUndo: stack.length > 0,
    top: stack.length > 0 ? (stack[stack.length - 1] as T) : undefined,
    depth: stack.length,
  };
}
```

> **Why `pop()` returns the entry instead of an `undo(apply)` that runs the effect inside the updater:** running a side effect inside a `setState` updater is the exact anti-pattern that bit `SeriesScopingPage` (an impure updater that StrictMode double-invoked minted two history entries — see `SeriesScopingPage.tsx:349-352`). Here the updater stays pure: `pop` reads the top off a ref (synchronous, no effect), schedules a pure `slice(0,-1)`, and RETURNS the entry. The consumer then runs the inverse effect outside React's state machinery: `const e = undoStack.pop(); if (e) e.undo();`. The Step-1 test proves a single `pop` removes a single entry; the consumer's effect runs exactly once because the consumer (not the updater) invokes it. **No StrictMode-conditional second variant is needed** — this IS the StrictMode-safe variant.

- [ ] **Step 4: Run test → PASS**

`npm test -- test/useUndoStack.test.ts` → PASS.

- [ ] **Step 5: Adopt in SeriesScopingPage**

Read `SeriesScopingPage.tsx` around the `HistoryEntry` declaration (anchor `:77`) and its `history`/`setHistory` state plus the `undoLast` and "Undo last change" call sites (anchors `:312`, `:492`). Replace the raw `useState<HistoryEntry[]>` with:
```ts
const undoStack = useUndoStack<HistoryEntry>();
```
Replace each `setHistory((h) => [...h, entry])` with `undoStack.push(entry)`; replace the undo handler body with `const e = undoStack.pop(); if (!e) return; /* the existing per-type switch on e.type, calling setRows/setOrder/setLoose */` (the per-type restore effects move OUT of the old `setHistory` updater into the handler body — which is correct: the live source already runs them inside the `setHistory` updater as a workaround, and `pop()` lets them run in the plain handler where they belong); replace `history.length > 0` gates with `undoStack.canUndo`; replace the clear-on-build/reset-effect `setHistory([])` with `undoStack.clear()`. Keep the `HistoryEntry` union type (the four variants `flag`/`reorder`/`value`/`add` at `SeriesScopingPage.tsx:77-87`) defined in the page (it is scoping-specific).

> **StrictMode parity check:** the live `undo` at `SeriesScopingPage.tsx:492-514` currently nests `setRows`/`setOrder`/`setLoose` INSIDE the `setHistory` updater (a deliberate but fragile workaround). After the swap, those restore calls run in the plain handler after `pop()`. Re-verify the existing scoping reorder-undo test still proves one reorder = one undo entry (Step 6).

- [ ] **Step 6: Verify scoping still green + type-check**

```bash
npm test -- test/seriesScoping  # run the existing scoping unit test(s); adjust glob to match
npx tsc --noEmit -p tsconfig.json
```
Expected: PASS, no type errors. (If a scoping test asserted on the inline `history` state via a testid, update it to the equivalent undo affordance — do not weaken the assertion.)

- [ ] **Step 7: Commit**
```bash
git add src/hooks/useUndoStack.ts src/print/pages/SeriesScopingPage.tsx test/useUndoStack.test.ts
git commit -m "refactor(undo): extract useUndoStack; adopt in SeriesScopingPage"
```

---

## Task 2: VERIFY E1's `Load` roll-up types + `listLoads` are importable (no redefinition)

Spec §8.8: the nested `Load`/`LoadSample`/`LoadExposure`/`GroupingFlag` types + `listLoads` are pinned **once in E1's `api.ts`** and imported by E2. **E2 must NOT declare any of them.** This task confirms they exist and match §8.8 so every downstream mutator/component/test can `import { type Load, type LoadSample, type LoadExposure, type GroupingFlag, listLoads } from "../src/api"`.

**Files:** none modified (a compile-only import test is optional).

- [ ] **Step 1: Verify E1 shipped the §8.8 shapes**

Run from `packages/HimalayaUI/frontend/`:
```bash
grep -n "export interface Load\b\|export interface LoadSample\|export interface LoadExposure\|export type GroupingFlag\|export const listLoads" src/api.ts
```
Confirm the shapes match spec §8.8 **exactly** (this is the contract every E2 test fixture must conform to):
```ts
export type GroupingFlag =
  | { kind: "merge"; merge_with_sample_id: number; merge_with_label: string }
  | { kind: "split"; split_at_index: number; jump_from: number; jump_to: number }
  | null;
export interface LoadExposure { id: number; filename: string; horizontal_position: number | null; timestamp: string | null; }
export interface LoadSample {
  sample_id: number; name: string; slot_index: number;
  grouping_source: string; name_source: string; merged_into_id: number | null;
  flag: GroupingFlag; exposures: LoadExposure[];
}
export interface Load {
  load_id: number; load_index: number; session_id: number | null;
  start_time: string | null; end_time: string | null;
  frame_count: number; note: string | null; samples: LoadSample[];
}
```
> **NOTE the exact field names** — every fixture in Tasks 4/5/10–15 uses these: an exposure leaf key is **`id`** (NOT `exposure_id`, NOT `frame_no`/`status`); a sample has **`merged_into_id`** and `grouping_source`/`name_source` typed as plain `string`; a load is keyed **`load_id`** (no top-level `id`/`experiment_id`); `queryKeys.loads(id)` = **`["experiment", id ?? "none", "loads"]`**. The earlier E2-draft fixtures used `exposure_id`/`frame_no`/`status`/top-level `Load.id`/`["loads", id]` — **all wrong**; the §8.8 shapes above are authoritative. The fixtures in this plan have been rewritten to §8.8; if you see a stray `exposure_id` or `["loads",` anywhere, treat it as a bug.
- If E1's shapes are FLAT (no `samples`) or keyed differently, STOP — that is an E1 contract bug to fix in E1 (see the dependency note up top), not to route around here.

This task produces no commit (verification-only).

---

## Task 3: Add the five structural `OpKind`s + the structural-edit fetchers (incl. dismiss)

**Files:**
- Modify: `src/lib/queue/types.ts` (`OpKind` union)
- Modify: `src/api.ts` (`renameSample`, `moveExposure`, `mergeSamples`, `splitSample`, `dismissGroupingFlag` fetchers)
- Test: `test/structuralApi.test.ts`

> **All fetchers use `request<T>(method, path, body?, opts?)`** — `api.ts` has NO `getJSON`/`postJSON`/`patchJSON` helper; `request` is the sole HTTP abstraction and is what threads `X-Username`/`X-Client-Id`/`X-Client-Op-Id` from `AuthOpts` (`api.ts:71-104`, header-build at `:83-87`). Mirror `updateSample` (`request<Sample>("PATCH", …, patch, opts)`).
>
> **Confirm Phase D route paths + the merge/split/dismiss request bodies + responses** before writing — Phase D is the contract. The spec resolves split's event burst as `sample_created` + `sample_split`; the route response for split is assumed to carry `{ new_sample_id }` and for merge `{ loser_id?, survivor_id? }`, but VERIFY against Phase D. The dismiss-flag route path (e.g. `POST /api/samples/{id}/dismiss-flag`) and body (`{ flag_kind, merge_with_sample_id? }`) must match Phase D's `grouping_flag_dismissed` emitter.

- [ ] **Step 1: Write the failing test**

Create `test/structuralApi.test.ts`:
```ts
import { describe, it, expect, vi, afterEach } from "vitest";
import { renameSample, moveExposure, mergeSamples, splitSample, dismissGroupingFlag } from "../src/api";

afterEach(() => vi.restoreAllMocks());

function mockOnce(status = 200, body: unknown = {}) {
  return vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(body), { status, headers: { "Content-Type": "application/json" } }),
  );
}

describe("structural-edit fetchers (use request<T>, thread X-Client-Op-Id)", () => {
  it("renameSample PATCHes /api/samples/:id/name with { name }", async () => {
    const f = mockOnce(200, { id: 1, name: "X", notes: null });
    await renameSample(1, "X", { username: "alice", clientOpId: "op1" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/samples/1/name");
    expect((init as RequestInit).method).toBe("PATCH");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({ name: "X" });
    expect((init as RequestInit).headers).toMatchObject({ "X-Client-Op-Id": "op1" });
  });

  it("moveExposure POSTs /api/exposures/:id/move with { sample_id }", async () => {
    const f = mockOnce();
    await moveExposure(100, 20, { clientOpId: "op2" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/exposures/100/move");
    expect((init as RequestInit).method).toBe("POST");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({ sample_id: 20 });
  });

  it("mergeSamples POSTs /api/samples/:loser/merge with { survivor_id }", async () => {
    const f = mockOnce(200, { loser_id: 5, survivor_id: 3 });
    await mergeSamples(5, 3, { clientOpId: "op3" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/samples/5/merge");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({ survivor_id: 3 });
  });

  it("splitSample POSTs /api/samples/:id/split with { exposure_ids, name }", async () => {
    const f = mockOnce(201, { new_sample_id: 9 });
    await splitSample(1, [101, 102], "HA85 (S01P01b)", { clientOpId: "op4" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/samples/1/split");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({
      exposure_ids: [101, 102], name: "HA85 (S01P01b)",
    });
  });

  it("dismissGroupingFlag POSTs the dismiss route with { flag_kind, merge_with_sample_id? }", async () => {
    const f = mockOnce(200, {});
    await dismissGroupingFlag(20, { flag_kind: "merge", merge_with_sample_id: 10 }, { clientOpId: "op5" });
    const [url, init] = f.mock.calls[0]!;
    expect(String(url)).toContain("/api/samples/20/dismiss-flag"); // confirm exact path vs Phase D
    expect((init as RequestInit).method).toBe("POST");
    expect(JSON.parse(String((init as RequestInit).body))).toEqual({ flag_kind: "merge", merge_with_sample_id: 10 });
  });
});
```

- [ ] **Step 2: Run test → FAIL**

`npm test -- test/structuralApi.test.ts` → FAIL (fetchers undefined).

- [ ] **Step 3: Add the `OpKind`s**

In `src/lib/queue/types.ts`, extend the `OpKind` union (after the existing kinds — confirm the live members at `types.ts:37-58`; do NOT remove any):
```ts
  // Ingestion grouping-review structural edits (Phase E2). `move_exposure`,
  // `rename_sample`, and `dismiss_grouping_flag` are single-entity (surgical
  // cache splice); `merge_samples` and `split_sample` are route-level
  // orchestrations (optimistic tree edit, invalidate-only confirm). The OpKind
  // names the user GESTURE; the wire event kinds are exposure_moved /
  // sample_renamed / sample_created+sample_split / grouping_flag_dismissed —
  // and a merge emits NO sample_merged kind (it fans out as exposure_moved).
  | "move_exposure" | "rename_sample" | "merge_samples" | "split_sample" | "dismiss_grouping_flag"
```

- [ ] **Step 4: Add the fetchers (all via `request<T>`)**

In `src/api.ts` (mirror `updateSample` at `:134`):
```ts
export const renameSample = (id: number, name: string, opts?: AuthOpts): Promise<Sample> =>
  request<Sample>("PATCH", `/api/samples/${id}/name`, { name }, opts);

export const moveExposure = (exposureId: number, sampleId: number, opts?: AuthOpts): Promise<Exposure> =>
  request<Exposure>("POST", `/api/exposures/${exposureId}/move`, { sample_id: sampleId }, opts);

export interface MergeSamplesResponse { loser_id: number; survivor_id: number }
export const mergeSamples = (loserId: number, survivorId: number, opts?: AuthOpts): Promise<MergeSamplesResponse> =>
  request<MergeSamplesResponse>("POST", `/api/samples/${loserId}/merge`, { survivor_id: survivorId }, opts);

export interface SplitSampleResponse { new_sample_id: number }
export const splitSample = (sampleId: number, exposureIds: number[], name: string, opts?: AuthOpts): Promise<SplitSampleResponse> =>
  request<SplitSampleResponse>("POST", `/api/samples/${sampleId}/split`, { exposure_ids: exposureIds, name }, opts);

/** "Keep separate" — durable dismissal of a backend-produced grouping flag
 *  (spec §9.3: grouping_flag_dismissed; suppressed in get_loads_rollup, so it
 *  stays gone across rescans). VERIFY the route path + body against Phase D. */
export interface DismissGroupingFlagBody { flag_kind: "merge" | "split"; merge_with_sample_id?: number }
export const dismissGroupingFlag = (sampleId: number, body: DismissGroupingFlagBody, opts?: AuthOpts): Promise<void> =>
  request<void>("POST", `/api/samples/${sampleId}/dismiss-flag`, body, opts);
```

- [ ] **Step 5: Run test → PASS + type-check**
```bash
npm test -- test/structuralApi.test.ts
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 6: Commit**
```bash
git add src/lib/queue/types.ts src/api.ts test/structuralApi.test.ts
git commit -m "feat(api): structural OpKinds + rename/move/merge/split/dismiss fetchers"
```

---

## Task 4: `moveExposure` + `renameSample` mutators (single-entity, surgical cache splice)

Both are single-entity: they splice the `loads(id)` cache directly in `onMutate` and confirm via SSE. `moveExposure` moves one `LoadExposure` between two `LoadSample`s in the cached tree; `renameSample` rewrites one `LoadSample.name`.

**Files:**
- Create: `src/lib/queue/mutators/grouping.ts` (the two mutators + a shared cache helper)
- Modify: `src/lib/queue/mutatorRegistry.ts` (register both, outbound + event)
- Test: `test/groupingMutators.move-rename.test.ts`

- [ ] **Step 1: Write the failing test**

Create `test/groupingMutators.move-rename.test.ts`:
```ts
import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../src/queries";
import type { Load } from "../src/api";
import { moveExposureMutator, renameSampleMutator } from "../src/lib/queue/mutators/grouping";

function seedLoads(qc: QueryClient, experimentId: number) {
  // §8.8 shape: exposures keyed `.id`; samples carry `merged_into_id`; load keyed `load_id`.
  const loads: Load[] = [{
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 8, note: null,
    samples: [
      { sample_id: 10, name: "A", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
        exposures: [{ id: 100, filename: "a1.tif", horizontal_position: 8, timestamp: null }] },
      { sample_id: 20, name: "B", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
        exposures: [{ id: 200, filename: "b1.tif", horizontal_position: 12, timestamp: null }] },
    ],
  }];
  qc.setQueryData(queryKeys.loads(experimentId), loads);
}

describe("moveExposureMutator", () => {
  it("onMutate moves the exposure between samples in the loads cache; restore reverts", () => {
    const qc = new QueryClient();
    seedLoads(qc, 7);
    // exposureId + sampleId are in the Input (not scope) — one hook per experiment
    const ctx = moveExposureMutator.onMutate(
      { kind: "move_exposure", clientOpId: "op", payload: { exposureId: 100, sampleId: 20 },
        experimentId: 7, exposureId: 100, sampleId: 20,
        username: "alice", clientId: "c" } as never, qc);
    const after = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    const a = after[0]!.samples.find((s) => s.sample_id === 10)!;
    const b = after[0]!.samples.find((s) => s.sample_id === 20)!;
    expect(a.exposures.map((e) => e.id)).toEqual([]);
    expect(b.exposures.map((e) => e.id)).toEqual([200, 100]);
    ctx.restore();
    const reverted = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    expect(reverted[0]!.samples.find((s) => s.sample_id === 10)!.exposures.map((e) => e.id)).toEqual([100]);
  });

  it("onSuccess invalidates loads(experimentId) (own-op reconcile — replay never calls applyRemoteToCache)", () => {
    const qc = new QueryClient();
    seedLoads(qc, 7);
    const inv = vi.spyOn(qc, "invalidateQueries");
    moveExposureMutator.onSuccess(
      { experimentId: 7, exposureId: 100, sampleId: 20 } as never,
      { id: 100 } as never, qc);
    const keys = inv.mock.calls.map((c) => JSON.stringify((c[0] as { queryKey: unknown }).queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });
});

describe("renameSampleMutator", () => {
  it("onMutate rewrites the sample name + sets name_source=user; restore reverts", () => {
    const qc = new QueryClient();
    seedLoads(qc, 7);
    // sampleId is in the Input (not scope) — one hook per experiment
    const ctx = renameSampleMutator.onMutate(
      { kind: "rename_sample", clientOpId: "op", payload: { sampleId: 10, name: "Renamed" },
        experimentId: 7, sampleId: 10, name: "Renamed",
        username: "alice", clientId: "c" } as never, qc);
    const after = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    const a = after[0]!.samples.find((s) => s.sample_id === 10)!;
    expect(a.name).toBe("Renamed");
    expect(a.name_source).toBe("user");
    ctx.restore();
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples.find((s) => s.sample_id === 10)!.name).toBe("A");
  });
});
```

- [ ] **Step 2: Run test → FAIL**

`npm test -- test/groupingMutators.move-rename.test.ts` → FAIL (module/mutators undefined).

- [ ] **Step 3: Implement the mutators**

Create `src/lib/queue/mutators/grouping.ts`:
```ts
import type { QueryClient } from "@tanstack/react-query";
import * as api from "../../../api";
import type { Load, LoadSample } from "../../../api";
import { queryKeys } from "../../../queries";
import { authOpts } from "../../authOpts";          // src/lib/authOpts.ts — the shared helper
import type { Mutator, RollbackContext } from "../types";

// The shared scope every grouping mutator carries (experiment + identity).
// merge/split overlap input and scope; see Task 5.
interface BaseScope {
  username: string | undefined;
  clientId: string;
  experimentId: number;
}

// authOpts(username, clientId, clientOpId) — same idiom as trivial.ts:29-35.
function buildAuthOpts(p: {
  username?: string | undefined;
  clientId?: string | undefined;
  clientOpId?: string | undefined;
}): api.AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

/** Invalidate the loads roll-up for an experiment. Called from EVERY grouping
 *  mutator's onSuccess: the replay coordinator's own-op path resolves the
 *  deferred + aborts the HTTP request and NEVER calls applyRemoteToCache
 *  (replayCoordinator.ts case-1), so the SSE arm alone would not reconcile an
 *  own-op edit — the saveSeries.ts `series_save` precedent invalidates in
 *  onSuccess for exactly this reason. E1's applyRemoteToCache structural arms
 *  (E1-owned; E2 verifies in Task 7) cover FOREIGN tabs. */
function invalidateLoads(qc: QueryClient, experimentId: number): void {
  qc.invalidateQueries({ queryKey: queryKeys.loads(experimentId) });
  qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) });
}

/** Snapshot the loads cache, run `mutate`, return a restore closure. Shared by
 *  every grouping mutator so rollback is uniform. */
function patchLoads(
  qc: QueryClient, experimentId: number, mutate: (loads: Load[]) => Load[],
): RollbackContext {
  const key = queryKeys.loads(experimentId);
  const prev = qc.getQueryData<Load[]>(key);
  if (prev) qc.setQueryData<Load[]>(key, mutate(structuredClone(prev)));
  return { restore: () => { if (prev !== undefined) qc.setQueryData(key, prev); } };
}

// ---------------------------------------------------------------------------
// move_exposure  (single-entity → exposure_moved)
// Entity ids in Input, NOT scope — one hook instance per experiment, not per row.
// ---------------------------------------------------------------------------
type MoveExposureInput = { exposureId: number; sampleId: number };
type MoveExposureScope = BaseScope;

export const moveExposureMutator: Mutator<MoveExposureInput, MoveExposureScope, api.Exposure> = {
  kind: "move_exposure",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      let moved: LoadSample["exposures"][number] | undefined;
      for (const ld of loads) {
        for (const s of ld.samples) {
          const i = s.exposures.findIndex((e) => e.id === p.exposureId); // §8.8: exposure key is `.id`
          if (i >= 0) { moved = s.exposures.splice(i, 1)[0]; break; }
        }
        if (moved) break;
      }
      if (moved) {
        for (const ld of loads) {
          const dest = ld.samples.find((s) => s.sample_id === p.sampleId);
          if (dest) { dest.exposures.push(moved); break; }
        }
      }
      return loads;
    }),
  request: (p) => api.moveExposure(p.exposureId, p.sampleId, buildAuthOpts(p)),
  // Own-op reconcile: the surgical onMutate reflects the end state, but a move
  // can flip a flag (re-derived server-side), and replay never calls
  // applyRemoteToCache for own ops — so invalidate loads(id) to refetch.
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};

// ---------------------------------------------------------------------------
// rename_sample  (single-entity → sample_renamed)
// Entity id in Input, NOT scope — one hook instance per experiment, not per row.
// ---------------------------------------------------------------------------
type RenameSampleInput = { sampleId: number; name: string };
type RenameSampleScope = BaseScope;

export const renameSampleMutator: Mutator<RenameSampleInput, RenameSampleScope, api.Sample> = {
  kind: "rename_sample",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      for (const ld of loads) {
        const s = ld.samples.find((x) => x.sample_id === p.sampleId);
        if (s) { s.name = p.name; s.name_source = "user"; break; }
      }
      return loads;
    }),
  request: (p) => api.renameSample(p.sampleId, p.name, buildAuthOpts(p)),
  onSuccess: (p, response, qc) => {
    // `response.name` may be a synthesized PARTIAL on the own-op SSE-wins path
    // (replayCoordinator resolves the deferred with a frame-derived stub, not
    // the HTTP body — feedback_queue_sse_wins_partial), so guard it before
    // splicing; then invalidate to re-derive flags / refresh the flat samples.
    patchLoads(qc, p.experimentId, (loads) => {
      for (const ld of loads) {
        const s = ld.samples.find((x) => x.sample_id === p.sampleId);
        if (s && response?.name) s.name = response.name;
      }
      return loads;
    });
    invalidateLoads(qc, p.experimentId);
  },
};
```

> **`structuredClone` availability:** Node 18+ and JSDOM provide `structuredClone` globally; Vitest runs on Node 18+. If the project's test env lacks it, add a tiny deep-clone helper instead — but verify first (`node -e "structuredClone({})"`).

- [ ] **Step 4: Register in `mutatorRegistry.ts`**

Add imports and `resolveMutator` cases (the live switch is at `mutatorRegistry.ts:58-134`; insert before `default:`):
```ts
import { moveExposureMutator, renameSampleMutator } from "./mutators/grouping";
// ...in resolveMutator switch:
    case "move_exposure": return moveExposureMutator;
    case "rename_sample": return renameSampleMutator;
```
And `resolveMutatorForEvent` cases (the live switch is at `:147-207`; event kind → mutator, for SSE-wins own-op confirmation). **Note: a merge emits NO `sample_merged` kind** — it fans out as `exposure_moved` frames, so `exposure_moved` is the only kind move/merge share on the wire:
```ts
    case "exposure_moved":  return moveExposureMutator;
    case "sample_renamed":  return renameSampleMutator;
```

- [ ] **Step 5: Run test → PASS + type-check + design guard**
```bash
npm test -- test/groupingMutators.move-rename.test.ts
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 6: Commit**
```bash
git add src/lib/queue/mutators/grouping.ts src/lib/queue/mutatorRegistry.ts test/groupingMutators.move-rename.test.ts
git commit -m "feat(queue): move_exposure + rename_sample mutators (surgical loads splice)"
```

---

## Task 5: `mergeSamples` + `splitSample` + `dismissGroupingFlag` mutators

`merge`/`split` are orchestrations: `onMutate` does an OPTIMISTIC tree edit (merge re-points the loser's exposures onto the survivor and retires the loser; split creates a new sample with a negative placeholder id and moves the chosen exposures into it), and confirmation is **invalidate-only** — the real survivor/new-sample rows aren't payload-derivable (spec §9.3), and the own-op path never calls `applyRemoteToCache`, so `onSuccess` invalidates `loads(id)` itself (the `saveSeries.ts` precedent). `dismissGroupingFlag` is single-entity (surgical splice clears `sample.flag`) and DURABLE (spec §9.3: it writes a `grouping_flag_dismissed` event, suppressed in `get_loads_rollup`, so it stays gone across rescans — NOT session-local `Set` state).

> **`MergeSamplesInput = { loserId; survivorId }` (decided here, ONCE).** Both ids arrive in the per-call mutate() input; the scope carries only experiment + identity (`MergeSamplesScope = BaseScope`). This removes the `as never` cast in the hook (Task 9). `splitSample`'s input is `{ sampleId; exposureIds; name }` — `sampleId` (the source sample to split) moves into Input, not scope, consistent with the hook-arity decision that all row-entity ids travel in the Input. `dismissGroupingFlag`'s input is `{ sampleId; flagKind; mergeWithSampleId? }` for the same reason.

**Files:**
- Modify: `src/lib/queue/mutators/grouping.ts`
- Modify: `src/lib/queue/mutatorRegistry.ts`
- Test: `test/groupingMutators.merge-split.test.ts`, `test/groupingMutators.dismiss.test.ts`

- [ ] **Step 1: Write the failing test**

Create `test/groupingMutators.merge-split.test.ts`:
```ts
import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../src/queries";
import type { Load } from "../src/api";
import { mergeSamplesMutator, splitSampleMutator } from "../src/lib/queue/mutators/grouping";

function seed(qc: QueryClient) {
  const loads: Load[] = [{
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [
      { sample_id: 10, name: "survivor", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
        exposures: [{ id: 100, filename: "a.tif", horizontal_position: 8, timestamp: null }] },
      { sample_id: 20, name: "loser", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
        flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "survivor" },
        exposures: [{ id: 200, filename: "b.tif", horizontal_position: 12, timestamp: null }] },
    ],
  }];
  qc.setQueryData(queryKeys.loads(7), loads);
}

describe("mergeSamplesMutator", () => {
  it("onMutate re-points loser exposures onto survivor and removes the loser; restore reverts", () => {
    const qc = new QueryClient(); seed(qc);
    const ctx = mergeSamplesMutator.onMutate(
      { kind: "merge_samples", clientOpId: "op", payload: { loserId: 20, survivorId: 10 },
        experimentId: 7, loserId: 20, survivorId: 10, username: "a", clientId: "c" } as never, qc);
    const after = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    expect(after[0]!.samples.map((s) => s.sample_id)).toEqual([10]);
    expect(after[0]!.samples[0]!.exposures.map((e) => e.id)).toEqual([100, 200]);
    ctx.restore();
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples.map((s) => s.sample_id)).toEqual([10, 20]);
  });
  it("onSuccess invalidates loads(experimentId)", () => {
    const qc = new QueryClient(); seed(qc);
    const inv = vi.spyOn(qc, "invalidateQueries");
    mergeSamplesMutator.onSuccess({ experimentId: 7, loserId: 20, survivorId: 10 } as never, { loser_id: 20, survivor_id: 10 } as never, qc);
    const keys = inv.mock.calls.map((c) => JSON.stringify((c[0] as { queryKey: unknown }).queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });
});

describe("splitSampleMutator", () => {
  it("onMutate creates a new placeholder sample with the chosen exposures; restore reverts", () => {
    const qc = new QueryClient(); seed(qc);
    // add a second exposure to sample 10 so we can split it
    const loads = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    loads[0]!.samples[0]!.exposures.push({ id: 101, filename: "a2.tif", horizontal_position: 30, timestamp: null });
    qc.setQueryData(queryKeys.loads(7), loads);

    // sampleId is in the Input (not scope) — one hook per experiment, not per row
    const ctx = splitSampleMutator.onMutate(
      { kind: "split_sample", clientOpId: "op", payload: { sampleId: 10, exposureIds: [101], name: "survivorb" },
        experimentId: 7, sampleId: 10, exposureIds: [101], name: "survivorb",
        username: "a", clientId: "c" } as never, qc);
    const after = qc.getQueryData<Load[]>(queryKeys.loads(7))!;
    const src = after[0]!.samples.find((s) => s.sample_id === 10)!;
    expect(src.exposures.map((e) => e.id)).toEqual([100]);
    const created = after[0]!.samples.find((s) => s.sample_id < 0)!;
    expect(created.name).toBe("survivorb");
    expect(created.exposures.map((e) => e.id)).toEqual([101]);
    ctx.restore();
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples.every((s) => s.sample_id > 0)).toBe(true);
  });
});
```

Create `test/groupingMutators.dismiss.test.ts`:
```ts
import { describe, it, expect, vi } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { queryKeys } from "../src/queries";
import type { Load } from "../src/api";
import { dismissGroupingFlagMutator } from "../src/lib/queue/mutators/grouping";

function seed(qc: QueryClient) {
  const loads: Load[] = [{
    load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [{ sample_id: 20, name: "B", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
      flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "A" }, exposures: [] }],
  }];
  qc.setQueryData(queryKeys.loads(7), loads);
}

describe("dismissGroupingFlagMutator (durable 'Keep separate')", () => {
  it("onMutate clears sample.flag optimistically; restore reverts", () => {
    const qc = new QueryClient(); seed(qc);
    // sampleId is in the Input (not scope) — one hook per experiment, not per row
    const ctx = dismissGroupingFlagMutator.onMutate(
      { kind: "dismiss_grouping_flag", clientOpId: "op",
        payload: { sampleId: 20, flagKind: "merge", mergeWithSampleId: 10 },
        experimentId: 7, sampleId: 20, flagKind: "merge", mergeWithSampleId: 10,
        username: "a", clientId: "c" } as never, qc);
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples[0]!.flag).toBeNull();
    ctx.restore();
    expect(qc.getQueryData<Load[]>(queryKeys.loads(7))![0]!.samples[0]!.flag).not.toBeNull();
  });
  it("onSuccess invalidates loads(experimentId)", () => {
    const qc = new QueryClient(); seed(qc);
    const inv = vi.spyOn(qc, "invalidateQueries");
    dismissGroupingFlagMutator.onSuccess({ experimentId: 7, sampleId: 20 } as never, undefined as never, qc);
    const keys = inv.mock.calls.map((c) => JSON.stringify((c[0] as { queryKey: unknown }).queryKey));
    expect(keys).toContain(JSON.stringify(queryKeys.loads(7)));
  });
});
```

- [ ] **Step 2: Run test → FAIL**

`npm test -- test/groupingMutators.merge-split.test.ts test/groupingMutators.dismiss.test.ts` → FAIL (mutators undefined).

- [ ] **Step 3: Implement, appending to `grouping.ts`**

```ts
import { nextOptimisticId } from "../optimisticId";   // negative placeholder id

// ---------------------------------------------------------------------------
// merge_samples  (orchestration; NO sample_merged event — the backend fans it
// out as exposure_moved frames. Optimistic tree edit, invalidate-only confirm.)
// ---------------------------------------------------------------------------
type MergeSamplesInput = { loserId: number; survivorId: number };
type MergeSamplesScope = BaseScope;

export const mergeSamplesMutator: Mutator<MergeSamplesInput, MergeSamplesScope, api.MergeSamplesResponse> = {
  kind: "merge_samples",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      // gather the loser's exposures, then drop the loser sample, then append
      // those exposures onto the survivor (preserving order).
      let loserExps: LoadSample["exposures"] = [];
      for (const ld of loads) {
        const li = ld.samples.findIndex((s) => s.sample_id === p.loserId);
        if (li >= 0) { loserExps = ld.samples[li]!.exposures; ld.samples.splice(li, 1); break; }
      }
      for (const ld of loads) {
        const surv = ld.samples.find((s) => s.sample_id === p.survivorId);
        if (surv) {
          surv.exposures.push(...loserExps);
          surv.grouping_source = "user_merged";
          break;
        }
      }
      return loads;
    }),
  request: (p) => api.mergeSamples(p.loserId, p.survivorId, buildAuthOpts(p)),
  // Invalidate-only confirm: the re-derived flags aren't response-derivable, and
  // the own-op path never calls applyRemoteToCache — so refetch loads(id) here.
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};

// ---------------------------------------------------------------------------
// split_sample  (orchestration → sample_created + sample_split; invalidate-only)
// sampleId (source) in Input, NOT scope — one hook per experiment, not per row.
// ---------------------------------------------------------------------------
type SplitSampleInput = { sampleId: number; exposureIds: number[]; name: string };
type SplitSampleScope = BaseScope;

export const splitSampleMutator: Mutator<SplitSampleInput, SplitSampleScope, api.SplitSampleResponse> = {
  kind: "split_sample",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      const ids = new Set(p.exposureIds);
      for (const ld of loads) {
        const src = ld.samples.find((s) => s.sample_id === p.sampleId);
        if (!src) continue;
        const moved = src.exposures.filter((e) => ids.has(e.id));        // §8.8: `.id`
        src.exposures = src.exposures.filter((e) => !ids.has(e.id));
        const created: LoadSample = {
          sample_id: nextOptimisticId(), // NEGATIVE placeholder until SSE confirms
          name: p.name, slot_index: src.slot_index, grouping_source: "manual",
          name_source: "user", merged_into_id: null, flag: null, exposures: moved,
        };
        const srcIdx = ld.samples.indexOf(src);
        ld.samples.splice(srcIdx + 1, 0, created);
        break;
      }
      return loads;
    }),
  request: (p) => api.splitSample(p.sampleId, p.exposureIds, p.name, buildAuthOpts(p)),
  // Invalidate-only: the real new sample id arrives only via the refetch; the
  // own-op path never calls applyRemoteToCache, so we must invalidate here or
  // the negative-id placeholder is never reconciled (404 on a follow-up op).
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};

// ---------------------------------------------------------------------------
// dismiss_grouping_flag  (single-entity, DURABLE → grouping_flag_dismissed)
// sampleId in Input, NOT scope — one hook per experiment, not per row.
// ---------------------------------------------------------------------------
type DismissFlagInput = { sampleId: number; flagKind: "merge" | "split"; mergeWithSampleId?: number };
type DismissFlagScope = BaseScope;

export const dismissGroupingFlagMutator: Mutator<DismissFlagInput, DismissFlagScope, void> = {
  kind: "dismiss_grouping_flag",
  onMutate: (p, qc): RollbackContext =>
    patchLoads(qc, p.experimentId, (loads) => {
      for (const ld of loads) {
        const s = ld.samples.find((x) => x.sample_id === p.sampleId);
        if (s) { s.flag = null; break; } // optimistic: clear the suggestion
      }
      return loads;
    }),
  request: (p) => api.dismissGroupingFlag(
    p.sampleId,
    p.mergeWithSampleId !== undefined
      ? { flag_kind: p.flagKind, merge_with_sample_id: p.mergeWithSampleId }
      : { flag_kind: p.flagKind },
    buildAuthOpts(p),
  ),
  onSuccess: (p, _response, qc) => invalidateLoads(qc, p.experimentId),
};
```

> **`exactOptionalPropertyTypes` note:** `merge_with_sample_id?` must be conditionally spread (as above), never set to `undefined`, or `request<void>` will type-error on the optional body field.

- [ ] **Step 4: Register in `mutatorRegistry.ts`**

```ts
import { mergeSamplesMutator, splitSampleMutator, dismissGroupingFlagMutator } from "./mutators/grouping";
// resolveMutator (before default:):
    case "merge_samples":        return mergeSamplesMutator;
    case "split_sample":         return splitSampleMutator;
    case "dismiss_grouping_flag": return dismissGroupingFlagMutator;
// resolveMutatorForEvent — there is NO sample_merged kind (merge fans out as
// exposure_moved, already mapped to moveExposureMutator in Task 4). Split emits
// sample_created (the new sample) + sample_split (the source). For own-op SSE
// confirmation, sample_split is the source-sample frame that confirms a split:
    case "sample_split":  return splitSampleMutator;
    case "grouping_flag_dismissed": return dismissGroupingFlagMutator;
    // sample_created has NO own-op mutator (it's the create half of split; the
    // foreign-tab arm in Task 7 invalidate-only refreshes loads) — leave it to
    // resolveMutatorForEvent's default (undefined) so it routes through the
    // applyRemoteToCache sample_created arm, not a mutator rerun.
```

> **Verify against `replayCoordinator.ts`:** for an own-op split, the backend emits TWO frames (`sample_created` then `sample_split`) sharing the op's `client_op_id`. The replay coordinator resolves the deferred on the FIRST matching `client_op_id` frame and aborts the HTTP request. Confirm which frame carries the `client_op_id` for split (Phase D), and that `resolveMutatorForEvent` returns `splitSampleMutator` for whichever frame the deferred resolves on. If both frames carry it, the second is handled foreign-style (it won't find a deferred) and falls to the `applyRemoteToCache` arm — which is invalidate-only, so it's harmless. **This is the one place the merge/split SSE wiring is subtle; trace it in Phase D before finalizing.**

- [ ] **Step 5: Run test → PASS + type-check**
```bash
npm test -- test/groupingMutators.merge-split.test.ts test/groupingMutators.dismiss.test.ts
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 6: Commit**
```bash
git add src/lib/queue/mutators/grouping.ts src/lib/queue/mutatorRegistry.ts test/groupingMutators.merge-split.test.ts test/groupingMutators.dismiss.test.ts
git commit -m "feat(queue): merge/split/dismiss mutators (optimistic tree, onSuccess invalidate-only)"
```

---

## Task 6: VERIFY E1's `queryKeys.loads(id)` + `useLoads` (E1-owned; verification-only)

`queryKeys.loads(id)` and `useLoads` are pinned in E1's `queries.ts` (spec §8.8/§9.6). E2 imports them; it does NOT add them. The canonical key is **`["experiment", id ?? "none", "loads"]`** (the live `queryKeys.samples` family is `["experiment", id ?? "none", "samples"]` — `loads` mirrors it so an experiment-prefix invalidation refreshes loads too). The earlier E2-draft `["loads", id]` is **dropped** (spec §8.8 explicit: "E2's earlier `["loads", id]` is dropped").

**Files:** none modified.

- [ ] **Step 1: Verify**

Run from `packages/HimalayaUI/frontend/`:
```bash
grep -n "loads:" src/queries.ts            # expect: loads: (id) => ["experiment", id ?? "none", "loads"]
grep -n "export function useLoads" src/queries.ts
```
Confirm `queryKeys.loads(7)` deep-equals `["experiment", 7, "loads"]` and `queryKeys.loads(undefined)` deep-equals `["experiment", "none", "loads"]`. If E1 shipped `["loads", id]`, STOP — that's an E1 contract bug (see the dependency note); every E2 cache arm + mutator + fixture in this plan assumes the experiment-prefixed key.

This task produces no commit (verification-only).

---

## Task 7: VERIFY E1's `applyRemoteToCache` structural arms (E1-owned RECEIVE path; verification-only)

**The `applyRemoteToCache` structural RECEIVE arms are E1-owned, NOT E2's.** E1's Task 19 declares ownership of "the SSE RECEIVE path for the structural event kinds" and adds all five arms (`exposure_moved`/`sample_renamed`/`sample_created`/`sample_split`/`grouping_flag_dismissed`) to `src/lib/queue/applyRemoteToCache.ts`. E1 lands first (the contract owner), the arms are pure invalidate-only refetches keyed on `payload.experiment_id` with **no dependency on E2's mutators**, and no structural SSE frame can fire until E2's mutations ship — so there is no ordering hazard. **E2 does NOT modify `applyRemoteToCache.ts` for structural events** (that would be a double-edit of the same arms). This task is a verification-only STOP-gate (mirrors Tasks 0/2/6).

> **What E2 still owns on the RECEIVE side:** nothing in `applyRemoteToCache`. E2's reconcile for OWN ops is the `onSuccess` `loads(experimentId)` invalidation in each of the five mutators (Tasks 4/5) — because `replayCoordinator.ts` resolves the deferred + aborts the HTTP and **never calls `applyRemoteToCache` for own ops**. E1's arms cover the FOREIGN-tab/cross-user case. The two paths are complementary and independently owned.

**Files:** none modified by E2 (verification only).

- [ ] **Step 1: Verify E1's Task 19 landed the five structural arms**

Run from `packages/HimalayaUI/frontend/`:
```bash
grep -n "exposure_moved\|sample_renamed\|sample_created\|sample_split\|grouping_flag_dismissed\|sample_merged" src/lib/queue/applyRemoteToCache.ts
```
Confirm ALL of the following against the live file:
- All FIVE kinds — `exposure_moved`, `sample_renamed`, `sample_created`, `sample_split`, `grouping_flag_dismissed` — have a `case` arm.
- **NO `sample_merged` arm exists** (there is no such event — a foreign merge arrives as an `exposure_moved` burst). If E1 added one, that is an E1 bug — escalate.
- Each arm is **invalidate-only** on `queryKeys.loads(payload.experiment_id)` (+ the flat `queryKeys.samples(payload.experiment_id)` listing), reads **`payload.experiment_id`** (NOT `remote.entity_id`, which `applyRemoteToCache.ts:105` reads unconditionally as the sample/exposure id), and is placed **BEFORE the `default:` arm** (`:334`, which fires the poisoning `peaks(id)`/`indices(id)` invalidations).
- If any kind is missing, shaped wrong (surgical instead of invalidate-only, or scoped on `entity_id`), or placed after `default:`: **STOP and escalate to the E1 owner.** Do not add or fix the arm here — it is E1's.

This task produces no commit (verification-only).

---

## Task 8: Add an `action` slot to the imperative toast API (undo)

The undo toast needs an `{label, onClick}` action. **`ToastContainer` is a NO-PROP imperative singleton** (`Toast.tsx:48` `export function ToastContainer(): JSX.Element` — it holds its own `useState<ToastItem[]>` and registers `showToast` via `setToastImpl` on mount). Toasts are pushed through the imperative API in **`src/lib/toast.ts`**: `showToast(msg: string, kind: ToastKind = "info")` (`toast.ts:12-14`), `setToastImpl(impl)` (`:18-31`), `ToastItem = { id: number; msg: string; kind: ToastKind }` (`Toast.tsx:5-9`), `ToastKind = "info" | "success" | "warning" | "error"` (`toast.ts:5`). There is **no `toasts` prop** to extend — the fix extends the imperative signature.

**Files:**
- Modify: `src/lib/toast.ts` (the `showToast`/`setToastImpl` signatures + an exported `ToastAction` type), `src/print/ui/Toast.tsx` (the `ToastItem` shape + render)
- Test: `test/toastAction.test.tsx`

- [ ] **Step 1: Write the failing test (imperative — render the singleton, then `showToast`)**

Create `test/toastAction.test.tsx`:
```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { ToastContainer } from "../src/print/ui/Toast";
import { showToast } from "../src/lib/toast";

describe("imperative toast action slot", () => {
  it("renders an action button and fires onClick", () => {
    render(<ToastContainer />); // singleton — registers showToast on mount
    const onUndo = vi.fn();
    act(() => { showToast("Moved a1.tif", "info", { label: "Undo", onClick: onUndo }); });
    const btn = screen.getByRole("button", { name: "Undo" });
    fireEvent.click(btn);
    expect(onUndo).toHaveBeenCalledTimes(1);
  });

  it("renders no action button when no action is passed", () => {
    render(<ToastContainer />);
    act(() => { showToast("Saved", "success"); });
    expect(screen.queryByRole("button")).toBeNull();
  });
});
```

> Read `src/lib/toast.ts` + `src/print/ui/Toast.tsx` first and align the test to the LIVE field names (`msg`/`kind`, not `message`/`tone`). The only NEW thing is the optional third `action` arg to `showToast` and its render.

- [ ] **Step 2: Run test → FAIL**

`npm test -- test/toastAction.test.tsx` → FAIL (`showToast` takes 2 args; no action renders).

- [ ] **Step 3: Extend the imperative API + render**

In `src/lib/toast.ts`:
```ts
export interface ToastAction { label: string; onClick: () => void }

let activeImpl: (msg: string, kind: ToastKind, action?: ToastAction) => void =
  (msg, kind) => { console.warn(`[toast:${kind}] ${msg}`); };

export function showToast(msg: string, kind: ToastKind = "info", action?: ToastAction): void {
  activeImpl(msg, kind, action);
}

export function setToastImpl(
  impl: ((msg: string, kind: ToastKind, action?: ToastAction) => void) | null,
): void {
  activeImpl = impl ?? ((msg, kind) => { console.warn(`[toast:${kind}] ${msg}`); });
}
```
In `src/print/ui/Toast.tsx`, extend `ToastItem` and thread `action` through the registered impl + the push:
```ts
import type { ToastAction } from "../../lib/toast"; // confirm relative path
interface ToastItem { id: number; msg: string; kind: ToastKind; action?: ToastAction }
```
The `setToastImpl` callback registered on mount becomes `(msg, kind, action) => setItems((xs) => [...xs, { id: nextId(), msg, kind, action }])` (mirror the live push, adding `action`). In the per-toast render, after the message, add (this file is `print/ui/**` → appearance is permitted; reuse the file's existing inner-button idiom if one exists):
```tsx
{t.action ? (
  <button type="button" className="ml-3 text-xs font-bold text-accent" onClick={t.action.onClick}>
    {t.action.label}
  </button>
) : null}
```

- [ ] **Step 4: Run test → PASS + design guard**
```bash
npm test -- test/toastAction.test.tsx
node scripts/check-design.mjs
```
→ PASS, guard exit 0.

- [ ] **Step 5: Commit**
```bash
git add src/lib/toast.ts src/print/ui/Toast.tsx test/toastAction.test.tsx
git commit -m "feat(ui): imperative showToast action slot (undo)"
```

---

## Task 9: Structural-edit + experiment-update query hooks (spec §9.6 names)

Thin `useQueueMutation` wrappers binding each mutator to its scope (experimentId + username + module-level `CLIENT_ID`), mirroring the live `useUpdateSample`/`useAddPeak` hooks. Plus `useUpdateExperiment` (E1 ships only READ hooks for experiments; the geometry-override mutation is E2's — spec §9.6).

> **Live plumbing (verified):** hooks read `const username = useAppState((s) => s.username)` and pass the **module-level `CLIENT_ID`** const (`queries.ts:42` `const CLIENT_ID = getClientId()`) — there is **no `useClientId()` hook**. `getClientId` lives in `src/lib/clientId.ts` (a function, sessionStorage-backed). Scope objects are passed inline to `useQueueMutation` (no `authOpts(...)` helper at the hook layer). See `useUpdateSample` (`queries.ts:520-526`) and `useAddPeak` (`:372-382`) as templates.

**Files:**
- Modify: `src/queries.ts`
- Test: `test/structuralHooks.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/structuralHooks.test.tsx`:
```tsx
import { describe, it, expect } from "vitest";
import { renderHook } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import {
  useRenameSample, useMoveExposure, useMergeSamples, useSplitSample,
  useDismissGroupingFlag, useUpdateExperiment,
} from "../src/queries";

function wrapper({ children }: { children: ReactNode }) {
  const qc = new QueryClient();
  return <QueryClientProvider client={qc}>{children}</QueryClientProvider>;
}

describe("structural-edit + experiment hooks expose mutate()", () => {
  // All row-scoped hooks take only experimentId; entity ids (sampleId/exposureId)
  // go in the mutate() input (one hook instance per experiment, not per row).
  it("useRenameSample", () => { expect(typeof renderHook(() => useRenameSample(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useMoveExposure", () => { expect(typeof renderHook(() => useMoveExposure(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useMergeSamples", () => { expect(typeof renderHook(() => useMergeSamples(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useSplitSample", () => { expect(typeof renderHook(() => useSplitSample(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useDismissGroupingFlag", () => { expect(typeof renderHook(() => useDismissGroupingFlag(7), { wrapper }).result.current.mutate).toBe("function"); });
  it("useUpdateExperiment", () => { expect(typeof renderHook(() => useUpdateExperiment(7), { wrapper }).result.current.mutate).toBe("function"); });
});
```

- [ ] **Step 2: Run test → FAIL**

`npm test -- test/structuralHooks.test.tsx` → FAIL (hooks undefined).

- [ ] **Step 3: Implement**

In `src/queries.ts` (mirror `useUpdateSample`/`useAddPeak`):
```ts
import {
  moveExposureMutator, renameSampleMutator, mergeSamplesMutator,
  splitSampleMutator, dismissGroupingFlagMutator,
} from "./lib/queue/mutators/grouping";
// CLIENT_ID, useAppState, useQueueMutation, useQueryClient, api are already in queries.ts.

// Row-scoped hooks take ONLY experimentId.
// The entity id (sampleId / exposureId) is carried in the mutate() call input,
// NOT in the hook's scope arg — one hook instance per experiment, not per row.
// This resolves the Tasks 9 ↔ 14 arity conflict: Task 14 can instantiate each hook
// once and pass the per-row id at call time, with no per-row hook construction.

export function useRenameSample(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(renameSampleMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useMoveExposure(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(moveExposureMutator, { experimentId, username, clientId: CLIENT_ID });
}

// merge: scope = experiment + identity ONLY; the caller passes { loserId, survivorId }
// as the mutate() input (MergeSamplesInput, decided in Task 5) — no `as never`.
export function useMergeSamples(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(mergeSamplesMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useSplitSample(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(splitSampleMutator, { experimentId, username, clientId: CLIENT_ID });
}

export function useDismissGroupingFlag(experimentId: number) {
  const username = useAppState((s) => s.username);
  return useQueueMutation(dismissGroupingFlagMutator, { experimentId, username, clientId: CLIENT_ID });
}

/** Geometry/name override for the Configuration tab (spec §9.6 — E1 ships only
 *  read hooks for experiments). Wraps E1's `updateExperiment(id, patch)` fetcher;
 *  this is a plain TanStack mutation (not a queue mutator — geometry override is
 *  not part of the structural-edit event log) that invalidates the experiment
 *  detail on success so the provenance chips re-render. */
export function useUpdateExperiment(experimentId: number) {
  const username = useAppState((s) => s.username);
  const qc = useQueryClient();
  return useMutation({
    mutationFn: (patch: api.ExperimentPatch) =>
      api.updateExperiment(experimentId, patch, { username, clientId: CLIENT_ID }),
    onSuccess: () => {
      qc.invalidateQueries({ queryKey: queryKeys.experiment(experimentId) });
    },
  });
}
```

> **Confirm `useMutation`/`useQueryClient` imports** are present in `queries.ts` (TanStack) — add to the existing `@tanstack/react-query` import if not. `api.ExperimentPatch` is E1's (DECISION: renamed from `ExperimentGeometryPatch`; covers name/description/geometry ×6/patterns ×3). `updateExperiment` is E1's widened fetcher.

- [ ] **Step 4: Run test → PASS + type-check (no `as never` left)**
```bash
npm test -- test/structuralHooks.test.tsx test/groupingMutators.merge-split.test.ts test/groupingMutators.dismiss.test.ts
npx tsc --noEmit -p tsconfig.json
```
→ PASS. `tsc` must accept all hooks WITHOUT an `as never` (all entity ids are in the mutate() inputs, not the scope).

- [ ] **Step 5: Commit**
```bash
git add src/queries.ts test/structuralHooks.test.tsx
git commit -m "feat(queries): structural-edit hooks (spec §9.6 names) + useUpdateExperiment"
```

---

## Task 10: `ExposureLeaf` component

The leaf row: detector thumbnail · filename (mono) · stage H-position · time · a Move… button. Composes `Thumbnail` + `IconButton`/`Menu` for Move. Presentational; the page passes the move targets + onMove.

**Files:**
- Create: `src/print/components/ExposureLeaf.tsx`
- Test: `test/ExposureLeaf.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/ExposureLeaf.test.tsx`:
```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { ExposureLeaf } from "../src/print/components/ExposureLeaf";
import type { LoadExposure } from "../src/api";

// §8.8: LoadExposure = { id, filename, horizontal_position, timestamp } — NO frame_no/status.
const exp: LoadExposure = {
  id: 100, filename: "HA_0241_001.tif", horizontal_position: 8.4, timestamp: "10:02:00",
};

describe("ExposureLeaf", () => {
  it("shows filename, H-position, and time", () => {
    render(<ExposureLeaf exposure={exp} thumbSrc={null} onMove={() => {}} />);
    expect(screen.getByText("HA_0241_001.tif")).toBeInTheDocument();
    expect(screen.getByText(/8\.4/)).toBeInTheDocument();
    expect(screen.getByText(/10:02/)).toBeInTheDocument();
  });
  it("Move button calls onMove with the exposure id", () => {
    const onMove = vi.fn();
    render(<ExposureLeaf exposure={exp} thumbSrc={null} onMove={onMove} />);
    fireEvent.click(screen.getByRole("button", { name: /move/i }));
    expect(onMove).toHaveBeenCalledWith(100);
  });
  it("renders an em-dash when horizontal_position is null", () => {
    render(<ExposureLeaf exposure={{ ...exp, horizontal_position: null }} thumbSrc={null} onMove={() => {}} />);
    expect(screen.getByTestId("exposure-leaf").textContent).toContain("—");
  });
});
```

- [ ] **Step 2: Run test → FAIL** (`npm test -- test/ExposureLeaf.test.tsx`) — component undefined.

- [ ] **Step 3: Implement**

Create `src/print/components/ExposureLeaf.tsx`:
```tsx
import type { JSX } from "react";
import type { LoadExposure } from "../../api";
import { Thumbnail } from "./Thumbnail";
import { IconButton } from "../ui/IconButton";

export interface ExposureLeafProps {
  exposure: LoadExposure;
  /** Detector thumbnail URL for this exposure, or null (Thumbnail draws a
   *  placeholder). The page derives it; Thumbnail takes `src`, NOT an id. */
  thumbSrc: string | null;
  /** Opens the Move… picker for this exposure. */
  onMove: (exposureId: number) => void;
  className?: string;
}

/** One exposure row in the grouping fold: thumb · filename · H · time · Move.
 *  Appearance lives in the composed primitives; this is layout only. */
export function ExposureLeaf({ exposure, thumbSrc, onMove, className }: ExposureLeafProps): JSX.Element {
  const h = exposure.horizontal_position;
  return (
    <div
      data-testid="exposure-leaf"
      className={`flex items-center gap-3 rounded-sm px-2.5 py-1.5${className ? ` ${className}` : ""}`}
    >
      <Thumbnail src={thumbSrc} size="xs" frameNo={exposure.id} />
      <span className="min-w-0 flex-1 truncate font-mono text-xs text-ink">{exposure.filename}</span>
      <span className="w-28 shrink-0 font-mono text-xs text-ink-soft">
        {h === null ? "—" : `H ${h.toFixed(1)} mm`}
      </span>
      <span className="w-20 shrink-0 font-mono text-xs text-ink-faint">{exposure.timestamp ?? "—"}</span>
      <IconButton label="Move to another sample" onClick={() => onMove(exposure.id)}>⋯</IconButton>
    </div>
  );
}
```

> **Primitive APIs (verified):** `Thumbnail` (`print/components/Thumbnail.tsx`) takes **`src: string | null`** + `size?: "xs" | "sm" | "lg"` + `frameNo` — it does NOT take an `exposureId`. The page must compute `thumbSrc` (e.g. the detector image URL for `exposure.id`) and pass it in. `IconButton` (`print/ui/IconButton.tsx`) requires **`label`** (NOT `aria-label` — it `Omit`s `aria-label` and renders `label` as the accessible name) and takes the glyph as `children`. The test asserts behaviour (filename text, onMove call, em-dash) + the Move button's accessible name (`/move/i` matches the `label`), not classes. A flex row (not `grid-cols-[…]`) keeps layout simple and guard-clean. Run `node scripts/check-design.mjs`.

- [ ] **Step 4: Run test → PASS + design guard**
```bash
npm test -- test/ExposureLeaf.test.tsx
node scripts/check-design.mjs
```
→ PASS / exit 0. (If the guard flags `grid-cols-[…]`, swap to a `<div className="flex …">` row layout or a shared grid const — the alignment doesn't need an arbitrary template; mirror however SheetTable lays out its columns.)

- [ ] **Step 5: Commit**
```bash
git add src/print/components/ExposureLeaf.tsx test/ExposureLeaf.test.tsx
git commit -m "feat(grouping): ExposureLeaf row"
```

---

## Task 11: `SampleFold` component

The sample accordion: edit-in-place serif name (`Input variant="title"`), meta (`N exposures · time`), a flag chip, a hover/open `Rename`/`Split…` action cluster, an optional merge-prompt (for `flag.kind==="merge"`), a split-divider (for `flag.kind==="split"`, before the split index), and the exposure leaves. A bulk-select `Checkbox`. Presentational; the page owns selection + all callbacks.

**Files:**
- Create: `src/print/components/SampleFold.tsx`
- Test: `test/SampleFold.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/SampleFold.test.tsx`:
```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { SampleFold } from "../src/print/components/SampleFold";
import type { LoadSample } from "../src/api";

function sample(over: Partial<LoadSample> = {}): LoadSample {
  return {
    sample_id: 10, name: "HA85 (S01P15)", slot_index: 15,
    grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
    exposures: [
      { id: 100, filename: "a1.tif", horizontal_position: 8, timestamp: "10:00" },
      { id: 101, filename: "a2.tif", horizontal_position: 36, timestamp: "10:01" },
    ],
    ...over,
  };
}

const noop = () => {};
const baseProps = {
  open: true, selected: false,
  onToggleOpen: noop, onToggleSelect: noop, onRename: noop,
  onSplit: noop, onMerge: noop, onDismissFlag: noop, onMoveExposure: noop,
  thumbSrcFor: () => null,
};

describe("SampleFold", () => {
  it("renders name + exposure count and the leaves when open", () => {
    render(<SampleFold sample={sample()} {...baseProps} />);
    expect(screen.getByText("HA85 (S01P15)")).toBeInTheDocument();
    expect(screen.getByText(/2 exposures/)).toBeInTheDocument();
    expect(screen.getAllByTestId("exposure-leaf")).toHaveLength(2);
  });

  it("checkbox toggles selection", () => {
    const onToggleSelect = vi.fn();
    render(<SampleFold sample={sample()} {...baseProps} onToggleSelect={onToggleSelect} />);
    fireEvent.click(screen.getByRole("checkbox"));
    expect(onToggleSelect).toHaveBeenCalledWith(10);
  });

  it("merge flag shows the merge prompt with Merge / Keep separate", () => {
    const onMerge = vi.fn(); const onDismiss = vi.fn();
    render(<SampleFold sample={sample({ flag: { kind: "merge", merge_with_sample_id: 5, merge_with_label: "HA85 (S01P08)" } })}
      {...baseProps} onMerge={onMerge} onDismissFlag={onDismiss} />);
    expect(screen.getByText(/HA85 \(S01P08\)/)).toBeInTheDocument();
    fireEvent.click(screen.getByRole("button", { name: /merge/i }));
    expect(onMerge).toHaveBeenCalledWith(10, 5);
    fireEvent.click(screen.getByRole("button", { name: /keep separate/i }));
    expect(onDismiss).toHaveBeenCalledWith(10);
  });

  it("split flag shows a divider with the position jump before the split index", () => {
    render(<SampleFold sample={sample({ flag: { kind: "split", split_at_index: 1, jump_from: 8, jump_to: 36 } })} {...baseProps} />);
    const divider = screen.getByTestId("split-divider");
    expect(divider.textContent).toMatch(/8\.0 → 36\.0/);
  });

  it("Rename action calls onRename", () => {
    const onRename = vi.fn();
    render(<SampleFold sample={sample()} {...baseProps} onRename={onRename} />);
    fireEvent.click(within(screen.getByTestId("sample-fold")).getByRole("button", { name: /^rename$/i }));
    expect(onRename).toHaveBeenCalledWith(10);
  });
});
```

- [ ] **Step 2: Run test → FAIL** (`npm test -- test/SampleFold.test.tsx`).

- [ ] **Step 3: Implement**

Create `src/print/components/SampleFold.tsx`:
```tsx
import type { JSX } from "react";
import type { LoadSample } from "../../api";
import { Checkbox } from "../ui/Checkbox";
import { Button } from "../ui/Button";
import { ExposureLeaf } from "./ExposureLeaf";

export interface SampleFoldProps {
  sample: LoadSample;
  open: boolean;
  selected: boolean;
  onToggleOpen: (sampleId: number) => void;
  onToggleSelect: (sampleId: number) => void;
  onRename: (sampleId: number) => void;
  onSplit: (sampleId: number) => void;
  /** Merge this sample (loser) into the flagged partner (survivor). */
  onMerge: (loserId: number, survivorId: number) => void;
  onDismissFlag: (sampleId: number) => void;
  onMoveExposure: (sampleId: number, exposureId: number) => void;
  /** Per-exposure thumbnail URL resolver (page supplies; Thumbnail takes a src). */
  thumbSrcFor: (exposureId: number) => string | null;
  className?: string;
}

export function SampleFold(props: SampleFoldProps): JSX.Element {
  const { sample: s, open, selected } = props;
  const flag = s.flag;
  const flagChip =
    flag?.kind === "split" ? "check split" :
    flag?.kind === "merge" ? "possible reshoot" : null;

  return (
    <div data-testid="sample-fold" className={`rounded-sm border border-hair bg-plate${flag ? " border-warning" : ""}${selected ? " ring-1 ring-accent" : ""}${props.className ? ` ${props.className}` : ""}`}>
      <div className="flex items-center gap-3 px-3.5 py-2.5">
        <Checkbox checked={selected} onChange={() => props.onToggleSelect(s.sample_id)} aria-label={`Select ${s.name}`} />
        <button type="button" className="min-w-0 flex-1 text-left" onClick={() => props.onToggleOpen(s.sample_id)}>
          <span className="text-headline text-ink">{s.name}</span>
          <span className="ml-2 font-mono text-xs text-ink-faint">{s.exposures.length} exposures · {s.exposures[0]?.timestamp ?? "—"}</span>
          {flagChip ? <span className="ml-2 text-xs font-bold uppercase text-warning">{flagChip}</span> : null}
        </button>
        <div className="flex items-center gap-1.5">
          <Button variant="outline" onClick={() => props.onRename(s.sample_id)}>Rename</Button>
          <Button variant="outline" onClick={() => props.onSplit(s.sample_id)}>Split…</Button>
        </div>
      </div>

      {open ? (
        <div className="border-t border-hair p-1">
          {flag?.kind === "merge" ? (
            <div data-testid="merge-prompt" className="m-2.5 rounded-sm border border-warning bg-paper-sunk px-3.5 py-2.5">
              <div className="text-xs text-ink-soft">
                Filename matches <b className="text-ink">{flag.merge_with_label}</b>. This looks like a reshoot in a later load.
              </div>
              <div className="mt-2 flex gap-2">
                <Button variant="outline" onClick={() => props.onMerge(s.sample_id, flag.merge_with_sample_id)}>
                  Merge into that sample
                </Button>
                <Button variant="ghost" onClick={() => props.onDismissFlag(s.sample_id)}>Keep separate</Button>
              </div>
            </div>
          ) : null}

          {s.exposures.map((e, i) => (
            <div key={e.id}>
              {flag?.kind === "split" && i === flag.split_at_index ? (
                <div data-testid="split-divider" className="mx-2.5 my-1 flex items-center gap-3 rounded-sm border border-warning bg-paper-sunk px-3 py-1.5">
                  <span className="flex-1 font-mono text-xs text-ink-soft">
                    position jumps <b className="text-warning">{flag.jump_from.toFixed(1)} → {flag.jump_to.toFixed(1)} mm</b> here — likely two samples
                  </span>
                  <Button variant="outline" onClick={() => props.onSplit(s.sample_id)}>Split here</Button>
                </div>
              ) : null}
              <ExposureLeaf exposure={e} thumbSrc={props.thumbSrcFor(e.id)} onMove={(eid) => props.onMoveExposure(s.sample_id, eid)} />
            </div>
          ))}
        </div>
      ) : null}
    </div>
  );
}
```

> **Add `thumbSrcFor: () => null` to `baseProps` in the test** (it's a required prop). The test fixtures already use §8.8 exposure shape (`.id`).

> **Verified primitive APIs:** `Checkbox` takes `checked?: boolean` + `onChange?: (checked: boolean) => void` + `aria-label?` and renders `role="checkbox"` (`Checkbox.tsx:9-23,68`) — `onChange={() => props.onToggleSelect(id)}` is fine (the boolean arg is ignored; the page toggles membership). `Button` (`Button.tsx`) takes `variant?: "solid"|"accent"|"success"|"ghost"|"danger"|"outline"|"ghostInverse"` and extends `ButtonHTMLAttributes` — **there is NO `size` prop** (removed from the JSX above). `IconButton` is not needed here (the Move button lives in `ExposureLeaf`).

> **Design-guard caution:** `border-warning`, `ring-accent`, `text-warning`, `bg-paper-sunk`, `text-headline`, `text-xs` are token/named-scale utilities — permitted. The `border` is a FULL border (not a banned `border-l` side-stripe) — correct per the spec's "full border + a leading word instead of a side-stripe". If you want the dotted-underline edit-in-place name affordance, render the name via `Input variant="title"` in edit mode (page toggles `renaming`) rather than styling a span. Run `node scripts/check-design.mjs`.

- [ ] **Step 4: Run test → PASS + design guard + type-check**
```bash
npm test -- test/SampleFold.test.tsx
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 5: Commit**
```bash
git add src/print/components/SampleFold.tsx test/SampleFold.test.tsx
git commit -m "feat(grouping): SampleFold (name/flag/merge-prompt/split-divider/leaves)"
```

---

## Task 12: `LoadFold` component

The load accordion: chevron · serif title · operator · meta (`N samples · M exposures · time`) · status (`✓ grouped cleanly` / `⚠ N to check`) · the sample folds. Presentational; the page passes the visible samples (post-filter) + a "shown N of M" subset hint + all sample callbacks.

**Files:**
- Create: `src/print/components/LoadFold.tsx`
- Test: `test/LoadFold.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/LoadFold.test.tsx`:
```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { LoadFold } from "../src/print/components/LoadFold";
import type { Load } from "../src/api";

function load(over: Partial<Load> = {}): Load {
  return {
    load_id: 1, load_index: 1, session_id: null, start_time: "10:02", end_time: "10:38", frame_count: 96, note: null,
    samples: [
      { sample_id: 10, name: "A", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
      { sample_id: 20, name: "B", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
        flag: { kind: "merge", merge_with_sample_id: 5, merge_with_label: "C" }, exposures: [] },
    ],
    ...over,
  };
}

const cb = {
  onToggleLoad: () => {}, openSamples: new Set<number>(), selected: new Set<number>(),
  onToggleSampleOpen: () => {}, onToggleSelect: () => {}, onRename: () => {},
  onSplit: () => {}, onMerge: () => {}, onDismissFlag: () => {}, onMoveExposure: () => {},
  thumbSrcFor: () => null,
};

describe("LoadFold", () => {
  it("shows title and an attn status with the flagged count when not clean", () => {
    render(<LoadFold load={load()} open visibleSamples={load().samples} {...cb} />);
    expect(screen.getByText("Load 1")).toBeInTheDocument();
    expect(screen.getByText(/1 to check/)).toBeInTheDocument();
  });
  it("shows a clean status when no sample is flagged", () => {
    const clean = load({ samples: [load().samples[0]!] });
    render(<LoadFold load={clean} open visibleSamples={clean.samples} {...cb} />);
    expect(screen.getByText(/grouped cleanly/i)).toBeInTheDocument();
  });
  it("renders the subset hint when fewer samples are visible than total", () => {
    const l = load();
    render(<LoadFold load={l} open visibleSamples={[l.samples[1]!]} {...cb} />);
    expect(screen.getByText(/1 of 2 shown/)).toBeInTheDocument();
  });
  it("header click toggles the load", () => {
    const onToggleLoad = vi.fn();
    render(<LoadFold load={load()} open={false} visibleSamples={[]} {...cb} onToggleLoad={onToggleLoad} />);
    fireEvent.click(screen.getByRole("button", { name: /load 1/i }));
    expect(onToggleLoad).toHaveBeenCalledWith(1);
  });
});
```

- [ ] **Step 2: Run test → FAIL** (`npm test -- test/LoadFold.test.tsx`).

- [ ] **Step 3: Implement**

Create `src/print/components/LoadFold.tsx`:
```tsx
import type { JSX } from "react";
import type { Load, LoadSample } from "../../api";
import { SampleFold } from "./SampleFold";

export interface LoadFoldProps {
  load: Load;
  open: boolean;
  /** Samples to render (already filtered by the page). */
  visibleSamples: LoadSample[];
  openSamples: Set<number>;
  selected: Set<number>;
  onToggleLoad: (loadId: number) => void;
  onToggleSampleOpen: (sampleId: number) => void;
  onToggleSelect: (sampleId: number) => void;
  onRename: (sampleId: number) => void;
  onSplit: (sampleId: number) => void;
  onMerge: (loserId: number, survivorId: number) => void;
  onDismissFlag: (sampleId: number) => void;
  onMoveExposure: (sampleId: number, exposureId: number) => void;
  thumbSrcFor: (exposureId: number) => string | null;
  className?: string;
}

export function LoadFold(p: LoadFoldProps): JSX.Element {
  const { load: l, open } = p;
  const flaggedCount = l.samples.filter((s) => s.flag).length;
  const clean = flaggedCount === 0;
  const totalExposures = l.samples.reduce((a, s) => a + s.exposures.length, 0);
  const time = l.start_time && l.end_time ? `${l.start_time}–${l.end_time}` : "";
  const hidden = l.samples.length - p.visibleSamples.length;

  return (
    <div data-testid="load-fold" className={`mb-3 overflow-hidden rounded-md border bg-plate${clean ? " border-hair-strong" : " border-warning"}${p.className ? ` ${p.className}` : ""}`}>
      <button
        type="button"
        aria-expanded={open}
        className="flex w-full items-center gap-3 px-4 py-3.5 text-left"
        onClick={() => p.onToggleLoad(l.load_id)}
      >
        <span className="min-w-0">
          <span className="text-headline text-ink">Load {l.load_index}</span>
          {hidden > 0 ? <span className="ml-2 text-xs text-ink-faint">· {p.visibleSamples.length} of {l.samples.length} shown</span> : null}
          <span className="block font-mono text-xs text-ink-faint">{l.samples.length} samples · {totalExposures} exposures{time ? ` · ${time}` : ""}</span>
        </span>
        <span className={`ml-auto text-xs font-semibold${clean ? " text-success" : " text-warning"}`}>
          {clean ? "✓ grouped cleanly" : `⚠ ${flaggedCount} to check`}
        </span>
      </button>

      {open ? (
        <div className="border-t border-hair bg-paper-sunk p-1.5">
          {p.visibleSamples.map((s) => (
            <div key={s.sample_id} className="m-1.5">
              <SampleFold
                sample={s}
                open={p.openSamples.has(s.sample_id) || !!s.flag}
                selected={p.selected.has(s.sample_id)}
                onToggleOpen={p.onToggleSampleOpen}
                onToggleSelect={p.onToggleSelect}
                onRename={p.onRename}
                onSplit={p.onSplit}
                onMerge={p.onMerge}
                onDismissFlag={p.onDismissFlag}
                onMoveExposure={p.onMoveExposure}
                thumbSrcFor={p.thumbSrcFor}
              />
            </div>
          ))}
        </div>
      ) : null}
    </div>
  );
}
```

> `text-success`/`text-warning`/`border-warning`/`border-hair-strong`/`bg-paper-sunk`/`text-headline` are tokens (`--color-success`/`--color-warning`/`--color-hair-strong`/`--color-paper-sunk` all exist in `styles.css`). The header padding uses `px-4` (a real step) — **`px-4.5` is NOT a default Tailwind spacing step**, so the JSX above uses `px-4`. `aria-expanded` satisfies the spec §11 accessibility requirement. Run the guard.

- [ ] **Step 4: Run test → PASS + guard + tsc**
```bash
npm test -- test/LoadFold.test.tsx
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 5: Commit**
```bash
git add src/print/components/LoadFold.tsx test/LoadFold.test.tsx
git commit -m "feat(grouping): LoadFold accordion"
```

---

## Task 13: `GroupingReviewPage` — filter + persistent selection (no edits yet)

The page wires `useLoads`, the `Needs review`/`All loads` filter, the fuzzy/glob `SearchInput`, and a `selected: Set<number>` that PERSISTS across filter changes. Editing actions come in Task 14; bulk-merge in Task 15. This task proves the filter/search/selection-persistence contract.

**Files:**
- Create: `src/print/components/GroupingReviewPage.tsx`
- Create: `src/lib/matchSample.ts` (the glob/fuzzy matcher — pure, separately testable)
- Test: `test/matchSample.test.ts`, `test/GroupingReviewPage.filter.test.tsx`

- [ ] **Step 1: Write the failing matcher test**

Create `test/matchSample.test.ts`:
```ts
import { describe, it, expect } from "vitest";
import { matchSample } from "../src/lib/matchSample";
import type { LoadSample } from "../src/api";

const s: LoadSample = {
  sample_id: 1, name: "HA85 (S01P15)", slot_index: 15,
  grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null,
  exposures: [{ id: 9, filename: "HA_0612_001.tif", horizontal_position: 1, timestamp: null }],
};

describe("matchSample", () => {
  it("empty query matches everything", () => { expect(matchSample(s, "")).toBe(true); });
  it("substring matches name (case-insensitive)", () => { expect(matchSample(s, "ha85")).toBe(true); });
  it("substring matches an exposure filename", () => { expect(matchSample(s, "0612")).toBe(true); });
  it("glob * matches a prefix", () => { expect(matchSample(s, "ha8*")).toBe(true); });
  it("glob ? matches one char", () => { expect(matchSample({ ...s, name: "JC C04" }, "JC C0?")).toBe(true); });
  it("non-matching glob fails", () => { expect(matchSample(s, "ZZ*")).toBe(false); });
});
```

- [ ] **Step 2: Run → FAIL.** `npm test -- test/matchSample.test.ts`

- [ ] **Step 3: Implement the matcher**

Create `src/lib/matchSample.ts` (matches the flat layout of existing `experimentFilter.ts`/`seriesRatio.ts`):
```ts
import type { LoadSample } from "../api";

/** Fuzzy/glob match against a sample's name and exposure filenames.
 *  `*` → any run, `?` → one char (anchored at start); otherwise substring. */
export function matchSample(s: LoadSample, query: string): boolean {
  const q = query.trim().toLowerCase();
  if (!q) return true;
  const fields = [s.name.toLowerCase(), ...s.exposures.map((e) => e.filename.toLowerCase())];
  if (/[*?]/.test(q)) {
    const re = new RegExp(
      "^" + q.replace(/[.+^${}()|[\]\\]/g, "\\$&").replace(/\*/g, ".*").replace(/\?/g, "."),
    );
    return fields.some((f) => re.test(f));
  }
  return fields.some((f) => f.includes(q));
}
```

- [ ] **Step 4: Run matcher test → PASS.**

- [ ] **Step 5: Write the page filter test**

Create `test/GroupingReviewPage.filter.test.tsx`:
```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

const LOADS: Load[] = [
  { load_id: 1, load_index: 1, session_id: null, start_time: "10:02", end_time: "10:38", frame_count: 0, note: null,
    samples: [{ sample_id: 10, name: "HA85 (S01P15)", slot_index: 15, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] }] },
  { load_id: 2, load_index: 2, session_id: null, start_time: "13:10", end_time: "13:51", frame_count: 0, note: null,
    samples: [{ sample_id: 20, name: "HA85 (S02P15)", slot_index: 15, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
      flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "HA85 (S01P15)" }, exposures: [] }] },
];

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return { ...actual, useLoads: () => ({ data: LOADS, isLoading: false }) };
});

function wrap(node: ReactNode) {
  const qc = new QueryClient();
  return render(<QueryClientProvider client={qc}>{node}</QueryClientProvider>);
}

describe("GroupingReviewPage filter + persistent selection", () => {
  it("defaults to 'Needs review' (only the flagged load 2 shows)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    expect(screen.queryByText("Load 1")).toBeNull();
    expect(screen.getByText("Load 2")).toBeInTheDocument();
  });

  it("'All loads' reveals every load", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
    expect(screen.getByText("Load 1")).toBeInTheDocument();
    expect(screen.getByText("Load 2")).toBeInTheDocument();
  });

  it("selection PERSISTS across filter changes (the cross-load merge gesture)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    // show all, select sample 10 in load 1
    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
    fireEvent.click(screen.getByRole("checkbox", { name: /select HA85 \(S01P15\)/i }));
    // search to sample 20's load only — sample 10 leaves the view
    fireEvent.change(screen.getByRole("textbox"), { target: { value: "S02" } });
    // bulk bar still reflects the 1 retained selection
    expect(screen.getByTestId("bulk-bar")).toHaveTextContent("1");
  });
});
```

- [ ] **Step 6: Run page test → FAIL.** `npm test -- test/GroupingReviewPage.filter.test.tsx`

- [ ] **Step 7: Implement the page (filter + selection only)**

First create the local bulk bar `src/print/components/GroupingBulkBar.tsx` (a dark floating action bar — NOT `CullBar`, whose props are hardcoded for corpus-culling: Drop/Keep/Restore/Clear with no parameterized primary label; see FIX rationale below):
```tsx
import type { JSX } from "react";
import { Button } from "../ui/Button";

export interface GroupingBulkBarProps {
  count: number;
  /** singular noun; pluralized with a naive +"s". */
  noun: string;
  primaryLabel: string;
  /** Primary action is disabled below this many selected (merge needs ≥2). */
  primaryEnabled: boolean;
  onPrimary: () => void;
  onClear: () => void;
  className?: string;
}

/** Dark floating bulk-action bar for the grouping selection. Mirrors the
 *  CullBar's dark `bg-ink` idiom but is parameterized for the merge gesture.
 *  Presentational; the page owns the selection + handlers. */
export function GroupingBulkBar(p: GroupingBulkBarProps): JSX.Element {
  return (
    <div
      data-testid="bulk-bar"
      className={`fixed inset-x-0 bottom-6 mx-auto flex w-fit items-center gap-4 rounded-md bg-ink px-5 py-3 text-paper shadow-lg${p.className ? ` ${p.className}` : ""}`}
    >
      <span className="font-mono text-sm">
        {p.count} {p.noun}{p.count === 1 ? "" : "s"} selected
      </span>
      <Button variant="accent" disabled={!p.primaryEnabled} onClick={p.onPrimary}>{p.primaryLabel}</Button>
      <Button variant="ghostInverse" onClick={p.onClear}>Clear</Button>
    </div>
  );
}
```
> **Why a local bar, not `CullBar` (FIX 7 decision):** `CullBar` (`print/components/CullBar.tsx`) is hardcoded to the corpus-culling verbs (`onReject`/`onKeep`/`onRestore`/`onClear` → "Drop"/"Keep"/"Restore"/"Clear"), has `data-testid="cull-bar"` baked in (no passthrough), and has no parameterized primary-label/noun. Parameterizing it to also serve a single-primary "Merge" bar would bloat its contract and risk the SheetTable caller. A small local `GroupingBulkBar` is the cleaner reuse-the-primitives (`Button`, the `bg-ink` dark idiom) path. `bg-ink`/`text-paper` are tokens; `ghostInverse` is the real dark-surface Button variant.

Now create `src/print/components/GroupingReviewPage.tsx`:
```tsx
import { useMemo, useState, type JSX } from "react";
import type { Load } from "../../api";
import { useLoads } from "../../queries";
import { LoadFold } from "./LoadFold";
import { GroupingBulkBar } from "./GroupingBulkBar";
import { SearchInput } from "../ui/SearchInput";
import { SegmentedControl } from "../ui/SegmentedControl";
import { EmptyState } from "../ui/EmptyState";
import { Kicker } from "../ui/Kicker";
import { matchSample } from "../../lib/matchSample";

type Filter = "attn" | "all";

export interface GroupingReviewPageProps {
  experimentId: number;
  onBack: () => void;
  className?: string;
}

export function GroupingReviewPage({ experimentId, onBack, className }: GroupingReviewPageProps): JSX.Element {
  const { data: loads = [], isLoading } = useLoads(experimentId);
  const [filter, setFilter] = useState<Filter>("attn");
  const [search, setSearch] = useState("");
  // ORDERED selection (first-selected = bulk-merge survivor — Task 15). Membership
  // checks use `.includes`; the LoadFold/SampleFold `selected` prop takes a Set,
  // so derive one. Selection PERSISTS across filter changes (never cleared here).
  const [selection, setSelection] = useState<number[]>([]);
  const selectedSet = useMemo(() => new Set(selection), [selection]);
  const [openLoads, setOpenLoads] = useState<Set<number>>(new Set());
  const [openSamples, setOpenSamples] = useState<Set<number>>(new Set());

  const q = search.trim();
  const flaggedTotal = useMemo(
    () => loads.reduce((a, l) => a + l.samples.filter((s) => s.flag).length, 0),
    [loads],
  );
  const totalSamples = useMemo(() => loads.reduce((a, l) => a + l.samples.length, 0), [loads]);

  // Visible (load, samples) pairs. A searched query overrides the attn filter
  // (search across all loads); a selected sample stays visible even when it
  // would be filtered out, so the persistent selection is verifiable.
  const visible = useMemo(() => {
    const base: Load[] = !q && filter === "attn" ? loads.filter((l) => l.samples.some((s) => s.flag)) : loads;
    return base
      .map((l) => ({
        load: l,
        samples: q
          ? l.samples.filter((s) => matchSample(s, q) || selectedSet.has(s.sample_id))
          : l.samples,
      }))
      .filter((x) => x.samples.length > 0);
  }, [loads, filter, q, selectedSet]);

  const toggleSelect = (id: number) =>
    setSelection((prev) => (prev.includes(id) ? prev.filter((x) => x !== id) : [...prev, id]));

  const toggleSet = (id: number, set: (u: (p: Set<number>) => Set<number>) => void) =>
    set((p) => { const n = new Set(p); n.has(id) ? n.delete(id) : n.add(id); return n; });

  const noop = () => {};

  return (
    <div className={`mx-auto max-w-[1180px] px-10 pb-32 pt-8${className ? ` ${className}` : ""}`}>
      <button type="button" className="mb-4 inline-flex items-center gap-1.5 text-xs font-semibold text-accent" onClick={onBack}>
        ← Back to corpus
      </button>
      <div className="flex items-end justify-between gap-8">
        <div>
          <Kicker>Grouping review</Kicker>
          <h1 className="text-display text-ink">Check the grouping</h1>
          <p className="mt-2 max-w-[66ch] text-sm text-ink-soft">
            Confirm every sample loaded, has all its exposures, and is split where it should be.
          </p>
        </div>
        <div className="shrink-0 text-right">
          <div className="text-display text-ink"><b className="text-accent">{flaggedTotal}</b> flagged</div>
          <div className="text-xs font-bold uppercase text-ink-faint">of {totalSamples} samples</div>
        </div>
      </div>

      <div className="my-4 flex items-center gap-2">
        <SegmentedControl<Filter>
          value={filter}
          onChange={setFilter}
          options={[{ value: "attn", label: "Needs review" }, { value: "all", label: "All loads" }]}
          aria-label="Grouping filter"
        />
        <div className="ml-auto w-80">
          <SearchInput
            value={search}
            onChange={setSearch}
            ariaLabel="Filter samples"
            placeholder="Filter samples — name or glob, e.g. HA8* or JC C0?"
          />
        </div>
      </div>

      {isLoading ? null : visible.length === 0 ? (
        <EmptyState title={`No samples match “${q}”`} />
      ) : (
        visible.map(({ load, samples }) => (
          <LoadFold
            key={load.load_id}
            load={load}
            open={openLoads.has(load.load_id) || !!q || load.samples.some((s) => s.flag)}
            visibleSamples={samples}
            openSamples={openSamples}
            selected={selectedSet}
            onToggleLoad={(id) => toggleSet(id, setOpenLoads)}
            onToggleSampleOpen={(id) => toggleSet(id, setOpenSamples)}
            onToggleSelect={toggleSelect}
            onRename={noop}
            onSplit={noop}
            onMerge={noop}
            onDismissFlag={noop}
            onMoveExposure={noop}
            thumbSrcFor={() => null /* Task 14 wires a real detector-URL resolver */}
          />
        ))
      )}

      {selection.length > 0 ? (
        <GroupingBulkBar
          count={selection.length}
          noun="sample"
          primaryLabel="Merge"
          primaryEnabled={selection.length >= 2}
          onPrimary={noop /* Task 15 */}
          onClear={() => setSelection([])}
        />
      ) : null}
    </div>
  );
}
```

> **Verified primitive APIs:** `SegmentedControl<T>` takes `value`/`onChange: (v: T) => void`/`options: SegmentOption<T>[]` + (required) `aria-label` — pass the generic `<Filter>` so `onChange={setFilter}` types cleanly (no `as Filter` cast). `SearchInput` takes `value` + **`onChange: (v: string) => void`** (NOT `onValueChange`) + optional **`ariaLabel`** + `placeholder` — it renders a `role="textbox"` input (the filter test's `getByRole("textbox")`). `EmptyState`/`Kicker` are exported primitives (confirm `EmptyState`'s `title` prop name against `print/ui/EmptyState.tsx`; adjust if it's `heading`/`message`). The bulk bar carries `data-testid="bulk-bar"` (the filter test asserts its text contains the count).

> **`max-w-[1180px]`/`w-80`/`pt-8`/`px-10` are layout utilities** (width/padding) — permitted outside `print/ui/`. `text-display`/`text-accent`/`text-ink*` are named roles/tokens. Run `node scripts/check-design.mjs`.

- [ ] **Step 8: Run page test → PASS + guard + tsc**
```bash
npm test -- test/GroupingReviewPage.filter.test.tsx test/matchSample.test.ts
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 9: Commit**
```bash
git add src/print/components/GroupingReviewPage.tsx src/print/components/GroupingBulkBar.tsx src/lib/matchSample.ts test/GroupingReviewPage.filter.test.tsx test/matchSample.test.ts
git commit -m "feat(grouping): GroupingReviewPage filter + persistent ordered selection + glob matcher + bulk bar"
```

---

## Task 14: Wire the single-entity edits (rename / move / split / merge-prompt / dismiss) into the page

Replace the `noop`s with the queue hooks + guarded confirm (`ModalShell`) + undo toast. Rename uses inline `Input variant="title"`; Split/Merge open a `ModalShell` preview-and-confirm; the merge-prompt's "Merge into that sample" and "Keep separate" wire to `useMergeSamples` / a local dismiss; per-exposure Move opens a `Menu` of same-load targets.

**Files:**
- Modify: `src/print/components/GroupingReviewPage.tsx`
- Test: `test/GroupingReviewPage.edits.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/GroupingReviewPage.edits.test.tsx`:
```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

const mergeMutate = vi.fn();
const renameMutate = vi.fn();

const LOADS: Load[] = [
  { load_id: 2, load_index: 2, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [{ sample_id: 20, name: "HA85 (S02P15)", slot_index: 15, grouping_source: "auto_position", name_source: "auto", merged_into_id: null,
      flag: { kind: "merge", merge_with_sample_id: 10, merge_with_label: "HA85 (S01P15)" }, exposures: [] }] },
];

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  // Hooks now take only experimentId (settled in Task 9 — entity ids go in mutate input).
  return {
    ...actual,
    useLoads: () => ({ data: LOADS, isLoading: false }),
    useMergeSamples: (_experimentId: number) => ({ mutate: mergeMutate, isPending: false }),
    useRenameSample: (_experimentId: number) => ({ mutate: renameMutate, isPending: false }),
    useMoveExposure: (_experimentId: number) => ({ mutate: vi.fn(), isPending: false }),
    useSplitSample: (_experimentId: number) => ({ mutate: vi.fn(), isPending: false }),
    useDismissGroupingFlag: (_experimentId: number) => ({ mutate: vi.fn(), isPending: false }),
  };
});

function wrap(node: ReactNode) {
  return render(<QueryClientProvider client={new QueryClient()}>{node}</QueryClientProvider>);
}

beforeEach(() => { mergeMutate.mockClear(); renameMutate.mockClear(); });

describe("GroupingReviewPage edits", () => {
  it("clicking the merge-prompt's Merge opens a confirm, then fires useMergeSamples with {loserId,survivorId}", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: /merge into that sample/i }));
    // confirm modal
    fireEvent.click(screen.getByRole("button", { name: /^merge$/i }));
    // assert on the FIRST arg only (the page may or may not pass a 2nd options arg)
    expect(mergeMutate).toHaveBeenCalledTimes(1);
    expect(mergeMutate.mock.calls[0]![0]).toEqual({ loserId: 20, survivorId: 10 });
  });
});
```

- [ ] **Step 2: Run → FAIL.** `npm test -- test/GroupingReviewPage.edits.test.tsx`

- [ ] **Step 3: Wire the hooks + confirm + toast**

In `GroupingReviewPage.tsx`:
- Import `useMergeSamples`, `useRenameSample`, `useMoveExposure`, `useSplitSample`, `useDismissGroupingFlag`, `ModalShell`, `Button`, `Menu`, `Input`, `showToast` (from `lib/toast`), and `useUndoStack`.
- Instantiate all five hooks scoped to `experimentId` only (SETTLED in Task 9 — row-scoped hooks take ONLY `experimentId`; entity ids travel in the `mutate()` input, so the page holds one instance of each hook and passes `sampleId`/`exposureId` at call time): `const renameMutate = useRenameSample(experimentId).mutate`, `const moveMutate = useMoveExposure(experimentId).mutate`, etc.
- Add local state: `confirm: { kind: "merge"|"split"; … } | null`, `renamingId: number | null`, `movePicker: { sampleId; exposureId; anchorEl } | null`, and `const undoStack = useUndoStack<{ label: string; undo: () => void }>()`.
- Replace `onMerge` noop: open a merge confirm; on confirm, `mergeMutate({ loserId, survivorId })` and `showToast(\`Merged into …\`, "info")` WITHOUT an undo action (merge is multi-row — one stamp can't reverse it server-side; spec §9.3 "undo is session-local" and v1 omits the merge undo action; flag to Jonathan in Self-Review).
- Replace `onRename` noop: set `renamingId`; render the sample name as `Input variant="title"` in `SampleFold` when `renamingId === sample_id` (thread `renaming`/`onCommitRename` props into `SampleFold`); on commit, capture `prev = sample.name`, call `renameMutate({ sampleId, name })`, then `undoStack.push({ label: "rename", undo: () => renameMutate({ sampleId, name: prev }) })` and `showToast("Renamed", "info", { label: "Undo", onClick: () => { const e = undoStack.pop(); if (e) e.undo(); } })`.
- Replace `onMoveExposure` noop: open a `Menu` anchored at the row listing same-load sibling samples; on pick, capture `fromSampleId`, call `moveMutate({ exposureId, sampleId: destId })`, push an undo that moves it back (`moveMutate({ exposureId, sampleId: fromSampleId })`), and show the toast with the Undo action.
- Replace `onSplit` noop: open a split confirm previewing the two-way split at `flag.split_at_index` (or the midpoint for a manual split); on confirm, `splitMutate({ sampleId, exposureIds, name })` + a no-action toast (multi-row).
- Replace `onDismissFlag` noop: **a DURABLE mutation** — `dismissMutate({ sampleId, flagKind: flag.kind, mergeWithSampleId: flag.kind === "merge" ? flag.merge_with_sample_id : undefined })` (spec §9.3: "Keep separate" writes `grouping_flag_dismissed`, suppressed in the roll-up across rescans — **NOT** a session-local `Set`). The optimistic mutator clears `sample.flag`, so the sample drops from "Needs review" immediately; the durable event keeps it gone. Show a toast with an Undo action that re-issues an inverse only if Phase D exposes an un-dismiss path; otherwise omit the action (the dismiss carries `undoes_event_id` server-side per §9.3, so a future un-dismiss is possible — for v1, omit the action and note it).

> **The "Needs review"/flagged counts derive from `sample.flag != null`** (the backend-produced flag), NOT a client-side `deriveFlag`. There is NO client-side flag derivation anywhere in this plan — `sample.flag` is read straight off the roll-up (spec §8.8/§9.3: the frontend is a pure consumer; the grouper computes the flag). Dismiss flips `flag` to null durably; it never recomputes it.

> Keep the confirm modal primary button **outline-accent** (`Button variant="outline"`) so Cancel is co-equal (spec §11: "Demote primary actions from filled-accent to outline-accent"). The live Button has no `out-accent` variant — use `variant="outline"`.

> **`thumbSrcFor`:** wire a real detector-thumbnail URL resolver here (the corpus/loupe already build detector image URLs from an exposure id — reuse that helper) and pass it down through `LoadFold` → `SampleFold` → `ExposureLeaf`. If no shared helper exists yet, `() => null` is acceptable for v1 (Thumbnail draws a placeholder); note it.

> **Undo via `useUndoStack` + the imperative toast:** for rename/move, `undoStack.push({ label, undo })` then `showToast(msg, "info", { label: "Undo", onClick: () => { const e = undoStack.pop(); if (e) e.undo(); } })`. The toast action slot is Task 8; `pop()` is the StrictMode-safe API (Task 1). **Verify under `<StrictMode>`** in this task's test that one rename → one history entry → one inverse on undo.

- [ ] **Step 4: Run → PASS + guard + tsc**
```bash
npm test -- test/GroupingReviewPage.edits.test.tsx test/GroupingReviewPage.filter.test.tsx
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 5: Commit**
```bash
git add src/print/components/GroupingReviewPage.tsx src/print/components/SampleFold.tsx test/GroupingReviewPage.edits.test.tsx
git commit -m "feat(grouping): wire rename/move/split/merge-prompt/dismiss with confirm + undo"
```

---

## Task 15: Bulk-merge of the persistent selection

The bulk bar's "Merge" (enabled at ≥2 selected) merges all selected into the first-selected survivor: issue one `mergeMutate` per loser. Confirm in a `ModalShell` previewing `N → survivor`.

**Files:**
- Modify: `src/print/components/GroupingReviewPage.tsx`
- Test: `test/GroupingReviewPage.bulk.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/GroupingReviewPage.bulk.test.tsx`:
```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, within } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import type { Load } from "../src/api";
import { GroupingReviewPage } from "../src/print/components/GroupingReviewPage";

const mergeMutate = vi.fn();
const LOADS: Load[] = [
  { load_id: 1, load_index: 1, session_id: null, start_time: null, end_time: null, frame_count: 0, note: null,
    samples: [
      { sample_id: 10, name: "A", slot_index: 1, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
      { sample_id: 11, name: "B", slot_index: 2, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
      { sample_id: 12, name: "C", slot_index: 3, grouping_source: "auto_position", name_source: "auto", merged_into_id: null, flag: null, exposures: [] },
    ] },
];
vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  // Hooks take only experimentId (settled in Task 9 — entity ids go in mutate input).
  return {
    ...actual, useLoads: () => ({ data: LOADS, isLoading: false }),
    useMergeSamples: (_experimentId: number) => ({ mutate: mergeMutate, isPending: false }),
    useRenameSample: (_experimentId: number) => ({ mutate: vi.fn() }),
    useMoveExposure: (_experimentId: number) => ({ mutate: vi.fn() }),
    useSplitSample: (_experimentId: number) => ({ mutate: vi.fn() }),
    useDismissGroupingFlag: (_experimentId: number) => ({ mutate: vi.fn() }),
  };
});
const wrap = (n: ReactNode) => render(<QueryClientProvider client={new QueryClient()}>{n}</QueryClientProvider>);
beforeEach(() => mergeMutate.mockClear());

describe("bulk merge", () => {
  it("merges all selected losers into the first-selected survivor (one mutate per loser)", () => {
    wrap(<GroupingReviewPage experimentId={7} onBack={() => {}} />);
    fireEvent.click(screen.getByRole("button", { name: /all loads/i }));
    fireEvent.click(screen.getByRole("checkbox", { name: /select A/i }));   // survivor
    fireEvent.click(screen.getByRole("checkbox", { name: /select B/i }));
    fireEvent.click(screen.getByRole("checkbox", { name: /select C/i }));
    // bulk-bar primary (scope the query to the bar so it doesn't collide with
    // the later modal confirm button, which also reads "Merge")
    const bar = screen.getByTestId("bulk-bar");
    fireEvent.click(within(bar).getByRole("button", { name: /^merge$/i }));
    fireEvent.click(screen.getByRole("button", { name: /^merge$/i }));        // confirm modal
    expect(mergeMutate).toHaveBeenCalledTimes(2);
    // assert first-arg inputs (the page may pass a 2nd options arg)
    const inputs = mergeMutate.mock.calls.map((c) => c[0]);
    expect(inputs).toContainEqual({ loserId: 11, survivorId: 10 });
    expect(inputs).toContainEqual({ loserId: 12, survivorId: 10 });
  });
});
```

> **Selection order:** to make "first-selected = survivor" deterministic, track selection as an ORDERED list (an array, or a `Map<number, order>`), not a bare `Set` — but the persistent-selection contract (Task 13) is about membership, not order. Refine `selected` to an array (or keep the Set for membership + a separate `selectionOrder: number[]`). The test asserts survivor=10 (A, first clicked). Use whichever keeps Task 13's tests green; an ordered array satisfies both.

- [ ] **Step 2: Run → FAIL.** `npm test -- test/GroupingReviewPage.bulk.test.tsx`

- [ ] **Step 3: Implement bulk-merge**

In `GroupingReviewPage.tsx`, give the bulk bar `onPrimary={() => openBulkMergeConfirm()}`; in the confirm's onConfirm, compute `[survivor, ...losers] = selectionOrder`, then `losers.forEach((loser) => mergeMutate({ loserId: loser, survivorId: survivor }, {}))`, clear the selection, show a toast (`Merged N samples into <survivor>`). Disable the bulk "Merge" when `< 2` selected (the bar can show but Merge is disabled; matching the prototype's `bulkMerge` early-return).

- [ ] **Step 4: Run → PASS + guard + tsc + re-run Tasks 13/14 tests**
```bash
npm test -- test/GroupingReviewPage.bulk.test.tsx test/GroupingReviewPage.filter.test.tsx test/GroupingReviewPage.edits.test.tsx
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 5: Commit**
```bash
git add src/print/components/GroupingReviewPage.tsx test/GroupingReviewPage.bulk.test.tsx
git commit -m "feat(grouping): bulk-merge of the persistent selection"
```

---

## Task 16: `GeometryLedger` component

Per-field rows: label · value (mono, edit-in-place via `Input mono`) · provenance chip (`PRP`/`setup files`/`edited`) · Override/Edit + Revert (when user-sourced); a "Undo last change" affordance in the card header; a multi-setup discrepancy banner when triggered. Presentational; `ConfigurationBody` (Task 19) owns the override/undo state.

**Files:**
- Create: `src/print/components/GeometryLedger.tsx`
- Test: `test/GeometryLedger.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/GeometryLedger.test.tsx`:
```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { GeometryLedger, type GeometryRow } from "../src/print/components/GeometryLedger";

const ROWS: GeometryRow[] = [
  { key: "beam_energy", label: "Beam energy", value: "9.00 keV", source: "prp" },
  { key: "beam_center_x", label: "Beam center X", value: "421.4 px", source: "setup" },
  { key: "flight_path", label: "Flight path", value: "1.81 m", source: "user" },
];

const cb = { onOverride: () => {}, onRevert: () => {}, onUndo: () => {}, canUndo: false };

describe("GeometryLedger", () => {
  it("renders each row with a provenance chip", () => {
    render(<GeometryLedger rows={ROWS} {...cb} />);
    expect(screen.getByText("Beam energy")).toBeInTheDocument();
    expect(screen.getByText("PRP")).toBeInTheDocument();
    expect(screen.getByText("setup files")).toBeInTheDocument();
    expect(screen.getByText("edited")).toBeInTheDocument();
  });
  it("user rows expose Revert; calls onRevert with the key", () => {
    const onRevert = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onRevert={onRevert} />);
    fireEvent.click(screen.getByRole("button", { name: /revert/i }));
    expect(onRevert).toHaveBeenCalledWith("flight_path");
  });
  it("Override calls onOverride with the key", () => {
    const onOverride = vi.fn();
    render(<GeometryLedger rows={ROWS} {...cb} onOverride={onOverride} />);
    fireEvent.click(screen.getAllByRole("button", { name: /override/i })[0]!);
    expect(onOverride).toHaveBeenCalledWith("beam_energy");
  });
  it("shows the discrepancy banner when discrepancies > 0", () => {
    render(<GeometryLedger rows={ROWS} {...cb} discrepancyCount={2} />);
    expect(screen.getByText(/geometry check found 2 issues/i)).toBeInTheDocument();
  });
  it("Undo last change is gated on canUndo", () => {
    const onUndo = vi.fn();
    const { rerender } = render(<GeometryLedger rows={ROWS} {...cb} canUndo={false} onUndo={onUndo} />);
    expect(screen.queryByRole("button", { name: /undo last change/i })).toBeNull();
    rerender(<GeometryLedger rows={ROWS} {...cb} canUndo onUndo={onUndo} />);
    fireEvent.click(screen.getByRole("button", { name: /undo last change/i }));
    expect(onUndo).toHaveBeenCalled();
  });
});
```

- [ ] **Step 2: Run → FAIL.** `npm test -- test/GeometryLedger.test.tsx`

- [ ] **Step 3: Implement**

Create `src/print/components/GeometryLedger.tsx`:
```tsx
import type { JSX } from "react";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";

export type GeometrySource = "prp" | "setup" | "user" | "default";

export interface GeometryRow {
  key: string;
  label: string;
  value: string;
  source: GeometrySource;
}

export interface GeometryLedgerProps {
  rows: GeometryRow[];
  onOverride: (key: string) => void;
  onRevert: (key: string) => void;
  onUndo: () => void;
  canUndo: boolean;
  discrepancyCount?: number;
  className?: string;
}

const SOURCE_LABEL: Record<GeometrySource, string> = {
  prp: "PRP", setup: "setup files", user: "edited", default: "unset",
};

export function GeometryLedger(p: GeometryLedgerProps): JSX.Element {
  return (
    <Card className={p.className}>
      <div className="flex items-center justify-between border-b border-hair px-4 py-3.5">
        <h3 className="text-headline text-ink">Geometry</h3>
        <div className="flex items-center gap-3.5">
          {p.canUndo ? (
            <button type="button" className="text-xs font-semibold text-accent" onClick={p.onUndo}>↺ Undo last change</button>
          ) : null}
          <span className="text-xs text-ink-faint">auto-derived</span>
        </div>
      </div>

      {p.discrepancyCount && p.discrepancyCount > 0 ? (
        <div className="m-4 rounded-sm border border-warning bg-paper-sunk px-3.5 py-2.5 text-xs text-ink-soft">
          ⚠ Geometry check found {p.discrepancyCount} issue{p.discrepancyCount > 1 ? "s" : ""} — constant fields varied across PRPs.
        </div>
      ) : null}

      <div className="px-4 pb-3.5">
        {p.rows.map((r) => (
          <div key={r.key} className="flex items-center gap-3 border-b border-hair py-2.5 last:border-b-0">
            <span className="w-32 shrink-0 text-xs text-ink-soft">{r.label}</span>
            <span className="font-mono text-sm font-medium text-ink">{r.value}</span>
            <SourceChip source={r.source} />
            {r.source === "user" ? (
              <Button variant="ghost" aria-label={`Revert ${r.label}`} onClick={() => p.onRevert(r.key)}>↺ Revert</Button>
            ) : null}
            <Button variant="ghost" aria-label={`Override ${r.label}`} onClick={() => p.onOverride(r.key)}>
              {r.source === "user" ? "Edit" : "Override"}
            </Button>
          </div>
        ))}
      </div>
    </Card>
  );
}
```

> **`Card`/`Button` props (verified):** `Card` accepts `className`/`children` (+ `as`/`elevated`/`padding`/etc., all optional) — `<Card className={…}>…</Card>` is valid. `Button` has **no `size` prop** (removed `size="sm"` above) and `variant="ghost"` exists. The chip's accent-wash is extracted into a primitive (`SourceChip`, below) so this design-guard-NON-exempt file carries no `bg-[color-mix(...)]`.

> **`SourceChip` — the accent-wash chip lives in a primitive or behind a real token (FIX 8).** `GeometryLedger` is in `print/components/` (NOT design-guard-exempt), and there is **no `--color-accent-wash` token** today (`styles.css` has `--color-accent` but not `--color-accent-wash`). Pick ONE:
> - **(a) Add the token** to `styles.css` `@theme`: `--color-accent-wash: oklch(from var(--color-accent) l c h / 0.12);` (or a flat oklch), then use `bg-accent-wash text-accent` on the chip. Tailwind 4 generates `bg-accent-wash`/`text-accent-wash` from the `@theme` token. **This file may then use `bg-accent-wash`** because it is a named token utility, not an arbitrary `bg-[…]`.
> - **(b) Render the chip via an existing accent-tinted primitive** — `Chip` (`print/ui/Chip.tsx`, exported, has `variant`s) or `Badge`/`NoticePill`. If one supports an accent tone, use it and drop the local chip entirely.
>
> The plan assumes (a). Add a tiny local presentational `SourceChip` ONLY if it lives in `print/ui/` (exempt); otherwise inline the token utilities:
> ```tsx
> // inline form (after adding --color-accent-wash):
> function SourceChip({ source }: { source: GeometrySource }) {
>   const user = source === "user";
>   return (
>     <span className={`ml-auto rounded-sm px-1.5 py-0.5 text-xs font-bold uppercase${user ? " bg-accent-wash text-accent" : " bg-paper-sunk text-ink-faint"}`}>
>       {SOURCE_LABEL[source]}
>     </span>
>   );
> }
> ```
> If you choose (b), delete `SourceChip` and render the primitive inline. Confirm with `node scripts/check-design.mjs` — it must pass with NO `bg-[…]`/`color-mix` literal in this file.

- [ ] **Step 4: Run → PASS + guard + tsc**
```bash
npm test -- test/GeometryLedger.test.tsx
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 5: Commit**
```bash
git add src/print/components/GeometryLedger.tsx src/styles.css test/GeometryLedger.test.tsx
git commit -m "feat(config): GeometryLedger with provenance chips + override/revert/undo + discrepancy banner"
```
(Stage `src/styles.css` only if you chose option (a) and added `--color-accent-wash`.)

---

## Task 17: `SourcesCard` component

Renders the directory + pattern rows. The 3 pattern fields (`image_pattern`, `metadata_pattern`, `integration_pattern`) are NOW DEFINITE typed columns on `experiments` (E1 adds them via additive migration — DECISION); render them as EDITABLE (edit-in-place, triggers a rescan — backend handles rescan on pattern write; the frontend just PATCHes via `ExperimentPatch`). `data_dir` and `analysis_dir` are READ-ONLY (DECISION: the directory is fixed at create; do NOT add an edit affordance for them). Plus the "edits apply on the next rescan" hint and a Rescan-now affordance slot. Presentational; the page owns the values + onEdit.

**Files:**
- Create: `src/print/components/SourcesCard.tsx`
- Test: `test/SourcesCard.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/SourcesCard.test.tsx`:
```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { SourcesCard, type SourceRow } from "../src/print/components/SourcesCard";

const ROWS: SourceRow[] = [
  { key: "data_dir", label: "Data directory", value: "/Volumes/data/ssrl/2026_04/1p7m" },
  { key: "image_pattern", label: "Image pattern", value: "{name}.tiff" },
];

describe("SourcesCard", () => {
  it("renders each source row", () => {
    render(<SourcesCard rows={ROWS} onEdit={() => {}} onRescan={() => {}} />);
    expect(screen.getByText("Data directory")).toBeInTheDocument();
    expect(screen.getByText("{name}.tiff")).toBeInTheDocument();
  });
  it("editing a value calls onEdit(key, newValue)", () => {
    const onEdit = vi.fn();
    render(<SourcesCard rows={ROWS} onEdit={onEdit} onRescan={() => {}} />);
    fireEvent.click(screen.getByText("{name}.tiff"));
    const input = screen.getByDisplayValue("{name}.tiff");
    fireEvent.change(input, { target: { value: "{name}.tif" } });
    fireEvent.keyDown(input, { key: "Enter" });
    expect(onEdit).toHaveBeenCalledWith("image_pattern", "{name}.tif");
  });
  it("Rescan now calls onRescan", () => {
    const onRescan = vi.fn();
    render(<SourcesCard rows={ROWS} onEdit={() => {}} onRescan={onRescan} />);
    fireEvent.click(screen.getByRole("button", { name: /rescan now/i }));
    expect(onRescan).toHaveBeenCalled();
  });
});
```

- [ ] **Step 2: Run → FAIL.** `npm test -- test/SourcesCard.test.tsx`

- [ ] **Step 3: Implement**

Create `src/print/components/SourcesCard.tsx` — a `Card` with rows; each value is a span that becomes an `Input` (mono) on click and commits onEnter/onBlur calling `onEdit(key, value)`; a footer hint + a `Button` "Rescan now" → `onRescan`. Use `useState` for the per-row editing id. Reuse `Input` (mono variant), `Card`, `Button`. Appearance via tokens only.

```tsx
import { useState, type JSX } from "react";
import { Card } from "../ui/Card";
import { Button } from "../ui/Button";
import { Input } from "../ui/Input";

export interface SourceRow { key: string; label: string; value: string }

export interface SourcesCardProps {
  rows: SourceRow[];
  onEdit: (key: string, value: string) => void;
  onRescan: () => void;
  className?: string;
}

export function SourcesCard(p: SourcesCardProps): JSX.Element {
  const [editing, setEditing] = useState<string | null>(null);
  const [draft, setDraft] = useState("");
  const begin = (r: SourceRow) => { setEditing(r.key); setDraft(r.value); };
  const commit = () => { if (editing) p.onEdit(editing, draft.trim()); setEditing(null); };

  return (
    <Card className={p.className}>
      <div className="flex items-center justify-between border-b border-hair px-4 py-3.5">
        <h3 className="text-headline text-ink">Sources</h3>
        <span className="text-xs text-ink-faint">edits apply on the next rescan</span>
      </div>
      <div className="px-4 pb-2">
        {p.rows.map((r) => (
          <div key={r.key} className="flex items-center gap-3 border-b border-hair py-2.5 last:border-b-0">
            <span className="w-40 shrink-0 text-xs text-ink-soft">{r.label}</span>
            {editing === r.key ? (
              <Input
                mono autoFocus value={draft} onValueChange={setDraft} aria-label={r.label}
                onKeyDown={(e) => { if (e.key === "Enter") commit(); if (e.key === "Escape") setEditing(null); }}
                onBlur={commit}
              />
            ) : (
              <button type="button" className="font-mono text-sm text-ink" onClick={() => begin(r)}>{r.value}</button>
            )}
          </div>
        ))}
      </div>
      <div className="flex items-center justify-between gap-4 border-t border-hair px-4 py-3.5">
        <span className="max-w-[64ch] text-xs text-ink-soft">
          Rescan checks for new files and ingests only those. It never re-reads existing exposures or overwrites your edits.
        </span>
        <Button variant="outline" onClick={p.onRescan}>Rescan now</Button>
      </div>
    </Card>
  );
}
```

> **Verified:** `Input` requires `value` + `onValueChange` and extends `InputHTMLAttributes` (omitting `size`), so `onKeyDown`/`onBlur`/`autoFocus` pass through; it has `mono`. `Button` has **no `size` prop** (removed `size="sm"`). `px-4.5` replaced with `px-4` (a real step). The edit-in-place value is a `<button>` that swaps to `Input` — appearance stays in the primitive. Run the guard.
>
> **DECISION: which rows are editable vs read-only.** `data_dir` and `analysis_dir` pass `editable={false}` (or are rendered as plain `<span>`s with no click-to-edit). The 3 pattern rows (`image_pattern`, `metadata_pattern`, `integration_pattern`) are editable (the DECISION adds them as typed `experiments` columns; editing a pattern writes via `ExperimentPatch` and the backend triggers a rescan). Do NOT render an edit affordance for `data_dir`/`analysis_dir` — the directory is fixed at create. If E1 has not yet landed the pattern columns, the implementer should STOP and confirm with the E1 owner before proceeding.

- [ ] **Step 4: Run → PASS + guard + tsc**
```bash
npm test -- test/SourcesCard.test.tsx
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 5: Commit**
```bash
git add src/print/components/SourcesCard.tsx test/SourcesCard.test.tsx
git commit -m "feat(config): SourcesCard editable source rows + rescan affordance"
```

---

## Task 18: `AcquisitionTimeline` — exposure-count bars on the trace-plot scaffolding

Spec §8.7 / E2 mandate: reuse the focus/series d3 declarative trace-plot infrastructure (`print/plot/`), NOT a new chart engine. The acquisition timeline is a per-session cluster of per-load bars (exposure-count height). Build it on the existing `PlotFrame`/`Axis` scaffolding (the plot engine's reusable axis/frame primitives) rather than authoring a standalone SVG chart.

**Files:**
- Create: `src/print/components/AcquisitionTimeline.tsx`
- Test: `test/AcquisitionTimeline.test.tsx`

- [ ] **Step 0: Plot scaffolding API (RESOLVED)**

`print/plot/` is a **trace-only declarative d3** engine — it exports `TracePlot`, `TraceLine`, `PlotPeaks` (line/mark renderers) plus **reusable** primitives `PlotFrame` (`print/plot/PlotFrame.tsx` — a measured `<svg>` container with a `render: (dims: PlotDims) => ReactNode` slot, types `Margins`/`PlotDims`), `Axis` (`print/plot/Axis.tsx` — `orientation: "bottom" | "left"`), and scale helpers `makeProjection`/`makeAxis`/`positiveExtent` (`print/plot/projection.ts`, `ScaleType` linear/log/sqrt). There is **NO reusable bar renderer** and **no `cleanFigureSvg` in `print/plot/`** (figure export lives in `lib/figure-export/`). 

**Resolution (FIX 12):** the acquisition bars are a thin hand-authored SVG in **`print/plot/AcquisitionChart.tsx`** — `print/plot/**` IS design-guard-exempt, so it may author `fill`/`stroke` appearance directly. It does NOT need the `PlotFrame`/`Axis` line-trace machinery (those assume a continuous trace, not a categorical bar cluster); a small linear y-scale + a categorical x-layout is simpler and correct. Reuse only the engine's *token language* (the axis hairline/label colours) so it reads as part of the plot family. **Do not** introduce Observable Plot or a second chart lib (the engine is the sole plot layer — `project_greenfield_trace_plot_engine`).

- [ ] **Step 1: Write the failing test**

Create `test/AcquisitionTimeline.test.tsx`:
```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { AcquisitionTimeline, type AcqSession } from "../src/print/components/AcquisitionTimeline";

const SESSIONS: AcqSession[] = [
  { label: "Apr 12", loadFrameCounts: [42, 30, 55] },
  { label: "Apr 13", loadFrameCounts: [48, 52] },
];

describe("AcquisitionTimeline", () => {
  it("renders one bar per load across all sessions", () => {
    render(<AcquisitionTimeline sessions={SESSIONS} />);
    expect(screen.getAllByTestId("acq-bar")).toHaveLength(5);
  });
  it("renders a session label per session", () => {
    render(<AcquisitionTimeline sessions={SESSIONS} />);
    expect(screen.getByText(/APR 12/i)).toBeInTheDocument();
    expect(screen.getByText(/APR 13/i)).toBeInTheDocument();
  });
  it("renders nothing meaningful (empty frame) for no sessions", () => {
    render(<AcquisitionTimeline sessions={[]} />);
    expect(screen.queryAllByTestId("acq-bar")).toHaveLength(0);
  });
});
```

- [ ] **Step 2: Run → FAIL.** `npm test -- test/AcquisitionTimeline.test.tsx`

- [ ] **Step 3: Implement**

Create `src/print/components/AcquisitionTimeline.tsx` — the placement wrapper (lives in `print/components/`, NON-exempt, so it carries NO appearance, only layout):
```tsx
import type { JSX } from "react";
import { AcquisitionChart } from "../plot/AcquisitionChart"; // render layer (appearance-exempt)

export interface AcqSession { label: string; loadFrameCounts: number[] }

export interface AcquisitionTimelineProps {
  sessions: AcqSession[];
  className?: string;
}

/** Per-session cluster of per-load exposure-count bars, drawn on the shared
 *  trace-plot family. Placement wrapper; appearance lives in the exempt render
 *  layer (print/plot/AcquisitionChart). */
export function AcquisitionTimeline({ sessions, className }: AcquisitionTimelineProps): JSX.Element {
  return (
    <div className={className} data-testid="acquisition-timeline">
      <AcquisitionChart sessions={sessions} />
    </div>
  );
}
```

Create `src/print/plot/AcquisitionChart.tsx` — the actual SVG bar chart (this file is design-guard-EXEMPT, so it authors `fill`/`stroke`/tokens directly). REAL code (not prose) — a linear y-scale over the max frame count, a categorical x-layout (session clusters of per-load bars), one `<rect data-testid="acq-bar">` per load, one `<text>` session label per session:
```tsx
import type { JSX } from "react";
import type { AcqSession } from "../components/AcquisitionTimeline";

const H = 132;          // chart height (px)
const PAD_TOP = 8;
const PAD_BOTTOM = 22;  // room for the session label row
const BAR_W = 14;
const BAR_GAP = 4;      // between bars in a cluster
const CLUSTER_GAP = 22; // between session clusters
const PLOT_H = H - PAD_TOP - PAD_BOTTOM;

export interface AcquisitionChartProps {
  sessions: AcqSession[];
  className?: string;
}

/** Categorical bar cluster (per session, per load) on the trace-plot token
 *  language. Appearance-exempt render layer — authors fill/stroke directly. */
export function AcquisitionChart({ sessions, className }: AcquisitionChartProps): JSX.Element {
  const maxCount = Math.max(1, ...sessions.flatMap((s) => s.loadFrameCounts));
  const yOf = (count: number): number => PLOT_H - (count / maxCount) * PLOT_H; // top edge of a bar

  // Walk left→right, accumulating x and recording per-cluster centers for labels.
  let x = 0;
  const clusters = sessions.map((s) => {
    const bars = s.loadFrameCounts.map((count) => {
      const bx = x;
      x += BAR_W + BAR_GAP;
      return { x: bx, count };
    });
    if (bars.length === 0) x += BAR_W; // empty session still occupies a slot
    const start = bars[0]?.x ?? x;
    const end = (bars[bars.length - 1]?.x ?? x) + BAR_W;
    const center = (start + end) / 2;
    x += CLUSTER_GAP;
    return { label: s.label, bars, center };
  });
  const width = Math.max(0, x - CLUSTER_GAP);

  return (
    <svg
      className={className}
      width="100%"
      viewBox={`0 0 ${Math.max(width, 1)} ${H}`}
      role="img"
      aria-label="Acquisition timeline — exposures per load, grouped by session"
    >
      {/* baseline hairline in the plot's frame-edge token */}
      <line x1={0} y1={PAD_TOP + PLOT_H} x2={width} y2={PAD_TOP + PLOT_H}
            stroke="var(--color-hair)" strokeWidth={1} />
      {clusters.map((c) => (
        <g key={c.label}>
          {c.bars.map((b, i) => (
            <rect
              key={i}
              data-testid="acq-bar"
              x={b.x}
              y={PAD_TOP + yOf(b.count)}
              width={BAR_W}
              height={PLOT_H - yOf(b.count)}
              rx={2}
              fill="var(--color-accent)"
            />
          ))}
          <text
            x={c.center}
            y={H - 6}
            textAnchor="middle"
            fontSize={10}
            fontWeight={700}
            fill="var(--color-ink-faint)"
            style={{ textTransform: "uppercase", letterSpacing: "0.04em" }}
          >
            {c.label}
          </text>
        </g>
      ))}
    </svg>
  );
}
```

> The test asserts (a) one `<rect data-testid="acq-bar">` per load across all sessions (5 for the fixture), (b) a session label per session (`getByText(/APR 12/i)` — the SVG `<text>` renders "Apr 12" upcased via CSS `text-transform`; **note `getByText` matches the DOM text content `"Apr 12"`, NOT the visual uppercasing**, so the test's `/APR 12/i` is case-insensitive and matches), (c) zero bars for an empty `sessions={[]}`. The `var(--color-*)` fills are permitted because this file is in the exempt `print/plot/` prefix. `frame_count` per load is the bar height input (the page maps `load.frame_count` → `loadFrameCounts`).

- [ ] **Step 4: Run → PASS + guard + tsc**
```bash
npm test -- test/AcquisitionTimeline.test.tsx
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS (the guard must stay green — `AcquisitionChart` lives in the exempt `print/plot/` prefix).

- [ ] **Step 5: Commit**
```bash
git add src/print/components/AcquisitionTimeline.tsx src/print/plot/AcquisitionChart.tsx test/AcquisitionTimeline.test.tsx
git commit -m "feat(config): AcquisitionTimeline on the trace-plot render layer"
```

---

## Task 19: `ConfigurationBody` — compose the Configuration tab internals

Fills the body of E1's `ExperimentConfigurationPage` shell. The Configuration tab shows, in order: **editable description** (E1 owns the `description TEXT` column, the `ExperimentShell` name/description edit-in-place, and the `ExperimentPatch.description` field — E2 RENDERS the description here via `useExperiment` and `updateMutate`, consuming E1's work; no longer "unowned/missing"), **Geometry** (Override/Revert ledger), **Acquisition** (read-only timeline), **Sources** (editable patterns + read-only dirs). Owns the geometry override/undo state (via `useUndoStack`), the geometry override → `PATCH /api/experiments/{id}` mutation (`updateExperiment`), and the field-revert (revert to the derived value). Derives the geometry rows from the experiment detail (`*_source` fields) and the acquisition sessions from `useLoads`.

**Files:**
- Create: `src/print/components/ConfigurationBody.tsx`
- Test: `test/ConfigurationBody.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/ConfigurationBody.test.tsx`:
```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { ConfigurationBody } from "../src/print/components/ConfigurationBody";

const updateMutate = vi.fn();

vi.mock("../src/queries", async (orig) => {
  const actual = await orig<typeof import("../src/queries")>();
  return {
    ...actual,
    // E1's Experiment carries PER-FIELD *_source columns (NOT a combined
    // beam_center_source): energy_kev_source, flight_path_m_source,
    // beam_center_x_source, beam_center_y_source, pixel_size_um_source,
    // q_units_source — plus ingest_status/scan_signature/last_scanned_at.
    useExperiment: () => ({
      data: {
        id: 7, name: "SSRL · 1p7m",
        description: "April 2026 beamtime at SSRL 1p7m",
        path: "/d", data_dir: "/d", analysis_dir: "/d/analysis",
        image_pattern: "{name}.tiff", metadata_pattern: "{name}.prp", integration_pattern: "{name}.dat",
        manifest_path: null, created_at: "2026-04-12", q_units: "A^-1",
        energy_kev: 9.0, energy_kev_source: "prp",
        flight_path_m: 1.81, flight_path_m_source: "setup",
        beam_center_x: 421.4, beam_center_x_source: "setup",
        beam_center_y: 836.9, beam_center_y_source: "setup",
        pixel_size_um: 172, pixel_size_um_source: "prp", q_units_source: "prp",
        last_scanned_at: "2026-04-12T10:00:00", scan_signature: "sig", ingest_status: "idle",
      },
      isLoading: false,
    }),
    useLoads: () => ({ data: [], isLoading: false }),
    useUpdateExperiment: () => ({ mutate: updateMutate, isPending: false }),
  };
});

const wrap = (n: ReactNode) => render(<QueryClientProvider client={new QueryClient()}>{n}</QueryClientProvider>);

describe("ConfigurationBody", () => {
  it("renders the description, Geometry, Acquisition, and Sources cards", () => {
    wrap(<ConfigurationBody experimentId={7} />);
    // description (DECISION: E1-owned column; E2 renders it here)
    expect(screen.getByText(/April 2026 beamtime/)).toBeInTheDocument();
    expect(screen.getByText("Geometry")).toBeInTheDocument();
    expect(screen.getByText("Acquisition")).toBeInTheDocument();
    expect(screen.getByText("Sources")).toBeInTheDocument();
  });
  it("derives a geometry row per typed field with its *_source provenance", () => {
    wrap(<ConfigurationBody experimentId={7} />);
    expect(screen.getByText("Beam energy")).toBeInTheDocument();
    expect(screen.getByText("setup files")).toBeInTheDocument(); // flight path + beam center source
  });
  it("renders the 3 pattern rows as editable (E1 columns: image/metadata/integration_pattern)", () => {
    wrap(<ConfigurationBody experimentId={7} />);
    expect(screen.getByText("{name}.tiff")).toBeInTheDocument();  // image_pattern
  });
});
```

> **Adapt to E1's actual `Experiment` shape (verify before writing).** The mock mirrors E1's `Experiment` extension: **per-field `*_source` columns** (`energy_kev_source`, `flight_path_m_source`, `beam_center_x_source`, `beam_center_y_source`, `pixel_size_um_source`, `q_units_source`) — there is **NO combined `beam_center_source`** and **NO `geometry_discrepancies` field**. The detail hook is `useExperiment` (E1). **Per the DECISION:** E1 adds `description TEXT` (nullable) and the 3 typed pattern columns (`image_pattern`, `metadata_pattern`, `integration_pattern`) via additive migration — E2 COUNTS on them existing. If E1 has not landed these columns yet, STOP and coordinate. `discrepancyCount` is `0` for v1 (DECISION: deferred, no backend field exists in Phase A–C — see Step 3).

- [ ] **Step 2: Run → FAIL.** `npm test -- test/ConfigurationBody.test.tsx`

- [ ] **Step 3: Implement**

Create `src/print/components/ConfigurationBody.tsx`:
- Read `useExperiment(experimentId)` + `useLoads(experimentId)`; instantiate `useUpdateExperiment(experimentId)`.
- Build `GeometryRow[]` from the typed fields, **each with its own per-field `*_source`** (spec §8.1/§8.7 mandate per-field provenance):
  - `energy_kev` → "Beam energy" (`energy_kev_source`), value `\`${e.energy_kev?.toFixed(2)} keV\``
  - `flight_path_m` → "Flight path" (`flight_path_m_source`)
  - `beam_center_x` → "Beam center X" (`beam_center_x_source`), value `\`${e.beam_center_x?.toFixed(1)} px\``
  - `beam_center_y` → "Beam center Y" (`beam_center_y_source`), value `\`${e.beam_center_y?.toFixed(1)} px\``
  - `pixel_size_um` → "Pixel pitch" (`pixel_size_um_source`)
  - `q_units` → "q units" (`q_units_source`)
- Use `useUndoStack<{ key; prevValue: number | string; prevSource }>` to back both per-field Revert AND the header "Undo last change". On Override commit: push `{ key, prevValue, prevSource }`, then `updateMutate({ <patchKey>: parsed })`. On Revert/Undo: `pop()` → re-`updateMutate` to the prev value (or, if Phase D supports it, clear the user source so the next scan re-derives — confirm).
- **`discrepancyCount`**: no discrepancy field exists in Phase A–C; pass `discrepancyCount={0}` for v1 and add a `// TODO(Phase-D/E1): source from geometry_discrepancy when the field lands` comment. The banner is deferred. Do NOT invent a `geometry_discrepancies` field on `Experiment` or compute it client-side — the spec §9.6 discrepancy detection lives in the backend scan path (not yet built in A–C) and E1's `Experiment` type carries no such column.
- Build `AcqSession[]` from `loads`: group by `start_time`'s date → session label; `loadFrameCounts` = each load's `frame_count`.
- Build `SourceRow[]` with `data_dir`/`analysis_dir` as read-only rows (no edit affordance — DECISION: directory is fixed at create) and `image_pattern`/`metadata_pattern`/`integration_pattern` as editable rows (DECISION: E1 adds these as typed columns; editing a pattern calls `updateMutate({ image_pattern: … })` via `ExperimentPatch` and the backend handles the rescan). The `ExperimentPatch` type (DECISION rename from `ExperimentGeometryPatch`) covers all three pattern fields; no separate PATCH route needed.
- **Render the editable description** (DECISION: E1 owns the `description TEXT` column and the `ExperimentShell` name/description edit-in-place; E2 CONSUMES it here): read `experiment.description` from `useExperiment` and render it as a plain editable text area (or reuse E1's pattern if it exposed a standalone description editor) at the TOP of the Configuration body — before Geometry. On commit, call `updateMutate({ description })` (no rescan). If `experiment.description` is null/empty, show a placeholder ("Add a description…"). This is a plain write with no rescan (DECISION).
- Compose, in order: (1) editable description at the top, (2) a two-column grid (`grid grid-cols-2 gap-4` — layout utilities) with `GeometryLedger` + a `Card` wrapping `<h3>Acquisition</h3>` + `AcquisitionTimeline`, (3) `SourcesCard` full-width below.

> **Override mutation:** `updateExperiment`'s patch type is E1's `ExperimentPatch` (DECISION rename from `ExperimentGeometryPatch`; covers name/description/geometry ×6/patterns ×3). Map each geometry row key → its patch key + parse the display string back to a number ("9.00 keV" → `9.0`). Tolerant parsing; if unparseable, don't mutate (keep the inline edit open / mark invalid).

- [ ] **Step 4: Run → PASS + guard + tsc**
```bash
npm test -- test/ConfigurationBody.test.tsx
node scripts/check-design.mjs
npx tsc --noEmit -p tsconfig.json
```
→ PASS.

- [ ] **Step 5: Wire into E1's `ExperimentConfigurationPage`**

In E1's `ExperimentConfigurationPage`, render `<ConfigurationBody experimentId={id} />` below the shared header. (If E1's page already renders a placeholder body, replace it.) Re-run E1's page test if one exists.

- [ ] **Step 6: Commit**
```bash
git add src/print/components/ConfigurationBody.tsx test/ConfigurationBody.test.tsx
git commit -m "feat(config): ConfigurationBody composes geometry/acquisition/sources with override+undo"
```

---

## Task 20: Mount `GroupingReviewPage` into the route + Playwright end-to-end (conditional on E1)

Wire `GroupingReviewPage` into E1's `/experiments/:id/grouping` route (the route element comes from E1; this task supplies the body element if E1 left it as a placeholder), and add a mocked Playwright spec for the merge / split / move / undo + filter-then-select flow.

**Files:**
- Modify: E1's grouping route registration (or `AppRoutes.tsx`) to render `GroupingReviewPage`
- Create: `e2e/grouping-review.spec.ts`
- Test: the Playwright spec

- [ ] **Step 1: Mount the page**

Read E1's `AppRoutes.tsx` `/experiments/:id/grouping` route. If it renders a placeholder, replace with a small wrapper that reads `:id` (via `useParams`) and renders `<GroupingReviewPage experimentId={Number(id)} onBack={() => navigate('../corpus')} />`. If E1 already mounts it, skip.

- [ ] **Step 2: Write the mocked Playwright spec**

Create `e2e/grouping-review.spec.ts` (follow the repo's mocked-E2E pattern in `e2e/` — route-mock `/api/experiments/:id/loads`, the structural-edit endpoints, and the SSE stream). Assert:
1. Default `/experiments/7/grouping` shows only flagged loads; "All loads" reveals all.
2. A glob filter (`HA8*`) narrows; selecting a sample then changing the filter keeps the selection (bulk bar count stays).
3. Clicking a merge-prompt "Merge into that sample" → confirm → the loser sample disappears and an undo toast shows.
4. A split divider's "Split here" → confirm → the sample splits into two.
5. A per-exposure Move opens the picker and moves the exposure.

Mock the structural endpoints to 200/201. The grouping mutators invalidate `loads` in `onSuccess` (own-op), so the optimistic edit + a second `/loads` GET (returning the post-edit roll-up) is sufficient for the mocked flow — you do NOT need to emit SSE for own-op reconcile. (Emit a structural SSE frame only if you also test the FOREIGN-tab path.) **Use the §9.3 event kinds if you do emit SSE** — `exposure_moved`/`sample_renamed`/`sample_created`/`sample_split`/`grouping_flag_dismissed` on the `curation` channel; **NOT `sample_merged`** (which does not exist — a merge is an `exposure_moved` burst).

> **Follow `e2e/AGENTS.md`** for selectors (use `getByRole`/`data-testid`, not class strings), port binding, and SSE-mock timing. The grouping flow's optimistic edits land before any refetch; assert on the optimistic DOM, then on the post-invalidate `/loads` refetch.

- [ ] **Step 3: Run the spec**

`npm run e2e -- grouping-review` (mocked; auto-starts Vite). Expected: PASS.

- [ ] **Step 4: Commit**
```bash
git add e2e/grouping-review.spec.ts src/print/shell/AppRoutes.tsx
git commit -m "test(e2e): grouping-review merge/split/move/undo + filter-then-select"
```

---

## Self-Review

**Verified against LIVE source (2026-06-18) — every API below was read, not assumed:**
- **`api.ts` HTTP layer:** `request<T>(method, path, body?, opts?)` is the SOLE fetch abstraction (`api.ts:77-104`) — there is **NO `getJSON`/`postJSON`/`patchJSON`**. `AuthOpts = { username?; clientId?; clientOpId? }` (`:71-75`); headers `X-Username`/`X-Client-Id`/`X-Client-Op-Id` are built at `:83-87`. All E2 fetchers use `request<T>` so the op-id threads. `Exposure` is exported (`:151-170`).
- **client id + hooks:** `getClientId()` is a FUNCTION in `src/lib/clientId.ts` (sessionStorage-backed) — **NO `useClientId()` hook**. `queries.ts:42` holds the module-level `const CLIENT_ID = getClientId()`; hooks pass `clientId: CLIENT_ID` inline (`useUpdateSample:520-526`, `useAddPeak:372-382`). No `authOpts(...)` at the hook layer; `authOpts(username, clientId, clientOpId)` is the MUTATOR-layer helper at `src/lib/authOpts.ts` (imported as `../../authOpts`).
- **Queue contract:** `Mutator<TInput,TScope,TResponse>` (`types.ts:193-223`) — `onMutate(payload, qc) → RollbackContext`, `request(payload, signal)`, `onSuccess(payload, response, qc)`, optional `synthesizeFromSse`. `RollbackContext = { restore }` (`:77-79`). `SseEvent` has required `entity_id`, `payload?: unknown` (`:146-157`). `optimisticId.ts:50` exports `nextOptimisticId()` (negative). `persistence.ts:14` `SCHEMA_VERSION = 4` (E1 bumps to 5).
- **`replayCoordinator.ts` (FIX 3 anchor):** the own-op SSE-wins path resolves the deferred + aborts the HTTP and `return`s — it **never calls `applyRemoteToCache`** (only `applyPostStateOnly`); the self-echo path also returns early. Only FOREIGN events call `applyRemoteToCache`. → every grouping mutator's `onSuccess` MUST invalidate `loads(id)` itself (the `saveSeries.ts` `series_save` precedent, which invalidates in `onSuccess`). DONE for all five mutators (Tasks 4/5).
- **`applyRemoteToCache.ts`:** `default:` arm (`:334-337`) fires `peaks(id)`/`indices(id)` — the poisoning E1's structural arms guard against by landing before it; **Task 7 only VERIFIES E1's arms** (the structural RECEIVE arms are E1-owned, per the cross-plan ownership reconciliation — E2 does not edit this file for structural events). The `update_sample` arm (`:314-318`) untyped-spreads `payload` — the `display_name` cross-version hazard, but its fix is **E1's** (Task 0 verification-only).
- **`mutatorRegistry.ts`:** `resolveMutator` (`:58-134`) + `resolveMutatorForEvent` (`:147-207`) dual maps; `saveSeries.ts` maps `series_save` ↔ `series_created`/`series_recipe_updated` and invalidates in `onSuccess` (the precedent).
- **`queryKeys`:** live `samples(id) = ["experiment", id ?? "none", "samples"]` — `loads(id)` mirrors it as `["experiment", id ?? "none", "loads"]` (E1-owned). The earlier E2 `["loads", id]` is dropped.
- **UI primitives (every prop verified):** `Toast`/`lib/toast.ts` is a NO-PROP imperative singleton (`showToast(msg, kind)`, item `{id, msg, kind}`) — Task 8 extends the imperative signature, NOT a `toasts` prop. `SearchInput` uses `onChange`+`ariaLabel` (NOT `onValueChange`/`label`). `IconButton` requires `label` (omits `aria-label`). `Checkbox` is `checked`/`onChange(checked)` + `aria-label`, renders `role="checkbox"`. `Button` variants `solid|accent|success|ghost|danger|outline|ghostInverse` — **no `size` prop, no `out-accent`**; `ghostInverse` is the dark-surface variant. `Input` has `variant="title"`/`mono`/`testId`/`onValueChange` + spreads `InputHTMLAttributes`. `Card` takes `className`/`children` (+ optionals). `Menu`/`SegmentedControl`/`Kicker`/`EmptyState`/`ModalShell` exported. `Thumbnail` takes `src: string|null` + `size: "xs"|"sm"|"lg"` — **NOT `exposureId`** (→ `ExposureLeaf` takes `thumbSrc`).
- **`CullBar` (FIX 7):** hardcoded corpus-culling verbs (`onReject/onKeep/onRestore/onClear`), `data-testid="cull-bar"` baked in, no parameterized primary — so E2 builds a small local `GroupingBulkBar` instead of bending CullBar.
- **plot engine (FIX 12):** `print/plot/` exports `PlotFrame`/`Axis`/`makeProjection`/`makeAxis`/`positiveExtent` (reusable scale/frame) but is **trace-only — NO bar renderer, NO `cleanFigureSvg`**. → `AcquisitionChart` is a thin hand-authored SVG in the exempt `print/plot/` prefix (real code written in Task 18).
- **design guard (`scripts/check-design.mjs`):** exempt prefixes are `print/ui/`, `print/plot/`, `print/detector/`, `print/comb/`, `print/export/` — **`print/components/` is NOT exempt** (every leaf in this plan lives there). Bans `text-[…]`/`rounded-[…]`/raw-colour-utility/raw-colour-literal/side-stripe(`border-l-*`). `--color-accent-wash` does NOT exist (FIX 8 → add the token or use a primitive).
- **`SeriesScopingPage.tsx`:** `HistoryEntry` 4-variant union at `:77-87`; the live undo (`:492-514`) nests restore `setState`s inside the `setHistory` updater (the workaround); the StrictMode-double-push lesson is commented at `:349-352`. `useUndoStack` adopts the `pop()` (effect-outside-updater) shape — the StrictMode-safe variant by construction.
- **spec §8.8 shapes (authoritative, replacing the earlier E2 draft):** exposure key `.id` (not `exposure_id`); no `frame_no`/`status`; sample has `merged_into_id`, `grouping_source`/`name_source: string`; load keyed `load_id` (no top-level `id`/`experiment_id`). ALL fixtures rewritten to this.
- **spec §9.3 event kinds (FIX 2):** **no `sample_merged`** (merge fans out as `exposure_moved`); split emits `sample_created` + `sample_split`; `grouping_flag_dismissed` is durable. The cache arms, the registry event map, and the dismiss-as-mutation all reflect this.

**Fixes applied (mapped to the team-lead list):** 1 (E2 imports E1's `Load`/`queryKeys.loads`, all fixtures on `["experiment",id,"loads"]`, nested §8.8 shape) · 2 (no `deriveFlag`; reads `sample.flag`; dismiss is a durable `useDismissGroupingFlag` mutator + arm) · 3 (all five mutators invalidate `loads` in `onSuccess` for own-op reconcile; the foreign-tab `applyRemoteToCache` arms are E1-owned, Task 7 verifies them) · 4 (`request<T>` + `CLIENT_ID`) · 5 (imperative `showToast` action slot) · 6 (SearchInput `onChange`/`ariaLabel`, IconButton `label`, Checkbox `role=checkbox`) · 7 (local `GroupingBulkBar`; coverage via page-level filter test `getByTestId('bulk-bar')`) · 8 (`--color-accent-wash` token or primitive; no `bg-[color-mix]`) · 9 (`useUpdateExperiment` uses `ExperimentPatch` (DECISION rename); per-field `*_source` rows; `discrepancyCount={0}` v1 deferred) · 10 (`MergeSamplesInput = {loserId,survivorId}`, no `as never`) · 11 (Task 0 verification-only; `update_sample` arm + SCHEMA bump are E1's) · 12 (real `AcquisitionChart.tsx` code, exempt prefix) · 13 (StrictMode-safe `pop()` hook; `px-4.5`→`px-4`; renamed `response.name` guard documented as own-op partial) · 14 (spec §9.6 hook names throughout). **DECISION revisions (2026-06-19):** description is E1-owned/E2-renders (Task 19); SourcesCard pattern rows DEFINITE editable (E1 adds columns); `data_dir`/`analysis_dir` READ-ONLY; `discrepancyCount={0}` v1 (deferred); `ExperimentGeometryPatch` renamed `ExperimentPatch` throughout. **Hook-arity fix (2026-06-19, MUST-FIX 1):** all row-scoped hooks (`useRenameSample`, `useMoveExposure`, `useSplitSample`, `useDismissGroupingFlag`) now take ONLY `experimentId`; entity ids (`sampleId`/`exposureId`) moved into the corresponding mutator Input types and the `mutate()` call — one hook instance per experiment, not per row. Task 14's re-decision paragraph removed; it now references the settled fact. **Optional guard tests removed (CUT 1):** Tasks 0/2/6/7 Step-2 guard test files (`labelCollapse.test.ts`, `loadsTypes.test.ts`, `loadsKey.test.ts`, `applyRemoteToCache.structural.e2guard.test.ts`) removed — the Step-1 grep/Read verification is the actual gate; these files tested E1's contracts, not E2 behavior.

**Live APIs that DIFFERED from the original E2 draft's assumptions (the substantive corrections):**
- `request<T>` is the only fetcher helper (draft assumed `getJSON`/`postJSON`/`patchJSON`).
- `getClientId()` function + module `CLIENT_ID` (draft assumed a `useClientId()` hook).
- Toast is a NO-PROP imperative singleton with `{id,msg,kind}` (draft assumed a `toasts`-prop container with `{message,tone}`).
- `Button` has no `size` prop; `SearchInput` uses `onChange`/`ariaLabel`; `IconButton` uses `label`; `Thumbnail` takes `src` not `exposureId` — all draft JSX corrected.
- `CullBar` is not parameterizable without contract bloat → local bar.
- `print/plot/` is trace-only (no bar frame) → hand-authored SVG.
- `--color-accent-wash` does not exist → token-or-primitive.
- The §8.8 Load shape (`.id`/`merged_into_id`/`load_id`) and the §9.3 kind set (no `sample_merged`; `sample_created`; durable dismiss) differ from the draft's `exposure_id`/`frame_no`/`status`/`["loads",id]`/`sample_merged` — all rewritten.

**Open questions the implementer must still resolve against E1/Phase-D:**
1. **E1 must ship the NESTED §8.8 `Load`** (its Task 1 currently shows a flat one) — flagged as an E1 bug to fix in E1, not to route around. Tasks 0/2/6 are verification gates that STOP if E1's contract is wrong.
2. **Phase D route paths + bodies** for merge/split/dismiss + the split event burst ordering (which frame of `sample_created`/`sample_split` carries the op's `client_op_id`) — trace before finalizing Tasks 3/5.
3. **Discrepancy banner** (Task 19) — **RESOLVED: pass `discrepancyCount={0}` for v1, deferred** (no backend discrepancy field exists in Phase A–C; the banner does not render). No open question.
4. **`thumbSrcFor`** — reuse the existing detector-image-URL helper if one exists; `() => null` (placeholder) acceptable for v1.
5. **Merge/split undo** — v1 shows their toasts WITHOUT an undo action (one stamp can't reverse a multi-row reassignment; spec §9.3); only rename/move get the undo action. Flag to Jonathan if product wants more.
6. **`structuredClone`** in the test env — assumed present (Node 18+/JSDOM); verify.

All referenced types/functions/props are defined: the five `OpKind`s (Task 3); `Load`/`LoadSample`/`LoadExposure`/`GroupingFlag`/`listLoads`/`queryKeys.loads` IMPORTED from E1 (Tasks 2/6, verification gates); the five mutators (Tasks 4–5); the six hooks (Task 9); the components (Tasks 10–19). The cross-plan dependency surface (E1) is gated by explicit verification tasks that STOP rather than redefine.
