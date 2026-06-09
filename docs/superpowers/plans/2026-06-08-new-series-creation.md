Now I have a complete, grounded picture. I'll produce the full implementation plan.

---

# New-Series Creation Flow Implementation Plan

> For agentic workers: read this file top-to-bottom before touching any file. Each task is self-contained: write the failing test first, run it, implement the minimum to pass, run again, then commit. Never skip the gate commands at the end of each task.

## Goal

Add an end-to-end "new series creation" flow (milestone M-A) to the HimalayaUI greenfield frontend. The flow is B-shaped:

1. A sample-grain checkbox column on the contact sheet lets the user check whole samples (distinct from the frame-grain CullBar).
2. A `ComposeBar` floating bar appears when one or more samples are checked, offering `+ New series` and `Clear`.
3. `+ New series` navigates to `/series/new` with the selected sample IDs carried as route location state (not Zustand — see Architecture).
4. `SeriesScopingPage` reads the seed and narrows its picker query to those samples; when no ordering variable can be proposed, it enters a cold-corpus path where the user names the variable and assigns per-sample values inline.
5. `Confirm & build` writes the scoping tags AND creates the series via `POST /api/series`, then navigates to `/series/:id`.

## Architecture Decision — Selection-Carry Mechanism

**Recommendation: React Router `location.state`**, not a new Zustand slice.

Rationale: the selected sample IDs are one-shot navigational intent — they are consumed once on arrival at `/series/new` and then discarded. They have no meaning to any other surface and must not survive a tab refresh or a bookmark. Zustand's `persist` middleware partialises to localStorage; even without partialising, Zustand state would leak across the back-button or a stale tab. `location.state` is exactly a React Router ephemeral navigation payload: it is scoped to the history entry, not visible until you arrive, and silently absent on a direct `/series/new` visit (which then falls back to the full corpus query — the existing code path).

Trade-off acknowledged: `location.state` is untyped at the call site without a helper. Mitigated in Task 3 by a typed `navigateToNewSeries(ids, navigate)` helper that encapsulates the cast and keeps the type contract in one place.

Comparison with Zustand slice: a new `newSeriesSeed: number[] | null` field in `state.ts` would require: a new action, a clear-on-arrive action, wiring into `partialize`, and a version bump. That is infrastructure work for a one-way channel — overkill. Location state is the conventional React Router answer for exactly this pattern.

## Tech Stack

- React 18 + React Router v6 (location.state), TypeScript strict + `exactOptionalPropertyTypes`
- `src/print/ui/` closed-look primitives (new `Checkbox` + new `ComposeBar`)
- `src/print/pages/SamplesPage.tsx` and `SeriesScopingPage.tsx` (page edits)
- `src/print/pages/scopingDerive.ts` (logic extension)
- `src/lib/queue/mutators/scopeSeries.ts` (extended to also create the series)
- `src/queries.ts` (`useCreateAndScopeSeries` new composite hook, or inline in page)
- Vitest unit tests + one Playwright E2E spec for the full create path
- Gate on every task: `npx tsc -p tsconfig.build.json --noEmit` + `npm run lint:design` + targeted `npm test` / `npm run e2e`

## File Structure

| File | Status | Responsibility |
|---|---|---|
| `src/print/ui/Checkbox.tsx` | CREATE | Closed-look checkbox primitive — indeterminate + checked states, `data-checked`, `data-indeterminate`, `aria-checked`. No inline appearance literals. |
| `src/print/ui/Checkbox.stories.tsx` | CREATE | Storybook coverage for the new primitive. |
| `src/print/ui/index.ts` | MODIFY | Add `Checkbox` export. |
| `src/print/ui/ComposeBar.tsx` | CREATE | Floating compose bar — `count`, `show`, `onNewSeries`, `onClear`. Dark bg-ink bar (parallel to `CullBar` pattern). Honest: `onNewSeries` button only renders / is enabled when `count > 0`. |
| `src/print/components/SheetTable.tsx` | MODIFY | Accept `checkboxColumn?: boolean` prop; when true, prepend a `36px` fixed track to `SAMPLE_TABLE_COLS` and render a column-header checkbox cell. The `SAMPLE_TABLE_COLS` constant becomes a function `sampleTableCols(withCheckbox)` so header and rows stay aligned. |
| `src/print/components/SampleTableRow.tsx` | MODIFY | Accept `checked?: boolean`, `indeterminate?: boolean`, `onCheck?: () => void`; when provided, render a leading checkbox cell aligned to the new grid track. `SAMPLE_TABLE_COLS` export updated to match. |
| `src/print/pages/SamplesPage.tsx` | MODIFY | Add `checkedSamples: Set<number>` state (sample-grain, distinct from exposure-grain `selected`). Wire `onCheck` on each row. Mount `ComposeBar` alongside `CullBar`. `+ New series` calls `navigateToNewSeries(checkedSamples, navigate)`. |
| `src/print/pages/samplesAdapters.ts` | NO CHANGE | No adapter change needed — `checked` is wired directly at the page. |
| `src/lib/series/newSeriesNav.ts` | CREATE | `navigateToNewSeries(sampleIds: Set<number>, navigate: NavigateFunction): void` — typed helper that encapsulates the `navigate("/series/new", { state: { seedSampleIds: [...sampleIds] } })` call. Also exports the `NewSeriesLocationState` type + `readNewSeriesSeed(location: Location): number[] | null` reader. |
| `src/print/pages/SeriesScopingPage.tsx` | MODIFY | Read `readNewSeriesSeed(location)` on mount; when a seed is present, pass a `seedSampleIds` prop to `proposeOrdering` (after filtering `pickerQ.data`). Handle the cold-corpus path (no `orderingKey`) with a new `AssignPanel` inline rather than the dead-end `EmptyState`. After `handleBuild` success, create the series and navigate to `/series/:id` (not `/series`). |
| `src/print/pages/scopingDerive.ts` | MODIFY | Add `buildColdAssignRows` (pure function: seed ids → `ColdAssignRow[]` with name + editable value). Add `canColdBuild(rows)` gate. |
| `src/print/components/ColdAssignPanel.tsx` | CREATE | L3 composite: ordered list of sample rows each with a name label and an `Input` field for the ordering-variable value. Props: `rows`, `variableKey`, `onKeyChange`, `onValueChange`. No appearance literals outside `src/print/ui/**`. |
| `src/lib/queue/mutators/scopeAndCreateSeries.ts` | CREATE | A new queue-routed mutator that (1) calls `batchSampleTags` (same as `scopeSeriesMutator`) and then (2) calls `api.saveSeries` (create, no id). Returns the new `Series`. On success, splices the series into `queryKeys.seriesList` cache and invalidates `corpusSampleTags` + `corpusPickerSamples`. |
| `src/queries.ts` | MODIFY | Export `useScopeAndCreateSeries()` hook backed by the new mutator. |
| `e2e/new-series-creation.spec.ts` | CREATE | Full Playwright E2E: check samples on `/samples` → compose bar appears → `+ New series` → `/series/new` loads pre-seeded → confirm-and-build POSTs tags + POSTs series → navigates to `/series/:id`. |
| `test/Checkbox.test.tsx` | CREATE | Vitest unit tests for the `Checkbox` primitive. |
| `test/ComposeBar.test.tsx` | CREATE | Vitest unit tests for the `ComposeBar` primitive. |
| `test/SamplesPage.samplePicker.test.tsx` | CREATE | Unit test: checking samples raises compose bar; compose bar hidden at zero count. |
| `test/scopingDerive.coldPath.test.ts` | CREATE | Unit tests for `buildColdAssignRows` + `canColdBuild`. |
| `test/newSeriesNav.test.ts` | CREATE | Unit tests for `navigateToNewSeries` / `readNewSeriesSeed`. |
| `test/queue/scopeAndCreateSeries.contract.test.ts` | CREATE | Contract test for the new mutator (mirrors `scopeSeries.contract.test.ts` pattern). |

---

## Build Sequence

### Task 1 — `Checkbox` primitive in `src/print/ui/`

**Why first:** `SheetTable` and `SampleTableRow` both depend on it. Nothing else can be checked before this closed-look primitive exists. Also satisfies the design-guard rule: appearance must live in `src/print/ui/`, not inline in a consumer.

**Files:** create `src/print/ui/Checkbox.tsx`, create `test/Checkbox.test.tsx`, modify `src/print/ui/index.ts`.

- [ ] Write `test/Checkbox.test.tsx`:

```tsx
// test/Checkbox.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { Checkbox } from "../src/print/ui/Checkbox";

describe("Checkbox", () => {
  it("renders unchecked by default with correct aria-checked", () => {
    render(<Checkbox aria-label="select sample" />);
    const cb = screen.getByRole("checkbox", { name: "select sample" });
    expect(cb).toHaveAttribute("aria-checked", "false");
    expect(cb).toHaveAttribute("data-checked", "false");
  });

  it("renders checked state", () => {
    render(<Checkbox checked aria-label="select sample" />);
    const cb = screen.getByRole("checkbox", { name: "select sample" });
    expect(cb).toHaveAttribute("aria-checked", "true");
    expect(cb).toHaveAttribute("data-checked", "true");
  });

  it("renders indeterminate state", () => {
    render(<Checkbox indeterminate aria-label="select sample" />);
    const cb = screen.getByRole("checkbox", { name: "select sample" });
    expect(cb).toHaveAttribute("aria-checked", "mixed");
    expect(cb).toHaveAttribute("data-indeterminate", "true");
  });

  it("calls onChange when clicked", async () => {
    const user = userEvent.setup();
    const onChange = vi.fn();
    render(<Checkbox onChange={onChange} aria-label="select sample" />);
    await user.click(screen.getByRole("checkbox", { name: "select sample" }));
    expect(onChange).toHaveBeenCalledTimes(1);
  });

  it("does not call onChange when disabled", async () => {
    const user = userEvent.setup();
    const onChange = vi.fn();
    render(<Checkbox disabled onChange={onChange} aria-label="select sample" />);
    await user.click(screen.getByRole("checkbox", { name: "select sample" }));
    expect(onChange).not.toHaveBeenCalled();
  });
});
```

- [ ] Run `node_modules/.bin/vitest run test/Checkbox.test.tsx` — expect "Cannot find module" / component-not-found failure.

- [ ] Implement `src/print/ui/Checkbox.tsx`:

```tsx
// src/print/ui/Checkbox.tsx
import type { HTMLAttributes } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface CheckboxProps extends Omit<HTMLAttributes<HTMLSpanElement>, "onChange"> {
  checked?: boolean;
  indeterminate?: boolean;
  disabled?: boolean;
  onChange?: (checked: boolean) => void;
  /** Accessible name — required when no visible label wraps this. */
  "aria-label"?: string;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * Checkbox — a closed-look tri-state checkbox primitive (src/print/ui/).
 *
 * Implemented as a `<span role="checkbox">` rather than a native `<input
 * type="checkbox">` so the visual ring is fully token-driven (bg-ink,
 * border-hair-strong, accent) and no browser reset is needed. The span
 * receives `tabIndex={0}`, responds to Space/Enter, and announces via
 * `aria-checked` (false / true / mixed for indeterminate). `data-checked`
 * and `data-indeterminate` are stable E2E selectors.
 *
 * Appearance is intentionally minimal: a 16px box, border-hair-strong resting
 * ring, bg-ink fill + paper SVG check when checked, accent fill + dash when
 * indeterminate. No inline colour literals or arbitrary sizes — all from @theme.
 */
export function Checkbox({
  checked = false,
  indeterminate = false,
  disabled = false,
  onChange,
  className,
  ...rest
}: CheckboxProps): JSX.Element {
  const isChecked = !indeterminate && checked;
  const ariaChecked: "true" | "false" | "mixed" = indeterminate ? "mixed"
    : checked ? "true" : "false";

  function activate(): void {
    if (!disabled) onChange?.(!checked);
  }

  return (
    <span
      role="checkbox"
      aria-checked={ariaChecked}
      aria-disabled={disabled || undefined}
      tabIndex={disabled ? -1 : 0}
      data-checked={isChecked ? "true" : "false"}
      data-indeterminate={indeterminate ? "true" : undefined}
      onClick={activate}
      onKeyDown={(e) => {
        if (e.key === " " || e.key === "Enter") { e.preventDefault(); activate(); }
      }}
      className={cx(
        "inline-flex items-center justify-center flex-shrink-0 w-4 h-4 rounded border cursor-pointer select-none transition-colors",
        disabled && "opacity-40 cursor-not-allowed",
        isChecked ? "bg-ink border-ink"
          : indeterminate ? "bg-accent border-accent"
          : "bg-plate border-hair-strong hover:border-ink",
        className,
      )}
      {...rest}
    >
      {isChecked && (
        <svg viewBox="0 0 16 16" className="w-full h-full" aria-hidden="true">
          <path
            d="M4 8.2l2.6 2.6 5.4-5.6"
            fill="none"
            stroke="var(--color-paper)"
            strokeWidth={1.8}
            strokeLinecap="round"
            strokeLinejoin="round"
          />
        </svg>
      )}
      {indeterminate && (
        <svg viewBox="0 0 16 16" className="w-full h-full" aria-hidden="true">
          <line x1="4" y1="8" x2="12" y2="8"
            stroke="var(--color-paper)" strokeWidth={2} strokeLinecap="round" />
        </svg>
      )}
    </span>
  );
}
```

- [ ] Add to `src/print/ui/index.ts`:
  ```ts
  export { Checkbox } from "./Checkbox";
  export type { CheckboxProps } from "./Checkbox";
  ```

- [ ] Run `node_modules/.bin/vitest run test/Checkbox.test.tsx` — all pass.
- [ ] Run `npx tsc -p tsconfig.build.json --noEmit` — clean.
- [ ] Run `npm run lint:design` — clean.
- [ ] Commit: `feat(ui): Checkbox primitive — tri-state, token-only appearance`

---

### Task 2 — `ComposeBar` primitive + unit tests

**Why here:** it is a closed-look primitive that `SamplesPage` will mount. Must be built before the page is modified. Mirrors the `CullBar` pattern: presentational, no local state, dark bar.

**Files:** create `src/print/ui/ComposeBar.tsx`, create `test/ComposeBar.test.tsx`, modify `src/print/ui/index.ts`.

- [ ] Write `test/ComposeBar.test.tsx`:

```tsx
// test/ComposeBar.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { ComposeBar } from "../src/print/ui/ComposeBar";

describe("ComposeBar", () => {
  it("is hidden (aria-hidden) when show=false", () => {
    render(<ComposeBar count={0} show={false} />);
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("aria-hidden", "true");
  });

  it("is visible when show=true", () => {
    render(<ComposeBar count={2} show />);
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");
    expect(screen.getByTestId("compose-bar")).not.toHaveAttribute("aria-hidden", "true");
  });

  it("displays the sample count", () => {
    render(<ComposeBar count={3} show />);
    expect(screen.getByTestId("compose-bar")).toHaveTextContent("3");
  });

  it("calls onNewSeries when + New series is clicked", async () => {
    const user = userEvent.setup();
    const onNewSeries = vi.fn();
    render(<ComposeBar count={2} show onNewSeries={onNewSeries} />);
    await user.click(screen.getByRole("button", { name: /new series/i }));
    expect(onNewSeries).toHaveBeenCalledTimes(1);
  });

  it("calls onClear when Clear is clicked", async () => {
    const user = userEvent.setup();
    const onClear = vi.fn();
    render(<ComposeBar count={1} show onClear={onClear} />);
    await user.click(screen.getByRole("button", { name: /clear/i }));
    expect(onClear).toHaveBeenCalledTimes(1);
  });

  it("+ New series button has tabIndex=-1 when hidden", () => {
    render(<ComposeBar count={0} show={false} />);
    expect(screen.getByRole("button", { name: /new series/i, hidden: true }))
      .toHaveAttribute("tabindex", "-1");
  });
});
```

- [ ] Run `node_modules/.bin/vitest run test/ComposeBar.test.tsx` — expect failure.

- [ ] Implement `src/print/ui/ComposeBar.tsx`:

```tsx
// src/print/ui/ComposeBar.tsx
import { Button } from "./Button";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ComposeBarProps {
  /** Number of checked samples. */
  count: number;
  /** true → visible + interactive. */
  show?: boolean;
  onNewSeries?: () => void;
  onClear?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * ComposeBar — floating bar that surfaces when ≥1 sample is checked on the
 * contact sheet. Mirrors the CullBar positioning contract: fixed, centred,
 * bottom-7, z-50. The "samples" grain is deliberately different from the
 * CullBar's "frames" grain — never say "frames" here.
 *
 * Honest: the bar only renders usefully (aria-hidden=false, buttons tabbable)
 * when `show` is true. `show` MUST equal `count > 0` at the call site.
 */
export function ComposeBar({
  count,
  show = false,
  onNewSeries,
  onClear,
  className,
}: ComposeBarProps): JSX.Element {
  return (
    <div
      data-testid="compose-bar"
      data-show={show ? "true" : "false"}
      aria-hidden={!show || undefined}
      className={cx(
        "fixed left-1/2 bottom-20 -translate-x-1/2 z-50 flex items-center gap-1 bg-ink text-paper rounded-md shadow-lg pl-4 pr-2 py-[7px] transition-opacity",
        show ? "opacity-100" : "opacity-0 pointer-events-none",
        className,
      )}
    >
      <span className="text-sm font-semibold mr-2.5">
        <b className="font-mono">{count}</b>{" "}
        {count === 1 ? "sample" : "samples"}
      </span>
      <Button
        variant="accent"
        onClick={onNewSeries}
        {...(show ? {} : { tabIndex: -1 })}
      >
        + New series
      </Button>
      <Button
        variant="ghostInverse"
        onClick={onClear}
        {...(show ? {} : { tabIndex: -1 })}
      >
        Clear
      </Button>
    </div>
  );
}
```

- [ ] Add to `src/print/ui/index.ts`:
  ```ts
  export { ComposeBar } from "./ComposeBar";
  export type { ComposeBarProps } from "./ComposeBar";
  ```

- [ ] Run `node_modules/.bin/vitest run test/ComposeBar.test.tsx` — all pass.
- [ ] Run `npx tsc -p tsconfig.build.json --noEmit` && `npm run lint:design` — clean.
- [ ] Commit: `feat(ui): ComposeBar primitive — sample-grain compose action bar`

---

### Task 3 — Grid alignment + checkbox column in `SheetTable` / `SampleTableRow`

**Why here:** the checkbox column must keep the shared `SAMPLE_TABLE_COLS` constant as the alignment invariant. This is the most structurally load-bearing change and must be done before `SamplesPage` can wire the checkbox state. The grid-template becomes a function keyed on a boolean; both `SheetTable` header and `SampleTableRow` body call it.

**Files:** modify `src/print/components/SheetTable.tsx`, modify `src/print/components/SampleTableRow.tsx`, create `test/SampleTableRow.checkbox.test.tsx`.

- [ ] Write `test/SampleTableRow.checkbox.test.tsx`:

```tsx
// test/SampleTableRow.checkbox.test.tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { SampleTableRow } from "../src/print/components/SampleTableRow";

const BASE_PROPS = {
  name: "DOPC-1:1",
  sampleId: "#10",
  exposures: [],
  kept: 2,
  total: 3,
  tags: [],
};

describe("SampleTableRow — checkbox column", () => {
  it("renders no checkbox cell when onCheck is absent (backward compat)", () => {
    render(<SampleTableRow {...BASE_PROPS} />);
    expect(screen.queryByRole("checkbox")).toBeNull();
  });

  it("renders a checkbox cell when onCheck is provided", () => {
    render(<SampleTableRow {...BASE_PROPS} onCheck={() => {}} />);
    expect(screen.getByRole("checkbox")).toBeInTheDocument();
  });

  it("checkbox has data-checked=false when unchecked", () => {
    render(<SampleTableRow {...BASE_PROPS} checked={false} onCheck={() => {}} />);
    expect(screen.getByRole("checkbox")).toHaveAttribute("data-checked", "false");
  });

  it("checkbox has data-checked=true when checked", () => {
    render(<SampleTableRow {...BASE_PROPS} checked onCheck={() => {}} />);
    expect(screen.getByRole("checkbox")).toHaveAttribute("data-checked", "true");
  });

  it("calls onCheck when checkbox is clicked", async () => {
    const user = userEvent.setup();
    const onCheck = vi.fn();
    render(<SampleTableRow {...BASE_PROPS} onCheck={onCheck} />);
    await user.click(screen.getByRole("checkbox"));
    expect(onCheck).toHaveBeenCalledTimes(1);
  });
});
```

- [ ] Run `node_modules/.bin/vitest run test/SampleTableRow.checkbox.test.tsx` — expect failure (no checkbox rendered yet).

- [ ] Modify `src/print/components/SampleTableRow.tsx`:
  - Change `SAMPLE_TABLE_COLS` from a `const string` to a function export and keep backward compat:
    ```ts
    export const SAMPLE_TABLE_COLS_BASE =
      "minmax(244px,1.4fr) minmax(360px,2fr) 96px minmax(168px,1fr) 150px";
    /** Returns the grid-template-columns string for a row/header.
     *  `withCheckbox=true` prepends a 36px checkbox track as the left-most column.
     *  The exported constant alias keeps existing callers (SheetTable, tests) working. */
    export function sampleTableCols(withCheckbox: boolean): string {
      return withCheckbox ? `36px ${SAMPLE_TABLE_COLS_BASE}` : SAMPLE_TABLE_COLS_BASE;
    }
    // Backward-compat alias — SheetTable currently imports this as a string constant.
    export const SAMPLE_TABLE_COLS = SAMPLE_TABLE_COLS_BASE;
    ```
  - Add to `SampleTableRowProps`:
    ```ts
    /** When provided, a checkbox column is rendered as the left-most cell. */
    checked?: boolean;
    indeterminate?: boolean;
    onCheck?: () => void;
    ```
  - In the render, swap `style={{ gridTemplateColumns: SAMPLE_TABLE_COLS }}` to `style={{ gridTemplateColumns: sampleTableCols(onCheck !== undefined) }}`.
  - Prepend a checkbox cell before the Sample cell when `onCheck !== undefined`:
    ```tsx
    {onCheck !== undefined && (
      <div className="flex items-center justify-center px-1 h-[92px]">
        <Checkbox
          checked={checked ?? false}
          {...(indeterminate ? { indeterminate: true } : {})}
          onChange={onCheck}
          aria-label="Select sample"
        />
      </div>
    )}
    ```
  - Import `Checkbox` from `"../ui"`.

- [ ] Modify `src/print/components/SheetTable.tsx`:
  - Add prop `checkboxColumn?: boolean` to `SheetTableProps`.
  - Import `sampleTableCols` from `./SampleTableRow`.
  - Replace `style={{ gridTemplateColumns: SAMPLE_TABLE_COLS }}` with `style={{ gridTemplateColumns: sampleTableCols(checkboxColumn ?? false) }}`.
  - Prepend header cell when `checkboxColumn` is true:
    ```tsx
    {checkboxColumn && (
      <div className="px-1 py-2.5 flex items-center justify-center" aria-hidden="true" />
    )}
    ```
    (blank header cell — the checkbox column has no text label; `aria-hidden` since the column semantics are on each row's checkbox.)

- [ ] Run `node_modules/.bin/vitest run test/SampleTableRow.checkbox.test.tsx` — all pass.
- [ ] Run `npx tsc -p tsconfig.build.json --noEmit` && `npm run lint:design` — clean.
- [ ] Commit: `feat(sheet): checkbox column — prepend 36px track to SheetTable + SampleTableRow (grid-aligned)`

---

### Task 4 — Selection-carry helper + `SamplesPage` sample-grain picker

This task wires the sample-grain `checkedSamples: Set<number>` state into `SamplesPage`, mounts `ComposeBar`, and implements `+ New series` navigation. It also creates the typed location-state helper.

**Files:** create `src/lib/series/newSeriesNav.ts`, create `test/newSeriesNav.test.ts`, modify `src/print/pages/SamplesPage.tsx`, create `test/SamplesPage.samplePicker.test.tsx`.

- [ ] Write `test/newSeriesNav.test.ts`:

```ts
// test/newSeriesNav.test.ts
import { describe, it, expect, vi } from "vitest";
import { navigateToNewSeries, readNewSeriesSeed } from "../src/lib/series/newSeriesNav";

describe("navigateToNewSeries", () => {
  it("calls navigate with /series/new and the seed as location state", () => {
    const navigate = vi.fn();
    navigateToNewSeries(new Set([10, 11, 12]), navigate);
    expect(navigate).toHaveBeenCalledWith("/series/new", {
      state: { seedSampleIds: expect.arrayContaining([10, 11, 12]) },
    });
    // Deterministic: always an array (not a Set) so it round-trips through history.
    const call = navigate.mock.calls[0]!;
    expect(Array.isArray((call[1] as { state: { seedSampleIds: unknown } }).state.seedSampleIds)).toBe(true);
  });
});

describe("readNewSeriesSeed", () => {
  it("returns the seed array when location state is well-formed", () => {
    const loc = { state: { seedSampleIds: [10, 11] } } as unknown as Location;
    expect(readNewSeriesSeed(loc)).toEqual([10, 11]);
  });

  it("returns null when location has no state", () => {
    expect(readNewSeriesSeed({ state: null } as unknown as Location)).toBeNull();
  });

  it("returns null when seedSampleIds is not an array", () => {
    expect(readNewSeriesSeed({ state: { seedSampleIds: "bad" } } as unknown as Location)).toBeNull();
  });

  it("returns null when seedSampleIds is missing", () => {
    expect(readNewSeriesSeed({ state: {} } as unknown as Location)).toBeNull();
  });
});
```

- [ ] Run `node_modules/.bin/vitest run test/newSeriesNav.test.ts` — expect failure.

- [ ] Implement `src/lib/series/newSeriesNav.ts`:

```ts
// src/lib/series/newSeriesNav.ts
import type { NavigateFunction, Location } from "react-router-dom";

export interface NewSeriesLocationState {
  seedSampleIds: number[];
}

/**
 * Navigate to /series/new seeded with the checked sample ids. The seed is
 * passed as React Router location state — ephemeral navigational intent,
 * scoped to the history entry, not stored in Zustand or localStorage.
 */
export function navigateToNewSeries(
  sampleIds: Set<number>,
  navigate: NavigateFunction,
): void {
  navigate("/series/new", {
    state: { seedSampleIds: [...sampleIds] } satisfies NewSeriesLocationState,
  });
}

/**
 * Read the seed from the current location's state (produced by
 * `navigateToNewSeries`). Returns the array of sample ids if present and
 * well-formed, otherwise null (direct /series/new visit → full-corpus path).
 */
export function readNewSeriesSeed(location: Pick<Location, "state">): number[] | null {
  const s = location.state;
  if (s == null || typeof s !== "object") return null;
  const seed = (s as Record<string, unknown>)["seedSampleIds"];
  if (!Array.isArray(seed)) return null;
  return seed as number[];
}
```

- [ ] Write `test/SamplesPage.samplePicker.test.tsx`:

```tsx
// test/SamplesPage.samplePicker.test.tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router-dom";
import { SamplesPage } from "../src/print/pages/SamplesPage";
import * as queries from "../src/queries";

// Minimal corpus: two samples, no exposures.
const SAMPLE_A = {
  id: 10, experiment_id: 1, name: "A", display_name: null,
  notes: null, tags: [], q_units: "A-1",
};
const SAMPLE_B = {
  id: 11, experiment_id: 1, name: "B", display_name: null,
  notes: null, tags: [], q_units: "A-1",
};

function setup() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  vi.spyOn(queries, "useCorpusSamples").mockReturnValue({
    data: [SAMPLE_A, SAMPLE_B], isLoading: false, isError: false,
  } as ReturnType<typeof queries.useCorpusSamples>);
  vi.spyOn(queries, "useCorpusExposures").mockReturnValue(
    { byId: new Map(), isLoading: false });
  vi.spyOn(queries, "useScreenedProgress").mockReturnValue({ screened: 0, total: 2 });
  vi.spyOn(queries, "useExperiments").mockReturnValue(
    { data: [], isLoading: false } as unknown as ReturnType<typeof queries.useExperiments>);
  vi.spyOn(queries, "useSetExposureStatusBatch").mockReturnValue(
    { mutate: vi.fn() } as unknown as ReturnType<typeof queries.useSetExposureStatusBatch>);

  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter>
        <SamplesPage />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("SamplesPage — sample-grain picker", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("compose bar is hidden when no samples are checked", () => {
    setup();
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
  });

  it("compose bar becomes visible when a sample is checked", async () => {
    const user = userEvent.setup();
    setup();
    const checkboxes = screen.getAllByRole("checkbox");
    await user.click(checkboxes[0]!);
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");
    expect(screen.getByTestId("compose-bar")).toHaveTextContent("1");
  });

  it("Clear resets selection and hides the compose bar", async () => {
    const user = userEvent.setup();
    setup();
    await user.click(screen.getAllByRole("checkbox")[0]!);
    await user.click(screen.getByRole("button", { name: /clear/i, hidden: false }));
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
  });

  it("cull bar remains independent of sample-grain selection", async () => {
    const user = userEvent.setup();
    setup();
    await user.click(screen.getAllByRole("checkbox")[0]!);
    // Frame-level CullBar should still be hidden (no exposure selected).
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
  });
});
```

- [ ] Run `node_modules/.bin/vitest run test/SamplesPage.samplePicker.test.tsx` — expect failure.

- [ ] Modify `src/print/pages/SamplesPage.tsx`:
  - Add import `ComposeBar` from `"../ui"` and `navigateToNewSeries` from `"../../lib/series/newSeriesNav"`.
  - Add state: `const [checkedSamples, setCheckedSamples] = useState<Set<number>>(() => new Set());`
  - Add `function toggleSampleCheck(sampleId: number): void { setCheckedSamples(prev => { const next = new Set(prev); if (next.has(sampleId)) next.delete(sampleId); else next.add(sampleId); return next; }); }`
  - Add `checkboxColumn` and `checked` / `onCheck` to each `SampleTableRow` and `SheetTable`:
    - `<SheetTable checkboxColumn ...>`
    - `<SampleTableRow ... checked={checkedSamples.has(s.id)} onCheck={() => toggleSampleCheck(s.id)} ...>`
  - Mount `ComposeBar` alongside `CullBar` (bottom of return, inside the `PageFrame`):
    ```tsx
    <ComposeBar
      count={checkedSamples.size}
      show={checkedSamples.size > 0}
      onNewSeries={() => navigateToNewSeries(checkedSamples, navigate)}
      onClear={() => setCheckedSamples(new Set())}
    />
    ```
  - `ComposeBar` sits at `bottom-20` (above `CullBar`'s `bottom-7`) so they don't stack when both are visible.

- [ ] Run `node_modules/.bin/vitest run test/SamplesPage.samplePicker.test.tsx` — all pass.
- [ ] Run `node_modules/.bin/vitest run test/newSeriesNav.test.ts` — all pass.
- [ ] Run `npx tsc -p tsconfig.build.json --noEmit` && `npm run lint:design` — clean.
- [ ] Commit: `feat(samples): sample-grain checkbox picker + ComposeBar + navigateToNewSeries`

---

### Task 5 — Scoping seed: filter picker by seeded sample ids

`SeriesScopingPage` must read the location state, filter `pickerQ.data` to only the seeded samples when a seed is present, and pass that filtered slice to `proposeOrdering`. No cold-corpus path yet — that is Task 6.

**Files:** modify `src/print/pages/SeriesScopingPage.tsx`, create `test/scopingDerive.seedFilter.test.ts`.

- [ ] Write `test/scopingDerive.seedFilter.test.ts`:

```ts
// test/scopingDerive.seedFilter.test.ts
import { describe, it, expect } from "vitest";
import { filterPickerBySeed } from "../src/print/pages/scopingDerive";
import type { PickerSampleRow } from "../src/api";

function row(id: number): PickerSampleRow {
  return {
    sample: { id, experiment_id: 1, name: `s${id}`, display_name: null, notes: null, tags: [] },
    indexing_exposure_id: null,
    all_exposures: [],
  };
}

describe("filterPickerBySeed", () => {
  it("returns all rows when seed is null (full-corpus path)", () => {
    const rows = [row(1), row(2), row(3)];
    expect(filterPickerBySeed(rows, null)).toEqual(rows);
  });

  it("filters to only the seeded ids when seed is provided", () => {
    const rows = [row(1), row(2), row(3)];
    const filtered = filterPickerBySeed(rows, [1, 3]);
    expect(filtered.map((r) => r.sample.id)).toEqual([1, 3]);
  });

  it("returns empty array when no seeded ids match (honest empty path)", () => {
    expect(filterPickerBySeed([row(1), row(2)], [99])).toEqual([]);
  });
});
```

- [ ] Run targeted test — expect failure.

- [ ] Add `filterPickerBySeed` to `src/print/pages/scopingDerive.ts`:

```ts
/** Filter the corpus picker rows to only those whose sample.id is in `seed`.
 *  When `seed` is null the full corpus is returned (no-filter path for a
 *  direct /series/new visit). */
export function filterPickerBySeed(
  rows: import("../../api").PickerSampleRow[],
  seed: number[] | null,
): import("../../api").PickerSampleRow[] {
  if (seed === null) return rows;
  const ids = new Set(seed);
  return rows.filter((r) => ids.has(r.sample.id));
}
```

- [ ] Modify `src/print/pages/SeriesScopingPage.tsx`:
  - Add imports: `useLocation` from `react-router-dom`, `readNewSeriesSeed` from `../../lib/series/newSeriesNav`, `filterPickerBySeed` from `./scopingDerive`.
  - Inside the component body (before any `useMemo`):
    ```ts
    const location = useLocation();
    const seed = useMemo(() => readNewSeriesSeed(location), [location]);
    ```
  - Change the `proposal` `useMemo` to filter the picker first:
    ```ts
    const seededPicker = useMemo(
      () => filterPickerBySeed(pickerQ.data ?? [], seed),
      [pickerQ.data, seed],
    );
    const proposal = useMemo(
      () => proposeOrdering(tagsQ.data ?? [], seededPicker),
      [tagsQ.data, seededPicker],
    );
    ```
  - `allSampleIds` continues to derive from `rows` / `loose` (already filtered via `proposal`). No other changes needed for this task.

- [ ] Run `node_modules/.bin/vitest run test/scopingDerive.seedFilter.test.ts` — pass.
- [ ] Run `npx tsc -p tsconfig.build.json --noEmit` && `npm run lint:design` — clean.
- [ ] Commit: `feat(scoping): read location-state seed + filter picker to seeded samples`

---

### Task 6 — Cold-corpus assign path in scoping

When `proposal.orderingKey === undefined` and the user arrived with a seed (i.e., they deliberately chose samples but none have a shared tag), the current dead-end `EmptyState` is replaced with:

1. An `Input` for naming the ordering variable key.
2. A `ColdAssignPanel` showing one row per seeded sample with an `Input` for the per-sample value.
3. The standard foot state gate + "Confirm & build" enabled only when key is non-empty and all samples have a value.

**Files:** create `src/print/components/ColdAssignPanel.tsx`, add helpers to `src/print/pages/scopingDerive.ts`, create `test/scopingDerive.coldPath.test.ts`, modify `src/print/pages/SeriesScopingPage.tsx`.

- [ ] Write `test/scopingDerive.coldPath.test.ts`:

```ts
// test/scopingDerive.coldPath.test.ts
import { describe, it, expect } from "vitest";
import {
  buildColdAssignRows, canColdBuild, buildColdScopePayload,
} from "../src/print/pages/scopingDerive";

const SEED = [
  { sampleId: 10, sampleName: "A" },
  { sampleId: 11, sampleName: "B" },
];

describe("buildColdAssignRows", () => {
  it("creates one row per seeded sample, all values empty", () => {
    const rows = buildColdAssignRows(SEED);
    expect(rows).toHaveLength(2);
    expect(rows[0]).toMatchObject({ sampleId: 10, sampleName: "A", value: "" });
    expect(rows[1]).toMatchObject({ sampleId: 11, sampleName: "B", value: "" });
  });
});

describe("canColdBuild", () => {
  it("returns false when key is empty", () => {
    expect(canColdBuild("", [{ sampleId: 10, sampleName: "A", value: "1.0" }])).toBe(false);
  });

  it("returns false when any sample has an empty value", () => {
    const rows = [
      { sampleId: 10, sampleName: "A", value: "1.0" },
      { sampleId: 11, sampleName: "B", value: "" },
    ];
    expect(canColdBuild("ratio", rows)).toBe(false);
  });

  it("returns true when key and all values are non-empty", () => {
    const rows = [
      { sampleId: 10, sampleName: "A", value: "1.0" },
      { sampleId: 11, sampleName: "B", value: "2.0" },
    ];
    expect(canColdBuild("ratio", rows)).toBe(true);
  });
});

describe("buildColdScopePayload", () => {
  it("returns the (sampleId, value) pairs for all rows", () => {
    const rows = [
      { sampleId: 10, sampleName: "A", value: "1.0" },
      { sampleId: 11, sampleName: "B", value: "2.0" },
    ];
    expect(buildColdScopePayload(rows)).toEqual([
      { sampleId: 10, value: "1.0" },
      { sampleId: 11, value: "2.0" },
    ]);
  });
});
```

- [ ] Run targeted test — expect failure.

- [ ] Add to `src/print/pages/scopingDerive.ts`:

```ts
export interface ColdAssignRow {
  sampleId: number;
  sampleName: string;
  value: string;
}

/** Seed the cold-assign worksheet from a list of (id, name) pairs.
 *  All values start empty — the user fills them in. */
export function buildColdAssignRows(
  seed: ReadonlyArray<{ sampleId: number; sampleName: string }>,
): ColdAssignRow[] {
  return seed.map((s) => ({ sampleId: s.sampleId, sampleName: s.sampleName, value: "" }));
}

/** Cold-assign build gate: key must be non-empty and every sample must have a
 *  non-empty value. Never blocks on a partial fill — controls-don't-lie. */
export function canColdBuild(key: string, rows: ColdAssignRow[]): boolean {
  if (key.trim() === "") return false;
  return rows.every((r) => r.value.trim() !== "");
}

/** Cold-assign scope payload: all rows unconditionally (gate is upstream). */
export function buildColdScopePayload(rows: ColdAssignRow[]): { sampleId: number; value: string }[] {
  return rows.map((r) => ({ sampleId: r.sampleId, value: r.value }));
}
```

- [ ] Implement `src/print/components/ColdAssignPanel.tsx`:

```tsx
// src/print/components/ColdAssignPanel.tsx
import { Input, Kicker, HintText, Field } from "../ui";
import type { ColdAssignRow } from "../pages/scopingDerive";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ColdAssignPanelProps {
  rows: ColdAssignRow[];
  /** The ordering-variable key the user is naming. */
  variableKey: string;
  onKeyChange: (key: string) => void;
  onValueChange: (sampleId: number, value: string) => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * ColdAssignPanel — the cold-corpus path inside SeriesScopingPage.
 *
 * When no tag key can be proposed (cold / unstructured selection), this panel
 * lets the user name the ordering variable once and fill in each sample's value.
 * Composed from `src/print/ui/` primitives only; no legacy imports.
 *
 * Layout: a labelled key-name Input at the top, then one row per sample showing
 * the sample name (text-body font-semibold text-ink) and a value Input.
 * controls-don't-lie: values are plain inputs — no fake "type-and-it-saves"
 * affordance; the foot-gate "Confirm & build" is the single commit action.
 */
export function ColdAssignPanel({
  rows,
  variableKey,
  onKeyChange,
  onValueChange,
  className,
}: ColdAssignPanelProps): JSX.Element {
  return (
    <div data-testid="cold-assign-panel" className={cx("space-y-4", className)}>
      <div>
        <Kicker tone="accent" className="mb-1">
          Name the ordering variable
        </Kicker>
        <p className="text-body text-ink-soft mb-2">
          These samples share no tag key yet. Name the variable (e.g. "lipid ratio") and assign each
          sample's value — then Confirm &amp; build.
        </p>
        <Input
          value={variableKey}
          onValueChange={onKeyChange}
          placeholder="e.g. lipid ratio"
          aria-label="Ordering variable name"
          data-testid="cold-key-input"
        />
      </div>
      <div className="space-y-1">
        {rows.map((r) => (
          <div
            key={r.sampleId}
            data-testid="cold-assign-row"
            data-sample-id={r.sampleId}
            className="flex items-center gap-3 rounded border border-hair px-3 py-2"
          >
            <div className="flex-1 min-w-0">
              <div className="text-body font-semibold text-ink truncate">{r.sampleName}</div>
              <HintText className="font-mono">smp_{r.sampleId}</HintText>
            </div>
            <Input
              value={r.value}
              onValueChange={(v) => onValueChange(r.sampleId, v)}
              placeholder="value"
              aria-label={`Value for ${r.sampleName}`}
              inputSize="sm"
              className="w-28 flex-shrink-0"
            />
          </div>
        ))}
      </div>
    </div>
  );
}
```

- [ ] Modify `src/print/pages/SeriesScopingPage.tsx` — add cold-corpus state + branching:
  - After the existing `const [history, setHistory] = ...` block, add:
    ```ts
    // ── Cold-corpus path (no proposable variable, user arrived with a seed) ──
    const [coldKey, setColdKey] = useState("");
    const [coldRows, setColdRows] = useState<ColdAssignRow[]>([]);
    // Seed cold rows once when the page loads into the cold path.
    useEffect(() => {
      if (!isLoading && proposal.orderingKey === undefined && seed !== null && rows.length === 0) {
        const names = new Map<number, string>();
        for (const s of pickerQ.data ?? []) names.set(s.sample.id, s.sample.display_name ?? s.sample.name ?? "");
        setColdRows(buildColdAssignRows(
          (seed ?? []).map((id) => ({ sampleId: id, sampleName: names.get(id) ?? `smp_${id}` })),
        ));
        setColdKey("");
      }
    }, [isLoading, proposal.orderingKey, seed, rows.length, pickerQ.data]);
    ```
  - Import `ColdAssignRow`, `buildColdAssignRows`, `canColdBuild`, `buildColdScopePayload` from `./scopingDerive`.
  - Import `ColdAssignPanel` from `../components/ColdAssignPanel`.
  - In the `proposal.orderingKey === undefined` branch inside `<Skeleton>`, replace the full dead-end `EmptyState` with:
    ```tsx
    {seed !== null ? (
      // Cold path: user arrived from the contact sheet with a selection but
      // no proposable variable. Let them name the variable + assign values.
      <ColdAssignPanel
        rows={coldRows}
        variableKey={coldKey}
        onKeyChange={setColdKey}
        onValueChange={(id, v) =>
          setColdRows((cur) => cur.map((r) => r.sampleId === id ? { ...r, value: v } : r))
        }
      />
    ) : (
      <EmptyState
        title={rows.length === 0 && loose.length === 0 ? "Nothing to scope yet" : "No shared ordering variable"}
        body={/* existing body */}
      />
    )}
    ```
  - The foot state + build gate for the cold path: add alongside the existing `canBuild` derivation:
    ```ts
    const canColdBuildNow = canColdBuild(coldKey, coldRows);
    const isColdPath = proposal.orderingKey === undefined && seed !== null;
    ```
  - In `handleBuild`, branch on `isColdPath`:
    ```ts
    const handleBuild = (): void => {
      if (isColdPath) {
        if (!canColdBuildNow) return;
        pendingBuildRef.current = true;
        scopeSeries.mutate({ key: coldKey.trim(), tags: buildColdScopePayload(coldRows) });
      } else {
        if (proposal.orderingKey === undefined) return;
        pendingBuildRef.current = true;
        scopeSeries.mutate({ key: proposal.orderingKey, tags: buildScopePayload(rows) });
      }
    };
    ```
  - The cold path's ScopePlate / foot state: when `isColdPath`, pass `buildDisabled={!canColdBuildNow}` and a foot-note describing the cold state. The existing `ScopePlate` rendering wraps around both paths — for the cold path, render only the `ColdAssignPanel` inside a minimal Card with the foot zone, not the full `ScopePlate` (which requires a `seriesName` + `orderedBy` + `rows`). Use a simpler layout matching the existing `EmptyState` wrapper width.

- [ ] Run `node_modules/.bin/vitest run test/scopingDerive.coldPath.test.ts` — all pass.
- [ ] Run `npx tsc -p tsconfig.build.json --noEmit` && `npm run lint:design` — clean.
- [ ] Commit: `feat(scoping): cold-corpus assign path — ColdAssignPanel + buildColdAssignRows/canColdBuild`

---

### Task 7 — `scopeAndCreateSeries` mutator + `Confirm & build` creates the series

Today `handleBuild` writes tags and navigates to `/series` (the folio). After this task, the successful tag write is immediately followed by a series create (`POST /api/series`), and the page navigates to `/series/:newId`. This requires a new mutator that sequences the two writes, or a two-step approach in the page. Architecture decision: **a new mutator `scopeAndCreateSeries`** that runs both API calls in its `request` function. This keeps the page logic minimal and the mutation queue contract clean (one `client_op_id` covers the atomic intent: scope-then-create). The mutator's `onSuccess` splices the new series into the listing cache.

**Files:** create `src/lib/queue/mutators/scopeAndCreateSeries.ts`, create `test/queue/scopeAndCreateSeries.contract.test.ts`, modify `src/queries.ts`, modify `src/print/pages/SeriesScopingPage.tsx`.

- [ ] Write `test/queue/scopeAndCreateSeries.contract.test.ts`:

```ts
// test/queue/scopeAndCreateSeries.contract.test.ts
import { describe, it, expect, vi, beforeEach } from "vitest";
import { QueryClient } from "@tanstack/react-query";
import { scopeAndCreateSeriesMutator } from "../../src/lib/queue/mutators/scopeAndCreateSeries";
import * as api from "../../src/api";

const FULL_SERIES: api.Series = {
  id: 7, title: "Series by ratio", description: null,
  state: "draft", ordering_variable: "ratio", order_rule: "asc",
  forked_from_id: null, forked_at_hash: null, content_hash: "abc",
  view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
  members: [], samples: [],
};

describe("scopeAndCreateSeriesMutator", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("kind is series_save (reuses existing series SSE event kind)", () => {
    expect(scopeAndCreateSeriesMutator.kind).toBe("series_save");
  });

  it("onMutate is a no-op", () => {
    const qc = new QueryClient();
    const ctx = scopeAndCreateSeriesMutator.onMutate(
      { key: "ratio", tags: [], title: "S", samples: [],
        username: "a", clientId: "c", clientOpId: "op1" } as never,
      qc,
    );
    expect(typeof ctx.restore).toBe("function");
  });

  it("request: calls batchSampleTags then saveSeries (create, no id)", async () => {
    const batchSpy = vi.spyOn(api, "batchSampleTags").mockResolvedValue([]);
    const saveSpy = vi.spyOn(api, "saveSeries").mockResolvedValue(FULL_SERIES);
    await scopeAndCreateSeriesMutator.request(
      { key: "ratio", tags: [{ sampleId: 10, value: "1:1" }],
        title: "Series by ratio",
        samples: [{ sample_id: 10 }],
        orderingVariable: "ratio",
        username: "a", clientId: "c", clientOpId: "op1" } as never,
      new AbortController().signal,
    );
    expect(batchSpy).toHaveBeenCalledWith(
      "ratio",
      [{ sample_id: 10, value: "1:1" }],
      "scoping",
      expect.anything(),
    );
    expect(saveSpy).toHaveBeenCalledWith(
      expect.objectContaining({ title: "Series by ratio", ordering_variable: "ratio" }),
      undefined, // no id → create
      expect.anything(),
    );
  });

  it("onSuccess: splices the new series into the series-list cache and invalidates scoping caches", () => {
    const qc = new QueryClient();
    qc.setQueryData(["series-list"], []);
    const inv = vi.spyOn(qc, "invalidateQueries");
    scopeAndCreateSeriesMutator.onSuccess(
      { key: "ratio", tags: [], title: "S", samples: [],
        username: "a", clientId: "c", clientOpId: "op1" } as never,
      FULL_SERIES,
      qc,
    );
    // Splices the full series into the series-list.
    expect(qc.getQueryData(["series-list"])).toEqual([FULL_SERIES]);
    // Invalidates the scoping proposal sources.
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-sample-tags"] });
    expect(inv).toHaveBeenCalledWith({ queryKey: ["corpus-picker-samples"] });
  });
});
```

- [ ] Run targeted test — expect failure.

- [ ] Implement `src/lib/queue/mutators/scopeAndCreateSeries.ts`:

```ts
// src/lib/queue/mutators/scopeAndCreateSeries.ts
import type { Mutator, RollbackContext } from "../types";
import type { AuthOpts, Series, SeriesSampleInput } from "../../../api";
import { authOpts } from "../../authOpts";
import { queryKeys } from "../../../queries";
import * as api from "../../../api";

export interface ScopeAndCreateSeriesInput {
  /** The ordering-variable key to write to sample_tags (source='scoping'). */
  key: string;
  tags: { sampleId: number; value: string }[];
  /** SaveSeriesBody fields for the new series. */
  title: string;
  description?: string | null;
  samples: SeriesSampleInput[];
  orderingVariable?: string | null;
}

interface ScopeAndCreateSeriesScope {
  username: string | undefined;
  clientId: string;
}

function buildAuthOpts(p: { username: string | undefined; clientId: string; clientOpId: string }): AuthOpts {
  return authOpts(p.username, p.clientId, p.clientOpId);
}

/**
 * scopeAndCreateSeriesMutator — atomic scope-then-create.
 *
 * `request` runs sequentially:
 *   1. batchSampleTags (source='scoping') — writes the ordering tags.
 *   2. saveSeries (no id → POST /api/series) — creates the series.
 *
 * Uses kind='series_save' so the SSE frame from the create (which the backend
 * emits as a `series_save` event) routes through the existing
 * `saveSeriesMutator.synthesizeFromSse` path for foreign-tab replay. The
 * scoping batch tags emit `add_tag` frames; those replay via the existing
 * `addSampleTagMutator.synthesizeFromSse` path (no change needed).
 *
 * onMutate is a no-op (mirrors scopeSeriesMutator; the user navigates away
 * on success, so there is no on-surface cache to optimistically patch).
 *
 * onSuccess: splices the new Series into the listing cache (so the folio
 * sees it immediately on the navigation to /series/:id) and invalidates the
 * two scoping caches.
 */
export const scopeAndCreateSeriesMutator: Mutator<
  ScopeAndCreateSeriesInput,
  ScopeAndCreateSeriesScope,
  Series
> = {
  kind: "series_save",
  onMutate: (): RollbackContext => ({ restore: () => {} }),
  request: async (p, _signal) => {
    const opts = buildAuthOpts(p);
    await api.batchSampleTags(
      p.key,
      p.tags.map((t) => ({ sample_id: t.sampleId, value: t.value })),
      "scoping",
      opts,
    );
    return api.saveSeries(
      {
        title: p.title,
        samples: p.samples,
        ...(p.description !== undefined ? { description: p.description } : {}),
        ...(p.orderingVariable !== undefined ? { ordering_variable: p.orderingVariable } : {}),
      },
      undefined, // no id → POST /api/series (create)
      opts,
    );
  },
  onSuccess: (_p, response, qc) => {
    // Splice new series into the listing.
    const existing = qc.getQueryData<Series[]>(queryKeys.seriesList) ?? [];
    if (api.isFullSeries(response)) {
      qc.setQueryData(queryKeys.seriesList, [...existing, response]);
      qc.setQueryData(queryKeys.series(response.id), response);
    } else {
      qc.invalidateQueries({ queryKey: queryKeys.seriesList });
    }
    // Invalidate scoping proposal sources (same as scopeSeriesMutator).
    qc.invalidateQueries({ queryKey: queryKeys.corpusSampleTags });
    qc.invalidateQueries({ queryKey: queryKeys.corpusPickerSamples });
  },
  synthesizeFromSse: (remote, base) => {
    const payload = (remote.payload as Record<string, unknown> | undefined) ?? {};
    return { ...base, ...payload, id: remote.entity_id } as unknown as Series;
  },
};
```

- [ ] Add `useScopeAndCreateSeries` to `src/queries.ts`:

```ts
export function useScopeAndCreateSeries() {
  const username = useAppState((s) => s.username);
  return useQueueMutation(scopeAndCreateSeriesMutator, { username, clientId: CLIENT_ID });
}
```

And add the import at the top of `queries.ts`:
```ts
import { scopeAndCreateSeriesMutator } from "./lib/queue/mutators/scopeAndCreateSeries";
```

- [ ] Modify `src/print/pages/SeriesScopingPage.tsx` — replace `useScopeSeries` with `useScopeAndCreateSeries` and thread the new series ID to navigation:
  - Replace `const scopeSeries = useScopeSeries();` with `const scopeAndCreate = useScopeAndCreateSeries();`.
  - The pending ref + success effect now reads `scopeAndCreate`:
    ```ts
    const pendingBuildRef = useRef(false);
    const handleBuild = (): void => {
      if (isColdPath) {
        if (!canColdBuildNow) return;
        pendingBuildRef.current = true;
        scopeAndCreate.mutate({
          key: coldKey.trim(),
          tags: buildColdScopePayload(coldRows),
          title: `Series by ${coldKey.trim()}`,
          samples: coldRows.map((r, i) => ({ sample_id: r.sampleId, position: i })),
          orderingVariable: coldKey.trim(),
        });
      } else {
        if (proposal.orderingKey === undefined) return;
        pendingBuildRef.current = true;
        scopeAndCreate.mutate({
          key: proposal.orderingKey,
          tags: buildScopePayload(rows),
          title: `Series by ${keyLabel}`,
          samples: sorted.map((r, i) => ({ sample_id: r.sampleId, position: i })),
          orderingVariable: proposal.orderingKey,
        });
      }
    };
    useEffect(() => {
      if (!scopeAndCreate.isSuccess || !pendingBuildRef.current) return;
      pendingBuildRef.current = false;
      // Navigate to the new series builder. scopeAndCreate.data is a Series.
      const newId = scopeAndCreate.data?.id;
      navigate(newId !== undefined ? `/series/${newId}` : "/series");
    }, [scopeAndCreate.isSuccess, scopeAndCreate.data, navigate]);
    useEffect(() => {
      if (scopeAndCreate.error) pendingBuildRef.current = false;
    }, [scopeAndCreate.error]);
    ```
  - Update `scopeSeries.error` references → `scopeAndCreate.error`.
  - Pass `buildDisabled` correctly for both paths:
    ```ts
    const canBuildForPath = isColdPath ? canColdBuildNow : canScopeBuild(rows, proposal.orderingKey);
    ```
  - Pass `{...(!canBuildForPath ? { buildDisabled: true } : {})}` to `ScopePlate`.

- [ ] Run `node_modules/.bin/vitest run test/queue/scopeAndCreateSeries.contract.test.ts` — all pass.
- [ ] Run `npx tsc -p tsconfig.build.json --noEmit` && `npm run lint:design` — clean.
- [ ] Commit: `feat(scoping): scopeAndCreateSeries mutator — sequential tag-write + series create, navigate to /series/:id`

---

### Task 8 — Playwright E2E: full new-series create path

One E2E spec covering the complete B-shaped flow: contact-sheet checkbox → compose bar → `+ New series` → scoping page pre-seeded → confirm-and-build → tags + series posted → navigate to `/series/:id`.

**Files:** create `e2e/new-series-creation.spec.ts`.

- [ ] Implement `e2e/new-series-creation.spec.ts`:

```ts
// e2e/new-series-creation.spec.ts
import { test, expect, type Page } from "@playwright/test";

const EXPERIMENT = {
  id: 1, name: "SSRL Test", path: "/p", data_dir: "/p/data",
  analysis_dir: "/p/analysis", manifest_path: null, created_at: "2026-05-01",
};
const SAMPLES = [
  { id: 10, experiment_id: 1, display_name: "DOPC-1:1", name: "run10",
    notes: null, tags: [{ id: 1, key: "ratio", value: "1:1", source: "manual" }], q_units: "A-1", phase: null },
  { id: 11, experiment_id: 1, display_name: "DOPC-2:1", name: "run11",
    notes: null, tags: [{ id: 2, key: "ratio", value: "2:1", source: "manual" }], q_units: "A-1", phase: null },
];
const EXPOSURES_10 = [
  { id: 1, sample_id: 10, filename: "f1.dat", kind: "file", selected: true,
    status: "accepted", image_path: null, image_version: "", tags: [], sources: [],
    trace_hash: null, analysis_inputs_hash: null },
];
const EXPOSURES_11 = [
  { id: 2, sample_id: 11, filename: "f2.dat", kind: "file", selected: true,
    status: "accepted", image_path: null, image_version: "", tags: [], sources: [],
    trace_hash: null, analysis_inputs_hash: null },
];
const NEW_SERIES = {
  id: 42, title: "Series by ratio", description: null, state: "draft",
  ordering_variable: "ratio", order_rule: "asc", forked_from_id: null,
  forked_at_hash: null, content_hash: "xyz",
  view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
  members: [], samples: [],
};

async function mockAll(
  page: Page,
  captured: { tagsBatch: unknown; seriesCreate: unknown },
): Promise<void> {
  await page.addInitScript(() =>
    localStorage.setItem(
      "himalaya-ui:state",
      JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 5 }),
    ),
  );
  await page.route("**/api/users", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  await page.route("**/api/experiments", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([EXPERIMENT]) }));
  await page.route("**/api/samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(SAMPLES) }));
  await page.route("**/api/samples/10/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(EXPOSURES_10) }));
  await page.route("**/api/samples/11/exposures*", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(EXPOSURES_11) }));
  // Scoping reads:
  await page.route("**/api/sample-tags", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([
        { key: "ratio", value: "1:1" },
        { key: "ratio", value: "2:1" },
      ]) }));
  await page.route("**/api/picker-samples", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify([
        { sample: SAMPLES[0], indexing_exposure_id: null, all_exposures: [] },
        { sample: SAMPLES[1], indexing_exposure_id: null, all_exposures: [] },
      ]) }));
  // Batch write — capture body.
  await page.route("**/api/samples/tags/batch", async (route) => {
    captured.tagsBatch = route.request().postDataJSON();
    await route.fulfill({ status: 201, contentType: "application/json",
      body: JSON.stringify([
        { id: 1, sample_id: 10, key: "ratio", value: "1:1", source: "scoping" },
        { id: 2, sample_id: 11, key: "ratio", value: "2:1", source: "scoping" },
      ]) });
  });
  // Series create — capture body, return the new series.
  await page.route("**/api/series", async (route) => {
    if (route.request().method() === "POST") {
      captured.seriesCreate = route.request().postDataJSON();
      await route.fulfill({ status: 201, contentType: "application/json",
        body: JSON.stringify(NEW_SERIES) });
    } else {
      await route.fulfill({ status: 200, contentType: "application/json", body: "[]" });
    }
  });
  // Builder reads (after navigation to /series/42):
  await page.route("**/api/series/42", (r) =>
    r.fulfill({ status: 200, contentType: "application/json",
      body: JSON.stringify(NEW_SERIES) }));
  await page.route("**/api/series/42/forks", (r) =>
    r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
  // SSE drain:
  await page.route("**/api/events", (r) =>
    r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));
}

test.describe("new-series creation flow (M-A)", () => {
  test("full B-path: check samples → compose bar → scoping → confirm → /series/:id", async ({ page }) => {
    const captured = { tagsBatch: null as unknown, seriesCreate: null as unknown };
    await mockAll(page, captured);
    await page.goto("/samples");

    // The contact sheet must load.
    await expect(page.getByTestId("samples-page")).toBeVisible();
    await expect(page.getByTestId("sample-table-row")).toHaveCount(2);

    // Compose bar is hidden initially.
    await expect(page.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");

    // Check the first sample row's checkbox.
    const firstRowCheckbox = page.getByTestId("sample-table-row").first().getByRole("checkbox");
    await firstRowCheckbox.click();
    await expect(page.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");
    await expect(page.getByTestId("compose-bar")).toContainText("1 sample");

    // Check the second sample.
    const secondRowCheckbox = page.getByTestId("sample-table-row").nth(1).getByRole("checkbox");
    await secondRowCheckbox.click();
    await expect(page.getByTestId("compose-bar")).toContainText("2 samples");

    // Cull bar must remain hidden (no exposure selected).
    await expect(page.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");

    // Click + New series.
    await page.getByRole("button", { name: /new series/i }).click();

    // Scoping page loads pre-seeded (only the 2 selected samples, not the full corpus).
    await expect(page).toHaveURL(/\/series\/new/);
    await expect(page.getByTestId("scoping-page")).toBeVisible();
    // Both seeded samples appear as scoping member rows.
    await expect(page.getByTestId("scope-sample-row")).toHaveCount(2);

    // Confirm & build is enabled.
    const buildBtn = page.getByRole("button", { name: /confirm & build/i });
    await expect(buildBtn).toBeEnabled();
    await buildBtn.click();

    // Tags batch posted with key="ratio".
    await expect.poll(
      () => (captured.tagsBatch as Record<string, unknown> | null)?.["key"],
      { timeout: 4000 },
    ).toBe("ratio");

    // Series POST was made.
    await expect.poll(
      () => (captured.seriesCreate as Record<string, unknown> | null)?.["title"],
      { timeout: 4000 },
    ).toMatch(/ratio/i);

    // Navigated to the new series builder.
    await expect(page).toHaveURL(/\/series\/42/);
  });

  test("Clear on compose bar resets selection", async ({ page }) => {
    const captured = { tagsBatch: null as unknown, seriesCreate: null as unknown };
    await mockAll(page, captured);
    await page.goto("/samples");
    await expect(page.getByTestId("sample-table-row")).toHaveCount(2);

    await page.getByTestId("sample-table-row").first().getByRole("checkbox").click();
    await expect(page.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");

    await page.getByTestId("compose-bar").getByRole("button", { name: /clear/i }).click();
    await expect(page.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
  });

  test("direct visit to /series/new shows full corpus (no seed, existing behaviour preserved)", async ({ page }) => {
    const captured = { tagsBatch: null as unknown, seriesCreate: null as unknown };
    await mockAll(page, captured);
    // Visit without going through the checker.
    await page.goto("/series/new");
    await expect(page.getByTestId("scoping-page")).toBeVisible();
    // Full corpus: both samples appear (the mock returns both in picker-samples).
    await expect(page.getByTestId("scope-sample-row")).toHaveCount(2);
  });

  test("cold-corpus path: no proposable variable → assign panel appears", async ({ page }) => {
    // Override to a corpus with no tags.
    await page.addInitScript(() =>
      localStorage.setItem(
        "himalaya-ui:state",
        JSON.stringify({ state: { username: "alice", tutorialSeen: true }, version: 5 }),
      ),
    );
    await page.route("**/api/users", (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
    await page.route("**/api/experiments", (r) =>
      r.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify([EXPERIMENT]) }));
    await page.route("**/api/samples", (r) =>
      r.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(SAMPLES) }));
    await page.route("**/api/samples/10/exposures*", (r) =>
      r.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(EXPOSURES_10) }));
    await page.route("**/api/samples/11/exposures*", (r) =>
      r.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify(EXPOSURES_11) }));
    // No tags in corpus — proposeOrdering returns orderingKey=undefined.
    await page.route("**/api/sample-tags", (r) =>
      r.fulfill({ status: 200, contentType: "application/json", body: "[]" }));
    // Picker: both samples, no tags.
    await page.route("**/api/picker-samples", (r) =>
      r.fulfill({ status: 200, contentType: "application/json",
        body: JSON.stringify([
          { sample: { ...SAMPLES[0], tags: [] }, indexing_exposure_id: null, all_exposures: [] },
          { sample: { ...SAMPLES[1], tags: [] }, indexing_exposure_id: null, all_exposures: [] },
        ]) }));
    await page.route("**/api/events", (r) =>
      r.fulfill({ status: 200, contentType: "text/event-stream", body: "" }));

    // Arrive via the picker (navigate with state carrying both sample ids).
    // Simulate by going to /samples, checking, and clicking + New series.
    await page.goto("/samples");
    await expect(page.getByTestId("sample-table-row")).toHaveCount(2);
    await page.getByTestId("sample-table-row").first().getByRole("checkbox").click();
    await page.getByTestId("sample-table-row").nth(1).getByRole("checkbox").click();
    await page.getByRole("button", { name: /new series/i }).click();

    // Cold assign panel must appear (no ordering variable proposable).
    await expect(page.getByTestId("cold-assign-panel")).toBeVisible();
    // Confirm & build must be disabled until key + values filled.
    await expect(page.getByRole("button", { name: /confirm & build/i })).toBeDisabled();
  });
});
```

- [ ] Run `npm run e2e -- --grep "M-A"` (or by spec file: `npx playwright test e2e/new-series-creation.spec.ts`) — all pass.
- [ ] Run full `npm run e2e` — no regressions (existing corpus-culling + series-scoping specs remain green).
- [ ] Run `npx tsc -p tsconfig.build.json --noEmit` && `npm run lint:design` — clean.
- [ ] Commit: `test(e2e): new-series creation flow M-A — full B-path Playwright spec`

---

## Self-Review: Spec Coverage Checklist

| Spec requirement | Covered | Where |
|---|---|---|
| Left-most checkbox column on contact sheet, sample-grain | Yes | Task 3 (`SheetTable` + `SampleTableRow`), Task 4 (`SamplesPage`) |
| Selection distinct from frame-grain CullBar | Yes | Task 4 — separate `checkedSamples` state; CullBar unchanged |
| Compose bar shows only when ≥1 sample checked | Yes | Task 2 (`ComposeBar` `show` prop), Task 4 call site; E2E Task 8 |
| `+ New series` navigates to `/series/new` with seed | Yes | Task 4 (`navigateToNewSeries`), Task 3/8 |
| Scoping reads seed + scopes to those samples | Yes | Task 5 (`filterPickerBySeed`) |
| Cold-corpus "No proposable variable" → honest assign path | Yes | Task 6 (`ColdAssignPanel` + `buildColdAssignRows/canColdBuild`) |
| `Confirm & build` writes tags AND creates series | Yes | Task 7 (`scopeAndCreateSeriesMutator` sequential write) |
| Navigates to `/series/:id` (not `/series`) after create | Yes | Task 7 (`scopeAndCreate.data.id`) |
| No-legacy-port: `src/print/**` imports only new or carried | Yes | All new components import from `src/print/ui/` or `src/queries`/`src/api`/`src/state`/`src/lib/**` |
| Closed-look / open-placement: appearance in `src/print/ui/` | Yes | `Checkbox` and `ComposeBar` are both `src/print/ui/` primitives |
| `lint:design` clean | Yes | Gate on every task |
| TDD: failing test first for every meaningful change | Yes | Each task has test-first steps |
| No assertions on Tailwind class strings | Yes | All tests use `data-*` / roles |
| Existing corpus-culling E2E not regressed | Yes | `CullBar` untouched; separate `show` props; E2E Task 8 runs full suite |

**Type-consistency notes:**

- `SampleTableRow` currently exports `SAMPLE_TABLE_COLS` as a `const string`. Task 3 keeps the old export alias (`export const SAMPLE_TABLE_COLS = SAMPLE_TABLE_COLS_BASE`) so all existing importers (`SheetTable`, any tests) compile without change. The `sampleTableCols(bool)` function is the new call-site.
- `scopeAndCreateSeriesMutator.request` is async and returns `Promise<Series>`. The `Mutator<I, S, Series>` generic matches `saveSeries`'s existing use — verified against `saveSeriesMutator` at `/Users/me/projects/Himalaya.jl/.claude/worktrees/greenfield-ui-rebuild/packages/HimalayaUI/frontend/src/lib/queue/mutators/saveSeries.ts:57`.
- `useLocation()` is already in scope (imported by `useGlobalShortcuts` and `useStateFromUrl` — no new dep). `SeriesScopingPage.tsx` currently only imports `useNavigate`; Task 5 adds `useLocation`.
- `exactOptionalPropertyTypes` is enforced: every conditional spread in `ColdAssignPanel` and `scopeAndCreateSeriesMutator.request` uses `...(x !== undefined ? { key: x } : {})` — no `{ key: undefined }` assignable to optional fields.
agentId: a03f2399f44e50b20 (use SendMessage with to: 'a03f2399f44e50b20' to continue this agent)
<usage>subagent_tokens: 113905
tool_uses: 60
duration_ms: 460238</usage>
