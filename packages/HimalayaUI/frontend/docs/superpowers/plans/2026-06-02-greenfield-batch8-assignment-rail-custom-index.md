# Greenfield Batch 8 — Focus Assignment Rail + Custom-Index Modal

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development. Fresh implementer per task + frontend-reviewer pass. Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`. Commit ONLY named files (never `git add -A`; never stage `src/bones/registry.ts`, `src/bones/contact-sheet.bones.json`, or untracked plan docs other than this one). Typecheck: `npx tsc --noEmit -p tsconfig.build.json`. Branch `worktree-greenfield-ui-rebuild` stays UNMERGED — do not push/PR.

**Goal:** Build the full Focus assignment right-rail (the "whole assignment sidebar") and the custom-index modal it opens — the last two Layer-3 clusters needed to make the Focus page assemblable.

**Architecture:** Two coupled clusters from the focus mockups. (1) `AssignmentRail` — an `<aside>` wrapping two `RailSection`s over the *already-built* `AssignmentCart` (extended with its `.as-foot` "+ custom index…" trigger) and `CandidateList`. (2) `CustomIndexModal` — `ModalShell` + `ModalHead` + two `ModalFieldRow`s (a `SegmentedControl` symmetry picker + a `LatticeParamControl`) + the already-built `CustomPreview` viz + `FitMetadata` + `ModalFooter`. All panels are **presentational** — every interaction datum is a prop with a consumer-owned handler (Batch 7 contract). Two tiny refactor-on-contact primitive extensions (`Input.mono`, `SegmentedControl.stretch`).

**Tech Stack:** React 18 + TS (`exactOptionalPropertyTypes: true`) · Tailwind v4 `@theme` tokens · Storybook (`@storybook/react-vite`) · Vitest/RTL.

---

## Fidelity source (read these regions)

- `docs/redesign-mockups/focus-workspace.html` — the `.rail` / `.rail-sec` / `.rail-h` markup (assignment + candidates sections).
- `docs/redesign-mockups/2026-05-29-focus-plot.html` lines ~489–547 — the `.assign`/`.as-foot` cart footer **and** the complete `.scrim`/`.custom-sheet` modal. Also lines ~704–709 for the `SYMS` symmetry table (story data only).

### Verified token / role map (use these — never `text-[…px]`, never raw `oklch(`/hex)

| Mockup decl | Use |
|---|---|
| `.rail` `background:var(--paper-sunk); border-left:1px solid var(--hair)` | `bg-paper-sunk border-l border-hair` |
| `.rail-h` 10.5px/700/uppercase/0.09em/ink-faint | `Kicker tone="faint"` (text-kicker-faint = 11.5px/700/uppercase/0.10em/ink-faint — nearest named role; snap) |
| `.rail-h .cnt` ink-soft | `<span className="text-ink-soft">` |
| `.rail-note` 10.5px/ink-faint/lh1.55 | `text-caption text-ink-faint leading-relaxed` (text-caption is ink-soft → override `text-ink-faint`) |
| `.as-foot` dashed top, `+` accent prefix, 11.5px/600/ink-soft, hover bg | `border-t border-dashed border-hair-strong … text-caption font-semibold text-ink-soft hover:bg-paper-sunk hover:text-ink` + `<span className="text-print-accent font-bold">+</span>` |
| `.cs-head` border-bottom hair | `border-b border-hair` |
| `.cs-kicker` accent uppercase | `Kicker tone="accent"` |
| `.cs-t` Newsreader 20px/500 | `text-headline` (serif 19px/500 — snap 20→19, per ledger registry) |
| `.cs-x` × close | `IconButton dismiss label="Close"` |
| `.cs-flabel` 11.5px/500/uppercase/0.06em/ink-soft/width:100px | `text-label` (EXACT match) + width 100 via `style={{width:100}} className="shrink-0"` |
| `.cs-flabel .mono` (param name) | `<span className="font-mono normal-case text-ink-faint">` |
| `.cs-previewwrap` paper-sunk + hair + radius | `bg-paper-sunk border border-hair rounded-sm` |
| `.cs-fit` 12px/ink-soft, `b`→ink/700 | `text-meta text-ink-soft` (override soft) + `<b className="text-ink font-bold">`; snapped → `<span className="text-print-accent">` |
| `.cs-foot` flex space-between | `flex items-center justify-between` |
| `.cs-note` 10.5px/ink-faint | `text-caption text-ink-faint` |
| `.tool-btn` (Cancel) | `Button variant="outline"` |
| `.cs-add` accent fill (Add) | `Button variant="accent"` |
| `.scrim`/`.custom-sheet` | `ModalShell variant="dialog" size="md"` (provides scrim+blur+focus-trap+esc+outside-click; 640px ≈ mockup 600px) |

**Spacing** (px → nearest Tailwind step; arbitrary `gap-[..]`/`p-[..]` geometry is guard-clean if no named step fits): rail `p-5 pt-5 pb-7 gap-5`, rail-sec `gap-2.5`, cart footer `px-4 py-2.5`, modal head `px-5 py-4`, modal body `px-5 pt-4 pb-2 gap-4`, field row `gap-3.5`, lattice control `gap-2.5`, modal foot `px-5 pt-3.5 pb-4`, actions `gap-2`.

### Guard reminder (`src/print/components/**` is NOT exempt)
Composites are PLACEMENT-ONLY: primitives + layout utils + token utils (`bg-plate`/`bg-paper-sunk`/`border-hair`/`text-ink-faint`/`text-print-accent`) + named roles + plain utils (`uppercase`/`font-bold`/`font-mono`/`border-dashed`/`rounded-sm`). Inline `style` for GEOMETRY only (`width`, never color). `exactOptionalPropertyTypes`: forward an optional prop to a renderer whose matching prop is not `| undefined` via conditional spread `{...(x !== undefined ? { x } : {})}`.

### Per-component cadence (TDD)
`test/print-components/<Name>.test.tsx` (assert via `data-*`/roles/text, never class strings) → run → fail → build `src/print/components/<Name>.tsx` (or `ui/` for the two primitive extensions) → `src/print/**/<Name>.stories.tsx` (one default export per file; `@storybook/react-vite`; loose `Meta<typeof X>`; `@storybook/test` NOT installed — use a plain `const noop = () => {}`) → green on `npm test -- <Name>` + `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json` → commit named files.

---

## Task 1: `Input` — `mono` prop (refactor-on-contact, `src/print/ui/`)

**Files:** Modify `src/print/ui/Input.tsx`; Test `test/print-ui/Input.test.tsx` (append).

The lattice number field shows a monospace value. The inner `<input>` is hardcoded `text-base` (sans). Add an opt-in mono.

- [ ] **Step 1 — failing test:** append to `test/print-ui/Input.test.tsx`:
```tsx
it("renders the value in mono when mono is set", () => {
  const { getByTestId } = render(
    <Input value="252" onValueChange={() => {}} mono testId="latt" />,
  );
  const inner = getByTestId("latt").querySelector("input")!;
  expect(inner.className).toContain("font-mono");
});
it("is sans by default (no mono)", () => {
  const { getByTestId } = render(<Input value="x" onValueChange={() => {}} testId="d" />);
  expect(getByTestId("d").querySelector("input")!.className).not.toContain("font-mono");
});
```
- [ ] **Step 2:** run `npm test -- Input` → FAIL.
- [ ] **Step 3 — implement:** add `mono?: boolean;` to `InputProps`; destructure `mono = false`; on the inner `<input>` className append `mono && "font-mono"` (use a `cx`/template — keep existing classes, add `font-mono` when mono). The outer field appearance is unchanged.
- [ ] **Step 4:** run `npm test -- Input` → PASS. `npm run lint:design` clean. `npx tsc --noEmit -p tsconfig.build.json` clean.
- [ ] **Step 5 — commit:** `git add src/print/ui/Input.tsx test/print-ui/Input.test.tsx && git commit` (trailer).

---

## Task 2: `SegmentedControl` — `stretch` prop (refactor-on-contact, `src/print/ui/`)

**Files:** Modify `src/print/ui/SegmentedControl.tsx`; Test `test/print-ui/SegmentedControl.test.tsx` (append).

The symmetry picker fills its row with equal-width segments (`.cs-syms{flex:1}` + `.cs-syms button{flex:1}`).

- [ ] **Step 1 — failing test:** append:
```tsx
it("stretch makes the container full-width and segments flex-1", () => {
  const { getByTestId, getAllByRole } = render(
    <SegmentedControl stretch testId="syms" aria-label="sym"
      options={[{ value: "a", label: "A" }, { value: "b", label: "B" }]}
      value="a" onChange={() => {}} />,
  );
  expect(getByTestId("syms").className).toContain("w-full");
  expect(getAllByRole("button")[0]!.className).toContain("flex-1");
});
```
- [ ] **Step 2:** run → FAIL.
- [ ] **Step 3 — implement:** add `stretch?: boolean;` to props; destructure `stretch = false`. Container: when stretch, append `flex w-full` (keep `containerClass[variant]`). Segment button: append `stretch && "flex-1"`.
- [ ] **Step 4:** `npm test -- SegmentedControl` PASS; lint:design + tsc clean.
- [ ] **Step 5 — commit** the two files (trailer).

---

## Task 3: `AssignmentCart` — `onCustomIndex` footer (extend, `src/print/components/`)

**Files:** Modify `src/print/components/AssignmentCart.tsx`; Modify `test/print-components/AssignmentCart.test.tsx` (append); Modify `src/print/components/AssignmentCart.stories.tsx` (add a story showing the footer).

The `.as-foot` "+ custom index…" trigger is the LAST child inside the `.assign` plate, with a dashed top border. Shown only when `onCustomIndex` is provided; present in BOTH empty and filled states (you can always open custom index). It is what opens the modal.

- [ ] **Step 1 — failing test:** append:
```tsx
it("renders the custom-index footer when onCustomIndex is given, and fires it", () => {
  const onCustomIndex = vi.fn();
  const { getByTestId } = render(
    <AssignmentCart onCustomIndex={onCustomIndex}>
      <div>blockA</div>
    </AssignmentCart>,
  );
  const foot = getByTestId("custom-index-trigger");
  expect(foot.textContent).toContain("custom index");
  fireEvent.click(foot);
  expect(onCustomIndex).toHaveBeenCalledTimes(1);
});
it("shows the footer in the empty state too", () => {
  const { getByTestId } = render(<AssignmentCart onCustomIndex={() => {}} />);
  expect(getByTestId("assignment-empty")).toBeTruthy();
  expect(getByTestId("custom-index-trigger")).toBeTruthy();
});
it("omits the footer when onCustomIndex is absent", () => {
  const { queryByTestId } = render(<AssignmentCart><div>x</div></AssignmentCart>);
  expect(queryByTestId("custom-index-trigger")).toBeNull();
});
```
(ensure `fireEvent`, `vi` are imported in the test file.)
- [ ] **Step 2:** run → FAIL.
- [ ] **Step 3 — implement:** add `onCustomIndex?: () => void;` to `AssignmentCartProps`. After the empty/filled branch (still inside the root `.assign` div), render when `onCustomIndex`:
```tsx
{onCustomIndex && (
  <button
    type="button"
    data-testid="custom-index-trigger"
    onClick={onCustomIndex}
    className="w-full text-left border-t border-dashed border-hair-strong px-4 py-2.5 text-caption font-semibold text-ink-soft transition-colors hover:bg-paper-sunk hover:text-ink"
  >
    <span className="text-print-accent font-bold">+</span> custom index…
  </button>
)}
```
- [ ] **Step 4 — story:** add an `WithCustomIndex` story (PhaseBlocks + `onCustomIndex={noop}`).
- [ ] **Step 5:** `npm test -- AssignmentCart` PASS; lint:design + tsc clean.
- [ ] **Step 6 — commit** the three files (trailer).

---

## Task 4: `RailSection` (new, `src/print/components/`)

**Files:** Create `src/print/components/RailSection.tsx`, `src/print/components/RailSection.stories.tsx`, `test/print-components/RailSection.test.tsx`.

Reusable `.rail-sec`: a kicker header row (faint label + optional count) + children + optional note. (BuilderRail will reuse it later.)

- [ ] **Step 1 — failing test:** `test/print-components/RailSection.test.tsx`:
```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { RailSection } from "../../src/print/components/RailSection";

describe("RailSection", () => {
  it("renders the label, optional count, children, and note", () => {
    const { getByTestId, getByText } = render(
      <RailSection label="Assignment" count={3} note="ranked note">
        <div>body</div>
      </RailSection>,
    );
    const head = getByTestId("rail-section-head");
    expect(head.textContent).toContain("Assignment");
    expect(head.textContent).toContain("3");
    expect(getByText("body")).toBeTruthy();
    expect(getByTestId("rail-section-note").textContent).toContain("ranked note");
  });
  it("omits count and note when not given", () => {
    const { queryByTestId, getByTestId } = render(
      <RailSection label="Candidates"><div>x</div></RailSection>,
    );
    expect(getByTestId("rail-section-head").textContent).toContain("Candidates");
    expect(queryByTestId("rail-section-note")).toBeNull();
  });
});
```
- [ ] **Step 2:** run → FAIL.
- [ ] **Step 3 — implement:**
```tsx
import type { ReactNode } from "react";
import { Kicker } from "../ui";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface RailSectionProps {
  label: ReactNode;
  count?: ReactNode;
  note?: ReactNode;
  children?: ReactNode;
  className?: string;
}

export function RailSection({ label, count, note, children, className }: RailSectionProps): JSX.Element {
  return (
    <div data-testid="rail-section" className={cx("flex flex-col gap-2.5", className)}>
      <div data-testid="rail-section-head" className="flex items-baseline justify-between">
        <Kicker tone="faint">
          {label}
          {count != null && <span className="text-ink-soft"> {count}</span>}
        </Kicker>
      </div>
      {children}
      {note != null && (
        <div data-testid="rail-section-note" className="text-caption text-ink-faint leading-relaxed">
          {note}
        </div>
      )}
    </div>
  );
}
```
- [ ] **Step 4 — story:** a `Default` showing label+count+note around a placeholder.
- [ ] **Step 5:** green (`npm test -- RailSection`, lint:design, tsc).
- [ ] **Step 6 — commit** the three files (trailer).

---

## Task 5: `AssignmentRail` (new, `src/print/components/`)

**Files:** Create `src/print/components/AssignmentRail.tsx`, `.stories.tsx`, `test/print-components/AssignmentRail.test.tsx`.

The full right sidebar. Presentational: the `assignment` (an `AssignmentCart`) and `candidates` (a `CandidateList`) are passed as slots; counts/notes as props. The story wires real `PhaseBlock`/`CandidateRow` children from fixtures.

- [ ] **Step 1 — failing test:**
```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { AssignmentRail } from "../../src/print/components/AssignmentRail";

describe("AssignmentRail", () => {
  it("renders an aside with Assignment (count) and Candidates sections wrapping the slots", () => {
    const { getByTestId, getByText, getAllByTestId } = render(
      <AssignmentRail
        assignmentCount={2}
        assignment={<div data-testid="slot-cart">CART</div>}
        candidates={<div data-testid="slot-cands">CANDS</div>}
        candidatesNote="ranked by fit"
      />,
    );
    expect(getByTestId("assignment-rail").tagName).toBe("ASIDE");
    const heads = getAllByTestId("rail-section-head").map((h) => h.textContent);
    expect(heads.some((t) => t?.includes("Assignment") && t?.includes("2"))).toBe(true);
    expect(heads.some((t) => t?.includes("Candidates"))).toBe(true);
    expect(getByTestId("slot-cart")).toBeTruthy();
    expect(getByTestId("slot-cands")).toBeTruthy();
    expect(getByText("ranked by fit")).toBeTruthy();
  });
});
```
- [ ] **Step 2:** run → FAIL.
- [ ] **Step 3 — implement:**
```tsx
import type { ReactNode } from "react";
import { RailSection } from "./RailSection";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface AssignmentRailProps {
  /** The AssignmentCart element (with its PhaseBlock children + onCustomIndex). */
  assignment: ReactNode;
  /** Count shown beside the "Assignment" kicker. */
  assignmentCount?: ReactNode;
  /** Optional status note under the cart (the mockup .rail-note#assign-note). */
  assignmentNote?: ReactNode;
  /** The CandidateList element (with its CandidateRow children). */
  candidates: ReactNode;
  /** The candidates rationale note (Bonnet explainer). */
  candidatesNote?: ReactNode;
  className?: string;
}

export function AssignmentRail({
  assignment, assignmentCount, assignmentNote,
  candidates, candidatesNote, className,
}: AssignmentRailProps): JSX.Element {
  return (
    <aside
      data-testid="assignment-rail"
      className={cx(
        "bg-paper-sunk border-l border-hair overflow-y-auto p-5 pb-7 flex flex-col gap-5",
        className,
      )}
    >
      <RailSection label="Assignment" {...(assignmentCount != null ? { count: assignmentCount } : {})} {...(assignmentNote != null ? { note: assignmentNote } : {})}>
        {assignment}
      </RailSection>
      <RailSection label="Candidates" {...(candidatesNote != null ? { note: candidatesNote } : {})}>
        {candidates}
      </RailSection>
    </aside>
  );
}
```
(Conditional spreads satisfy `exactOptionalPropertyTypes` against `RailSection`'s optional `count`/`note`.)
- [ ] **Step 4 — story:** wire a real composition — `assignment={<AssignmentCart onCustomIndex={noop}><PhaseBlock .../><PhaseBlock .../></AssignmentCart>}`, `candidates={<CandidateList><CandidateRow .../>…</CandidateList>}`, using `realMembers`/phase data already used by existing PhaseBlock/CandidateRow stories. Give the aside a fixed width wrapper (`style={{ width: 340 }}` or a `w-[340px]`) so it reads like the page rail. Note copy from the mockup: `Ranked by fit; Bonnet bumps a coexisting cubic whose lattice matches the Gauss–Bonnet ratio. Add the ones that make sense.`
- [ ] **Step 5:** green; visually verify in Storybook against `focus-workspace.html`.
- [ ] **Step 6 — commit** the three files (trailer).

---

## Task 6: `ModalHead` (new, `src/print/components/`)

**Files:** Create `src/print/components/ModalHead.tsx`, `.stories.tsx`, `test/print-components/ModalHead.test.tsx`. (Reused by the future ExportSheet.)

- [ ] **Step 1 — failing test:** assert kicker text, an `<h2>` (or the `titleId`-bearing element) with the title text, and a close button (`aria-label="Close"`) that fires `onClose`.
- [ ] **Step 2:** FAIL.
- [ ] **Step 3 — implement:**
```tsx
import type { ReactNode } from "react";
import { Kicker, IconButton } from "../ui";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface ModalHeadProps {
  kicker: ReactNode;
  title: ReactNode;
  /** id wired to the dialog's aria-labelledby. */
  titleId?: string;
  onClose: () => void;
  className?: string;
}

export function ModalHead({ kicker, title, titleId, onClose, className }: ModalHeadProps): JSX.Element {
  return (
    <div data-testid="modal-head" className={cx("flex items-start justify-between px-5 py-4 border-b border-hair", className)}>
      <div>
        <Kicker tone="accent">{kicker}</Kicker>
        <h2 {...(titleId ? { id: titleId } : {})} className="text-headline mt-0.5">{title}</h2>
      </div>
      <IconButton dismiss label="Close" tone="ghost" onClick={onClose} className="shrink-0 -mr-1" />
    </div>
  );
}
```
- [ ] **Step 4 — story** (Default). **Step 5** green. **Step 6** commit (trailer).

---

## Task 7: `ModalFieldRow` (new, `src/print/components/`)

**Files:** Create `src/print/components/ModalFieldRow.tsx`, `.stories.tsx`, `test/print-components/ModalFieldRow.test.tsx`.

`.cs-row`: a fixed-width uppercase label + a flex control slot. Used for both the Symmetry and Lattice rows.

- [ ] **Step 1 — failing test:** assert the label text renders with `text-label` role context (assert text + that a `data-testid="modal-field-row"` exists and the control slot child renders), label width is fixed (assert `style.width`).
- [ ] **Step 2:** FAIL.
- [ ] **Step 3 — implement:**
```tsx
import type { ReactNode } from "react";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface ModalFieldRowProps {
  label: ReactNode;
  /** Optional element appended after the label text (e.g. the mono param name). */
  labelSuffix?: ReactNode;
  children: ReactNode;        // the control (segmented control / lattice control)
  className?: string;
}

export function ModalFieldRow({ label, labelSuffix, children, className }: ModalFieldRowProps): JSX.Element {
  return (
    <div data-testid="modal-field-row" className={cx("flex items-center gap-3.5", className)}>
      <span data-testid="modal-field-label" className="text-label shrink-0" style={{ width: 100 }}>
        {label}
        {labelSuffix != null && <span className="font-mono normal-case text-ink-faint"> {labelSuffix}</span>}
      </span>
      <div className="flex-1 min-w-0 flex items-center">{children}</div>
    </div>
  );
}
```
- [ ] **Step 4 — story.** **Step 5** green. **Step 6** commit (trailer).

---

## Task 8: `LatticeParamControl` (new, `src/print/components/`)

**Files:** Create `src/print/components/LatticeParamControl.tsx`, `.stories.tsx`, `test/print-components/LatticeParamControl.test.tsx`. Depends on `Slider` + `Input` (with the Task-1 `mono`).

`.cs-param`: bare slider (flex-1) + mono number field (64px) + unit. Presentational: `value` (string) + `onValueChange`, `min`/`max`/`step`, `unit`.

- [ ] **Step 1 — failing test:**
```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { LatticeParamControl } from "../../src/print/components/LatticeParamControl";

describe("LatticeParamControl", () => {
  it("renders a slider + a mono number field + unit, both reflecting value", () => {
    const { getByTestId, getByText } = render(
      <LatticeParamControl value="252" min={120} max={360} step={1} onValueChange={() => {}} unit="Å" />,
    );
    const slider = getByTestId("slider") as HTMLInputElement;
    expect(slider.value).toBe("252");
    const num = getByTestId("lattice-num").querySelector("input") as HTMLInputElement;
    expect(num.value).toBe("252");
    expect(num.className).toContain("font-mono");
    expect(getByText("Å")).toBeTruthy();
  });
  it("fires onValueChange from both controls", () => {
    const onValueChange = vi.fn();
    const { getByTestId } = render(
      <LatticeParamControl value="252" min={120} max={360} step={1} onValueChange={onValueChange} unit="Å" />,
    );
    fireEvent.change(getByTestId("slider"), { target: { value: "260" } });
    fireEvent.change(getByTestId("lattice-num").querySelector("input")!, { target: { value: "300" } });
    expect(onValueChange).toHaveBeenNthCalledWith(1, "260");
    expect(onValueChange).toHaveBeenNthCalledWith(2, "300");
  });
});
```
- [ ] **Step 2:** FAIL.
- [ ] **Step 3 — implement:** (Slider emits `number` → stringify; Input emits string.)
```tsx
import { Slider, Input } from "../ui";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface LatticeParamControlProps {
  value: string;            // current lattice value as a string (mono display)
  min: number;
  max: number;
  step?: number;
  onValueChange: (v: string) => void;
  unit?: string;            // e.g. "Å"
  /** a11y name for the slider/field (defaults below). */
  label?: string;
  className?: string;
}

export function LatticeParamControl({
  value, min, max, step, onValueChange, unit = "Å", label = "lattice parameter", className,
}: LatticeParamControlProps): JSX.Element {
  return (
    <div data-testid="lattice-param" className={cx("flex items-center gap-2.5 flex-1 min-w-0", className)}>
      <Slider
        value={Number(value)}
        min={min}
        max={max}
        {...(step !== undefined ? { step } : {})}
        onChange={(n) => onValueChange(String(n))}
        label={label}
        className="flex-1 min-w-0"
      />
      <Input
        testId="lattice-num"
        type="number"
        mono
        inputSize="sm"
        value={value}
        onValueChange={onValueChange}
        min={min}
        max={max}
        {...(step !== undefined ? { step } : {})}
        aria-label={`${label} value`}
        className="w-16 shrink-0"
      />
      <span className="text-caption text-ink-faint font-mono shrink-0">{unit}</span>
    </div>
  );
}
```
Note: `Slider` renders its `label` into a visible label row; here we want the BARE range (the mockup has no inline label — the row label lives in `ModalFieldRow`). Pass the label for a11y but suppress the visible row: if `Slider` always renders the label row when `label` is set, instead pass NO `label` and add `aria-label` via a wrapping technique — **verify `Slider`'s render**: it only shows the row `if (label || valueDisplay)`. To keep the slider bare AND accessible, pass `label` (acceptable: a tiny faint label row) OR omit `label` (no a11y name). **Decision:** omit the visible row — pass no `label`, and rely on `ModalFieldRow`'s label for sighted users; for AT, wrap the slider's `aria-label` is only set when `label` is passed. To get an a11y name without the visible row, the cleanest is a `Slider` that always sets `aria-label` — but it doesn't. So: pass `label` (the faint row is acceptable and matches no mockup element but is harmless) **OR** extend `Slider` with an `aria-label`-only prop. To avoid a third primitive edit, pass `label` and accept the small label row. If the reviewer flags the extra row as a fidelity miss, the follow-up is a `Slider` `hideLabel`/`aria-label` prop — note it, don't block.
- [ ] **Step 4 — story** (interactive `useState` value). **Step 5** green. **Step 6** commit (trailer).

---

## Task 9: `FitMetadata` (new, `src/print/components/`)

**Files:** Create `src/print/components/FitMetadata.tsx`, `.stories.tsx`, `test/print-components/FitMetadata.test.tsx`.

`.cs-fit`: "Lands on **N** of M observed peaks · **{param}** = **{value} {unit}**[ · snapped]".

- [ ] **Step 1 — failing test:** assert the rendered text contains `Lands on`, the `landed`/`total` numbers, the `paramName`, `paramValue`, unit, and that `snapped` text appears only when `snapped` is true (assert a `data-testid="fit-snapped"` present/absent).
- [ ] **Step 2:** FAIL.
- [ ] **Step 3 — implement:**
```tsx
function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface FitMetadataProps {
  landed: number;
  total: number;
  paramName: string;     // "a" | "d"
  paramValue: string;    // e.g. "252"
  unit?: string;         // "Å"
  snapped?: boolean;
  className?: string;
}

export function FitMetadata({ landed, total, paramName, paramValue, unit = "Å", snapped = false, className }: FitMetadataProps): JSX.Element {
  return (
    <div data-testid="fit-metadata" className={cx("text-meta text-ink-soft", className)}>
      Lands on <b className="text-ink font-bold">{landed}</b> of {total} observed peaks
      {" · "}
      <b className="text-ink font-bold">{paramName}</b> = <b className="text-ink font-bold">{paramValue} {unit}</b>
      {snapped && (
        <> · <span data-testid="fit-snapped" className="text-print-accent">snapped</span></>
      )}
    </div>
  );
}
```
- [ ] **Step 4 — story** (snapped + not-snapped). **Step 5** green. **Step 6** commit (trailer).

---

## Task 10: `ModalFooter` (new, `src/print/components/`)

**Files:** Create `src/print/components/ModalFooter.tsx`, `.stories.tsx`, `test/print-components/ModalFooter.test.tsx`. (Reused by ExportSheet.)

`.cs-foot`: a faint note slot (left) + an actions slot (right).

- [ ] **Step 1 — failing test:** assert `note` and `actions` children both render inside `data-testid="modal-foot"`.
- [ ] **Step 2:** FAIL.
- [ ] **Step 3 — implement:**
```tsx
import type { ReactNode } from "react";

function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface ModalFooterProps {
  note?: ReactNode;
  actions: ReactNode;
  className?: string;
}

export function ModalFooter({ note, actions, className }: ModalFooterProps): JSX.Element {
  return (
    <div data-testid="modal-foot" className={cx("flex items-center justify-between px-5 pt-3.5 pb-4", className)}>
      {note != null ? <span className="text-caption text-ink-faint">{note}</span> : <span />}
      <div className="flex gap-2">{actions}</div>
    </div>
  );
}
```
- [ ] **Step 4 — story.** **Step 5** green. **Step 6** commit (trailer).

---

## Task 11: `CustomIndexModal` (new, `src/print/components/`)

**Files:** Create `src/print/components/CustomIndexModal.tsx`, `.stories.tsx`, `test/print-components/CustomIndexModal.test.tsx`. Depends on `ModalShell`, `ModalHead`, `ModalFieldRow`, `SegmentedControl` (+ Task-2 stretch), `LatticeParamControl`, `CustomPreview` (`../comb`), `FitMetadata`, `ModalFooter`, `Button`.

Fully presentational. `SymmetrySelector` is folded in as a `SegmentedControl` (no separate file — it IS a segmented control; the `SYMS` domain table is the consumer's/story's concern). The `CustomPreview` SVG is wrapped in the `.cs-previewwrap` plate.

- [ ] **Step 1 — failing test:**
```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { CustomIndexModal } from "../../src/print/components/CustomIndexModal";
import { PN3M } from "../../src/print/comb/comb.fixtures";

const base = {
  open: true,
  symmetries: ["Pn3m", "Im3m", "Ia3d", "Lamellar", "Hexagonal"],
  symmetry: "Im3m",
  paramName: "a",
  paramValue: "252",
  paramMin: 120, paramMax: 360, paramStep: 1,
  previewSeries: PN3M,
  observed: [0.03, 0.045, 0.06],
  fit: { landed: 3, total: 5, snapped: true },
};

describe("CustomIndexModal", () => {
  it("renders head, symmetry seg, lattice control, preview, fit line and footer when open", () => {
    const { getByTestId, getByText, getAllByRole } = render(
      <CustomIndexModal {...base} onClose={() => {}} onCancel={() => {}} onAdd={() => {}}
        onSymmetryChange={() => {}} onParamChange={() => {}} />,
    );
    expect(getByText("Custom index")).toBeTruthy();
    expect(getByText("Speculative")).toBeTruthy();
    // symmetry segmented control: the active one is Im3m
    expect(getAllByRole("button").some((b) => b.textContent === "Im3m")).toBe(true);
    expect(getByTestId("lattice-param")).toBeTruthy();
    expect(getByTestId("custom-preview")).toBeTruthy();
    expect(getByTestId("fit-metadata").textContent).toContain("Lands on");
    expect(getByText("Add to assignment")).toBeTruthy();
    expect(getByText("Cancel")).toBeTruthy();
  });
  it("renders nothing when open=false", () => {
    const { queryByText } = render(
      <CustomIndexModal {...base} open={false} onClose={() => {}} onCancel={() => {}} onAdd={() => {}}
        onSymmetryChange={() => {}} onParamChange={() => {}} />,
    );
    expect(queryByText("Custom index")).toBeNull();
  });
  it("wires the actions", () => {
    const onAdd = vi.fn(), onCancel = vi.fn(), onSymmetryChange = vi.fn();
    const { getByText } = render(
      <CustomIndexModal {...base} onClose={() => {}} onCancel={onCancel} onAdd={onAdd}
        onSymmetryChange={onSymmetryChange} onParamChange={() => {}} />,
    );
    fireEvent.click(getByText("Add to assignment"));
    fireEvent.click(getByText("Cancel"));
    fireEvent.click(getByText("Pn3m"));
    expect(onAdd).toHaveBeenCalledTimes(1);
    expect(onCancel).toHaveBeenCalledTimes(1);
    expect(onSymmetryChange).toHaveBeenCalledWith("Pn3m");
  });
});
```
- [ ] **Step 2:** FAIL.
- [ ] **Step 3 — implement:**
```tsx
import { ModalShell, Button } from "../ui";
import { SegmentedControl } from "../ui";
import { CustomPreview } from "../comb";
import type { CombSeries } from "../comb/combModel";
import { ModalHead } from "./ModalHead";
import { ModalFieldRow } from "./ModalFieldRow";
import { LatticeParamControl } from "./LatticeParamControl";
import { FitMetadata } from "./FitMetadata";
import { ModalFooter } from "./ModalFooter";

export interface CustomIndexFit {
  landed: number;
  total: number;
  snapped?: boolean;
}

export interface CustomIndexModalProps {
  open: boolean;
  onClose: () => void;
  onCancel: () => void;
  onAdd: () => void;
  symmetries: readonly string[];
  symmetry: string;
  onSymmetryChange: (s: string) => void;
  paramName: string;          // "a" | "d"
  paramValue: string;
  paramMin: number;
  paramMax: number;
  paramStep?: number;
  onParamChange: (v: string) => void;
  unit?: string;
  previewSeries: CombSeries;
  observed: number[];
  fit: CustomIndexFit;
  className?: string;
}

export function CustomIndexModal(props: CustomIndexModalProps): JSX.Element | null {
  const {
    open, onClose, onCancel, onAdd,
    symmetries, symmetry, onSymmetryChange,
    paramName, paramValue, paramMin, paramMax, paramStep, onParamChange, unit = "Å",
    previewSeries, observed, fit, className,
  } = props;

  return (
    <ModalShell
      open={open}
      onClose={onClose}
      size="md"
      testId="custom-index-modal"
      aria-labelledby="cix-title"
      {...(className ? { className } : {})}
    >
      <ModalHead kicker="Speculative" title="Custom index" titleId="cix-title" onClose={onClose} />

      <div className="px-5 pt-4 pb-2 flex flex-col gap-4">
        <ModalFieldRow label="Symmetry">
          <SegmentedControl
            stretch
            aria-label="Crystal symmetry"
            options={symmetries.map((s) => ({ value: s, label: s }))}
            value={symmetry}
            onChange={onSymmetryChange}
          />
        </ModalFieldRow>

        <ModalFieldRow label="Lattice" labelSuffix={paramName}>
          <LatticeParamControl
            value={paramValue}
            min={paramMin}
            max={paramMax}
            {...(paramStep !== undefined ? { step: paramStep } : {})}
            onValueChange={onParamChange}
            unit={unit}
          />
        </ModalFieldRow>

        <div className="bg-paper-sunk border border-hair rounded-sm px-2 py-1.5">
          <CustomPreview series={previewSeries} observed={observed} className="block w-full" />
        </div>

        <FitMetadata
          landed={fit.landed}
          total={fit.total}
          paramName={paramName}
          paramValue={paramValue}
          unit={unit}
          {...(fit.snapped !== undefined ? { snapped: fit.snapped } : {})}
        />
      </div>

      <ModalFooter
        note="Drag the lattice until the teeth land on your peaks."
        actions={
          <>
            <Button variant="outline" onClick={onCancel}>Cancel</Button>
            <Button variant="accent" onClick={onAdd}>Add to assignment</Button>
          </>
        }
      />
    </ModalShell>
  );
}
```
**Note:** confirm the `ui` barrel exports `ModalShell`, `SegmentedControl`, `Button`; confirm `../comb` barrel exports `CustomPreview` and `combModel` exports `CombSeries` — adjust imports to the real barrels (verify before writing).
- [ ] **Step 4 — story:** an interactive `Open` story with a `useState` harness: a small in-story `SYMS` table (param defaults/min/max per symmetry, copied from `2026-05-29-focus-plot.html` ~704–709) driving `symmetry`/`paramValue`, a fixed `previewSeries` from `comb.fixtures` (`PN3M`/`IM3M`) + `observed` peaks, and a `fit` object. A trigger `<button>` toggles `open`. Also a static `Closed` story is unnecessary (open=false renders null) — instead add a `Pn3m`/`Im3m` variant story showing two symmetries.
- [ ] **Step 5:** green (`npm test -- CustomIndexModal`, lint:design, tsc). Visually verify in Storybook against the `.custom-sheet` mockup.
- [ ] **Step 6 — commit** the three files (trailer).

---

## Final verification + ledger

- [ ] Full gate: `npm test` (full Vitest) green · `npm run lint:design` clean · `npx tsc --noEmit -p tsconfig.build.json` exit 0 · `npm run build-storybook` exit 0.
- [ ] Visual: build-storybook → serve → Playwright screenshot `AssignmentRail` (vs `focus-workspace.html`) and `CustomIndexModal` (vs `.custom-sheet`).
- [ ] Update `docs/greenfield-component-ledger.md`: flip `AssignmentRail`(new row), the `CustomIndexModal` cluster rows (`CustomIndexModal`/`ModalHead`/`SymmetrySelector`→folded/`LatticeParamControl`/`FitMetadata`/`ModalFooter`), and the `.as-foot` footer note → ✅; add `RailSection`/`ModalFieldRow` rows; record the Batch-8 decisions in the registry (SymmetrySelector folded into SegmentedControl; Input.mono + SegmentedControl.stretch refactors-on-contact; AssignmentRail/CustomIndexModal presentational; Slider visible-label-row follow-up if flagged). Bump the coverage summary. Commit the ledger.
- [ ] Update memory `project_greenfield_composite_layer.md` + `MEMORY.md` with a Batch-8 paragraph.

## Decisions locked (don't re-litigate)
- **AssignmentRail + CustomIndexModal are presentational** — all symmetry/param/preview/fit state + handlers are props (Batch 7 contract). The `SYMS` table, snap math, and predicted-teeth computation live in the future Focus page, not here.
- **SymmetrySelector is folded into a `SegmentedControl`** (with the new `stretch`) — no separate file; it carries no logic the segmented control doesn't.
- **The `.as-foot` "custom index…" trigger lives inside `AssignmentCart`** (deferred from Batch 6) — it's the mockup's last child of `.assign`, shown whenever `onCustomIndex` is provided.
- **`ModalHead`/`ModalFooter`/`ModalFieldRow` are generic** — built here, reused by the later `ExportSheet` slice.
