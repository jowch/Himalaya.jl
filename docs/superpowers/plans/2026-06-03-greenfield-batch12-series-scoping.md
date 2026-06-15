# Greenfield Batch 12 — Series-scoping slice (`FlagButton` + `ScopeSampleRow` + `ScopeCandidateRow` + `ScopePlate`)

> Execute with superpowers:subagent-driven-development (fresh implementer per task + two-stage review).
> Commit trailer (every commit, exact last line): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Worktree frontend dir: `packages/HimalayaUI/frontend`. Typecheck: `npx tsc --noEmit -p tsconfig.build.json`.
> Commit ONLY named files. NEVER `git add -A`/`.`. NEVER stage `src/bones/*` or anything under `docs/superpowers/plans/`.

## Context

The greenfield "The Print" rebuild (branch `worktree-greenfield-ui-rebuild`, NOT merged). Layers 0/1/2 complete; the build frontier is the last two Layer-3 panels. This batch builds the **series-scoping** surface, derived top-to-bottom from `docs/redesign-mockups/series-scoping.html` (read it — it is the source of truth).

The scoping surface is a **review worksheet**, not a form: the user arrives with samples already multi-selected on the contact sheet, Himalaya has grouped them and parsed an ordering value from each sample's name, and the human glances, fixes the one flagged parse, and confirms. Confirming is what records the metadata.

**Decomposition (4 tasks, dependency-ordered):**

1. **`FlagButton`** (new Layer-0 primitive, `ui/`) — the `.s-val` parsed-value control. Two states (`flagged` / `ok`). This is the slice's single refactor-on-contact.
2. **`ScopeSampleRow`** (`components/`) — the `.srow`: grip + trace sparkline + name/id + a `FlagButton`. Needs Task 1.
3. **`ScopeCandidateRow`** (`components/`) — the `.cand-row`: a dashed *sample* candidate (sparkline + name + "why" + Add button). Independent of Tasks 1–2.
4. **`ScopePlate`** (`components/`) + a `ScopingAssembly` page-sim story — the whole `.scope-plate` worksheet, presentational, composing Tasks 2–3 + existing primitives. Needs Tasks 1–3.

**Why `ScopeCandidateRow` is NOT the existing `CandidateRow`:** the Focus `CandidateRow` (`components/CandidateRow.tsx`) is a *phase* candidate (PhaseChip + ScoreBar + numeric score). The scoping `.cand-row` is a *sample* candidate (sparkline + sample name + a one-line "why" + "Add to series"). Structurally distinct — same precedent as the `SeriesMemberRow`-split-from-`MemberRow` decision. Do NOT reuse or modify `CandidateRow`.

**Layer-4 deferral (consistent with Batch 7/9/10/11):** the topbar (`Wordmark`/`StageTabs`/Discard), the centred `.stage` scroller, the dropdown MENU behind the Ordered-by field, the editable-title text input, and the durable confirm/undo wiring are **page assembly**. `ScopePlate` is presentational (every interaction datum is a prop with a consumer-owned handler — the Batch 7 contract). The `ScopingAssembly` story simulates the page (owns the `SERIES`/`CANDIDATES` arrays + flag-toggle + undo + derived preview/foot) for fidelity without building the Layer-4 page.

### Verified reuse APIs (checked against live source 2026-06-03 — do NOT re-derive)

- **`Sparkline`** — import `from "../plot/Sparkline"` (it is **NOT** re-exported from `../plot/index.ts`, verified 2026-06-03). Props `{ trace: { q: number[]; I: number[] }; phase?: string | null; className? }`. Renders a 76×28 `bg-paper-sunk` framed mini area-trace; `data-testid="sparkline"`. **`className` is PLACEMENT ONLY** — to dim the candidate spark to the mockup's `opacity:0.7`, pass `className="opacity-70"` (opacity is allowed in `components/`). Feed it a real fixture trace, NOT the mockup's synthetic intensity model.
- **`GripHandle`** (`../ui`) — `{ className? }`. Renders `⋮⋮`, `aria-hidden`, `text-hair-strong group-hover:text-ink-faint cursor-grab`. The consumer row MUST carry `group` for the hover reveal. `data-testid="grip-handle"`. In scoping the grip is the mockup's visual affordance; drag-reorder wiring is page-deferred (do not invent it).
- **`PhaseStrip`** (`../ui`) — `{ segments: PhaseSegment[]; size?: "sm"|"md"; emptyLabel?; orientation?; className? }`. `PhaseSegment = { phase: string|null; coexistWith?: string[]|null; state?: "form_factor"|"null" }`. **Use `size="sm"`** — its doc explicitly calls `sm` "the legacy 8px Scoping bar." Renders cells + a derived "A → B / throughout" caption beneath; `data-testid` `ps-seg`/`ps-cap`. The scoping preview's coexistence cells (a sample with 2 phases) use `coexistWith: ["Lamellar"]` on a `phase:"Pn3m"` segment.
- **`AutoGroup`** (`./AutoGroup`) — `{ variant?: "summary"|"compose"; title?; children: ReactNode; actions?; className? }`. **Use `variant="summary"`** (recessed `bg-paper-sunk`, the scoping look) with the grouping sentence as `children` (embed `<strong>` for the bold emphasis the mockup shows on "6 samples" / "LL37 : lipid ratio"). It already renders the accent star. `data-testid="auto-group"`.
- **`Kicker`** (`../ui`) — `{ tone?: "accent"|"faint"; as? }`. `tone="accent"` = the "New series" kicker; `tone="faint"` = the uppercase/700/tracked/ink-faint **section/field labels** ("Ordered by", "The series", "Himalaya also found", "Preview — phase across the series"). The mockup's field labels are 10.5px/0.09em; snap to the named `Kicker` role (the "snap to named role" decision — do NOT author `text-[10.5px]`).
- **`Dot`** (`../ui`) — `{ tone?: "accent"|"success"|"muted"|"neutral"; size?; label?; bordered? }`. Foot state: warn → `tone="accent"`, ready → `tone="success"` (the system success token replaces the mockup's raw `--im3m`). Give it a `label` (it gets `role="img"`) only if it is the sole signal; here the adjacent text carries the meaning, so `aria-hidden` is fine.
- **`Button`** (`../ui`) — `variant?: "solid"|"accent"|"ghost"|"danger"|"outline"|"ghostInverse"`. The dark **Confirm & build** button = `variant="solid"` (`bg-ink text-paper`); pass `disabled` for the gated state (Button forwards `disabled` + `...props`). The candidate **"+ Add to series"** button = `variant="outline"` (the mockup's `bg-plate` + hairline). `data-variant` reflects the variant.
- **`Card`** (`../ui`) — `elevated` → bg-plate + hairline + radius + soft warm shadow = the `.scope-plate` chrome.
- **`cx` helper** — NO shared export; copy the house idiom into each new file: `function cx(...parts: Array<string|false|null|undefined>): string { return parts.filter(Boolean).join(" "); }`.
- **Fixtures** — `src/print/fixtures/realTraces.ts` exports `realTraces: Record<number, Trace>` keyed `37/65/66/67/93` (`Trace = { q, I, sigma }`). Stories feed these to `Sparkline`. (The scoping mockup's sample ids `smp_04…` are page data; in the story, map each to a real fixture trace + a phase.)

### Design-guard contract (`npm run lint:design` MUST stay green)

`src/print/ui/**` **IS** guard-exempt (FlagButton may author appearance there). `src/print/components/**` is **NOT**. In `components/`:
- **ALLOWED:** token color classes (`bg-plate`, `bg-paper-sunk`, `bg-accent`, `bg-accent/5`, `text-ink`/`text-ink-soft`/`text-ink-faint`, `border-hair`/`border-hair-strong`, `text-paper`), named type roles (`text-xs`/`text-sm`/`text-base`/`text-body`/`text-meta`/`text-caption`/`text-data`/`text-display`/`text-headline`), `font-mono`/`font-semibold`/`font-bold`/`uppercase`/`tracking-*`, `opacity-70`, arbitrary SPACING/geometry (`min-w-[92px]`, `max-w-[760px]`, `gap-3`, `px-2`), `border-dashed`, `rounded`/`rounded-sm`/`rounded-md`/`rounded-full`.
- **BANNED:** arbitrary-value APPEARANCE (`text-[15px]`, `rounded-[7px]`, `bg-[…]`, `border-[…color…]`), raw `oklch(`/hex/`color-mix(` literals, side-stripe borders. If appearance is needed that a primitive doesn't expose → extend the primitive in `ui/` (only Task 1 does this).
- **Tests assert via `data-*`/roles/text — NEVER class strings.**

---

## Task 1: `FlagButton` primitive (`ui/`)

**Mockup region:** `series-scoping.html` `.s-val`, `.s-val.flag`, `.s-val.ok`, `.s-val .v`, `.s-val .check` (lines 195–222), and the `rowHTML` value branch (lines 452–460).

**The control:** a right-aligned, mono, clickable parsed value. Two states:
- **`flagged`** — `color: var(--accent)`, `cursor: pointer`; the value `.v` carries a **dashed** accent-tinted bottom border; below it a `.check` caption "▸ check the read" (9px/700/uppercase/accent, **sans** font not mono). Clicking it accepts the read (the consumer flips state).
- **`ok`** (default) — `color: var(--ink)`, `cursor: pointer`, `title="Click to re-open this value"`; the value `.v` has a **dotted transparent** bottom border that becomes `--hair-strong` on hover (the re-openable affordance). No caption.

It is a `<button type="button">`. Right-aligned text, `min-width: 92px`, value font 13px/700 mono.

**Files:**
- Create: `src/print/ui/FlagButton.tsx`
- Create: `src/print/ui/FlagButton.stories.tsx`
- Modify: `src/print/ui/index.ts` (add `export { FlagButton } from "./FlagButton"; export type { FlagButtonProps } from "./FlagButton";`)
- Test: `test/print-ui/FlagButton.test.tsx`

**Props:**
```ts
export interface FlagButtonProps {
  /** The parsed value, e.g. "1 : 0.25". Rendered mono. */
  value: string;
  /** True → the uncertain-parse look (accent, dashed, "check the read" caption).
   *  False/omitted → the resolved, re-openable look (ink, dotted-on-hover). */
  flagged?: boolean;
  /** Click toggles the read (resolve a flag / re-open a value). */
  onClick?: () => void;
  /** PLACEMENT ONLY. */
  className?: string;
}
```

- [ ] **Step 1 — failing test** (`test/print-ui/FlagButton.test.tsx`):
```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { FlagButton } from "../../src/print/ui/FlagButton";

describe("<FlagButton>", () => {
  it("renders the value and is a button", () => {
    render(<FlagButton value="1 : 0.25" />);
    const btn = screen.getByTestId("flag-button");
    expect(btn.tagName).toBe("BUTTON");
    expect(btn).toHaveTextContent("1 : 0.25");
  });

  it("defaults to the ok state with a re-open title", () => {
    render(<FlagButton value="1 : 1" />);
    const btn = screen.getByTestId("flag-button");
    expect(btn.getAttribute("data-state")).toBe("ok");
    expect(btn.getAttribute("title")).toBe("Click to re-open this value");
    expect(screen.queryByText(/check the read/i)).toBeNull();
  });

  it("renders the flagged state with the check-the-read caption", () => {
    render(<FlagButton value="1 : 0" flagged />);
    const btn = screen.getByTestId("flag-button");
    expect(btn.getAttribute("data-state")).toBe("flagged");
    expect(screen.getByText(/check the read/i)).toBeInTheDocument();
  });

  it("fires onClick", () => {
    const onClick = vi.fn();
    render(<FlagButton value="1 : 0" flagged onClick={onClick} />);
    fireEvent.click(screen.getByTestId("flag-button"));
    expect(onClick).toHaveBeenCalledOnce();
  });
});
```
Run: `npx vitest run test/print-ui/FlagButton.test.tsx` → FAIL (module not found).

- [ ] **Step 2 — implement** `FlagButton.tsx`. `ui/` is guard-exempt, so author appearance directly with tokens (prefer `var(--color-*)`; the dashed/dotted underlines can be `border-b border-dashed`/`border-dotted` with `border-accent`/`border-hair-strong`/`border-transparent`, and the accent-tinted dashed for the flagged underline may use `border-accent` — do NOT chase the mockup's exact `color-mix` blend, the token accent is close enough and on-system). Structure:
```tsx
<button type="button" data-testid="flag-button" data-state={flagged ? "flagged" : "ok"}
  onClick={onClick} title={flagged ? undefined : "Click to re-open this value"}
  className={cx("group/fb block text-right font-mono cursor-pointer min-w-[92px] flex-shrink-0",
    flagged ? "text-accent" : "text-ink", className)}>
  <span className={cx("text-data font-bold border-b pb-px",
    flagged ? "border-dashed border-accent"
            : "border-dotted border-transparent group-hover/fb:border-hair-strong")}>
    {value}
  </span>
  {flagged && (
    <span className="block text-[9px] font-bold uppercase tracking-wide text-accent mt-0.5">
      ▸ check the read
    </span>
  )}
</button>
```
NOTE the `text-[9px]` is permitted here because **`ui/` is exempt** (the caption is genuinely a one-off micro size with no named role). Use `cx` (copy the house helper).
Run the test → PASS.

- [ ] **Step 3 — story** `FlagButton.stories.tsx`: a `Default` (ok), a `Flagged`, and an interactive `Toggling` story (local `useState` flipping `flagged` on click) so the two states + the hover affordance are visible side by side. Title `"ui/FlagButton"`.

- [ ] **Step 4 — gate + commit.**
```
npx vitest run test/print-ui/FlagButton.test.tsx
npm run lint:design
npx tsc --noEmit -p tsconfig.build.json
git add src/print/ui/FlagButton.tsx src/print/ui/FlagButton.stories.tsx src/print/ui/index.ts test/print-ui/FlagButton.test.tsx
git commit -m "feat(print/ui): FlagButton — scoping parsed-value control (flag / re-open)"
```
(commit body ends with the required trailer line)

---

## Task 2: `ScopeSampleRow` (`components/`)

**Mockup region:** `series-scoping.html` `.srow`, `.srow.flagged`, `.grip`, `.spark`, `.s-text`/`.s-name`/`.s-id`, and `rowHTML` (lines 172–200, 452–470).

**The row:** `grip` + 76×28 trace `spark` + (`s-name` 13px/600/ink, `s-id` 10.5px mono ink-faint) flex-grown + the trailing `FlagButton`. When `flagged`, the whole row gets a faint accent wash (`.srow.flagged { background: color-mix(accent 5%) }` → `bg-accent/5`). Bottom hairline divider (the *container* owns last-row divider removal — pass divider as a positional concern; here keep `border-b border-hair` on the row and let `ScopePlate` strip the last one, OR rely on the parent — match the mockup `.srow:last-of-type{border-bottom:none}` via the parent in Task 4; the row itself carries `border-b border-hair`).

**Files:** Create `src/print/components/ScopeSampleRow.tsx` + `.stories.tsx`; Test `test/print-components/ScopeSampleRow.test.tsx`.

**Props:**
```ts
import type { Trace } from "../../api";
export interface ScopeSampleRowProps {
  name: string;
  sampleId: string;
  /** Measured trace for the sparkline. */
  trace: { q: number[]; I: number[] };
  /** Dominant phase → sparkline hue; null/undefined → unindexed. */
  phase?: string | null;
  /** The parsed ordering value, e.g. "1 : 0.25". */
  value: string;
  /** Himalaya is unsure it parsed this right (drives the FlagButton + row wash). */
  flagged?: boolean;
  /** Toggle the flag (resolve / re-open) — forwarded to the FlagButton onClick. */
  onToggleFlag?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}
```

- [ ] **Step 1 — failing test** asserting via data/role/text (NOT classes):
```tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ScopeSampleRow } from "../../src/print/components/ScopeSampleRow";

const TRACE = { q: [0.02, 0.04, 0.08, 0.16], I: [10, 80, 30, 5] };

describe("<ScopeSampleRow>", () => {
  it("renders name, id, sparkline, grip, and the value", () => {
    render(<ScopeSampleRow name="Lipid 1-2 + LL37 1:0.25" sampleId="smp_07"
      trace={TRACE} phase="Pn3m" value="1 : 0.25" />);
    const row = screen.getByTestId("scope-sample-row");
    expect(row).toHaveTextContent("Lipid 1-2 + LL37 1:0.25");
    expect(row).toHaveTextContent("smp_07");
    expect(screen.getByTestId("sparkline")).toBeInTheDocument();
    expect(screen.getByTestId("grip-handle")).toBeInTheDocument();
    expect(screen.getByTestId("flag-button")).toHaveTextContent("1 : 0.25");
  });

  it("marks the flagged state on the row and the value", () => {
    render(<ScopeSampleRow name="Lipid 1-2, no LL37" sampleId="smp_04"
      trace={TRACE} phase="Pn3m" value="1 : 0" flagged />);
    expect(screen.getByTestId("scope-sample-row").getAttribute("data-flagged")).toBe("true");
    expect(screen.getByTestId("flag-button").getAttribute("data-state")).toBe("flagged");
  });

  it("forwards the flag toggle", () => {
    const onToggleFlag = vi.fn();
    render(<ScopeSampleRow name="x" sampleId="smp_04" trace={TRACE} value="1 : 0"
      flagged onToggleFlag={onToggleFlag} />);
    fireEvent.click(screen.getByTestId("flag-button"));
    expect(onToggleFlag).toHaveBeenCalledOnce();
  });
});
```
Run → FAIL (module not found).

- [ ] **Step 2 — implement** `ScopeSampleRow.tsx`. Root `<div data-testid="scope-sample-row" data-flagged={flagged ? "true" : "false"}>` with `group flex items-center gap-3 px-2 py-2.5 border-b border-hair` + `flagged && "bg-accent/5"`. Children in order: `<GripHandle />`, `<Sparkline trace={trace} {...(phase != null ? { phase } : {})} />`, a `<div className="flex-1 min-w-0">` holding name (`text-body font-semibold text-ink`) + id (`text-caption font-mono text-ink-faint`), and `<FlagButton value={value} {...(flagged ? { flagged: true } : {})} {...(onToggleFlag ? { onClick: onToggleFlag } : {})} />`. Use conditional spreads for the optional props (`exactOptionalPropertyTypes: true`). Copy `cx`.
Run test → PASS.

- [ ] **Step 3 — story** using a real fixture trace (import `{ realTraces } from "../fixtures/realTraces"`): an `Ok`, a `Flagged`, a `Coexistence` (phase shown), and a `List` story stacking 3–4 rows so the divider + flagged wash read in context. Title `"components/ScopeSampleRow"`.

- [ ] **Step 4 — gate + commit** (test + `lint:design` + tsc), commit ONLY the 3 named files with the trailer:
`git commit -m "feat(print/components): ScopeSampleRow — scoping series row (grip + spark + flag-value)"`

---

## Task 3: `ScopeCandidateRow` (`components/`)

**Mockup region:** `series-scoping.html` `.candidates`/`.cand-row`/`.cand-info`/`.cand-name`/`.cand-why`/`.cand-why .warn`/`.cand-add`/`.cand-empty` (lines 224–248, `candHTML` 472–481).

**The row:** a **dashed** hairline-bordered, rounded row: dimmed `spark` (opacity-70) + (`cand-name` 12.5px/600/ink-soft, `cand-why` 11px ink-faint — a ReactNode so the caller embeds the accent `.warn` emphasis) flex-grown + a trailing "+ Add to series" `Button variant="outline"` (the `Button` primitive has **no** `size` prop, verified 2026-06-03 — do not pass one).

**Files:** Create `src/print/components/ScopeCandidateRow.tsx` + `.stories.tsx`; Test `test/print-components/ScopeCandidateRow.test.tsx`.

**Props:**
```ts
import type { ReactNode } from "react";
export interface ScopeCandidateRowProps {
  name: string;
  /** One-line rationale; ReactNode so the caller can wrap an accent <strong>. */
  why: ReactNode;
  trace: { q: number[]; I: number[] };
  phase?: string | null;
  onAdd?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}
```

- [ ] **Step 1 — failing test:** asserts `data-testid="scope-candidate-row"`, renders name + why text + the sparkline (dimmed), and an "Add to series" button whose click fires `onAdd`. (Assert the dimming via `data-testid="sparkline"` presence + that the row exists; do not assert the opacity class string.)

- [ ] **Step 2 — implement.** Root `<div data-testid="scope-candidate-row" className="flex items-center gap-3 px-2 py-2.5 border border-dashed border-hair-strong rounded">`. `<Sparkline ... className="opacity-70" />`, `<div className="flex-1 min-w-0">` with name (`text-meta font-semibold text-ink-soft`) + why (`text-caption text-ink-faint mt-px`), `<Button variant="outline" onClick={onAdd}>+ Add to series</Button>` (no `size` prop exists; the outline variant is the load-bearing part). Copy `cx`.

- [ ] **Step 3 — story:** a `Default` (with an accent `<strong className="text-accent font-semibold">` embedded in `why` to show the `.warn` emphasis), using a fixture trace. Title `"components/ScopeCandidateRow"`. Note: the empty-state (`.cand-empty` "Nothing else in the corpus matches this grouping.") is rendered by `ScopePlate`/the page, NOT this row.

- [ ] **Step 4 — gate + commit** (named files only, trailer):
`git commit -m "feat(print/components): ScopeCandidateRow — scoping sample-candidate row (dashed, Add-to-series)"`

---

## Task 4: `ScopePlate` (`components/`) + `ScopingAssembly` story

**Mockup region:** the entire `.scope-plate` (lines 96–286 CSS, 305–366 markup) + the `render()` derivations (482–539).

**The plate:** `Card elevated` worksheet, `max-w-[760px]`, padded. Sections top→bottom:
1. `Kicker tone="accent"` "New series".
2. Serif title (the series name). Render as a `text-display` (snap from the mockup's 29px) serif `<h1>` whose inner text-span carries the dotted-underline editable affordance (`border-b border-dotted border-hair-strong`); editing is page-deferred — accept an optional `onEditTitle` but the input itself is NOT built here.
3. `<AutoGroup variant="summary">` with the `grouping` ReactNode.
4. "Ordered by" — a `Kicker tone="faint"` field label + an `.order-field` button (bordered `border border-hair-strong rounded bg-plate`, hover `border-ink-faint`) showing `orderedBy` (`text-body font-semibold`? mockup 15px → snap to `text-body`/`text-headline`; pick the nearest named role) + a `▾` chevron (`text-ink-faint`) on the right; click → `onChangeOrder` (page opens the menu). Below: the `orderNote` (`text-caption text-ink-faint`).
5. "The series" head: a `Kicker tone="faint"` label on the left; on the right an **Undo** text button (accent, shown only when `onUndo` provided — `undoLabel` as its text/title) + a mono count string (`count`, `text-caption font-mono text-ink-faint`).
6. `rows` slot (the `ScopeSampleRow`×N). Wrap so the **last row's bottom border is removed** (`[&>*:last-child]:border-b-0` is an arbitrary VARIANT selector, not appearance — allowed; or simply document that rows carry the divider and the section clips it). Prefer a `rows` ReactNode slot (children-slotting, the SheetTable/Gallery pattern) — the page maps data→rows.
7. `.candidates` section: `border-t border-hair pt-4`, a `Kicker tone="faint"` "Himalaya also found" header, then the `candidates` ReactNode slot. If the page passes the empty node, it renders the italic `.cand-empty` line — that empty node is the page's responsibility; `ScopePlate` just renders whatever `candidates` is.
8. Preview: `Kicker tone="faint"` "Preview — phase across the series" + `<PhaseStrip segments={preview} size="sm" />`.
9. `.scope-foot`: `border-t border-hair pt-4 flex items-center justify-between gap-5`. Left = a state column: a `fs-line` (a `<Dot tone={footState.kind === "warn" ? "accent" : "success"} />` + `footState.text`, the line tinted accent when warn / ink when ready, `font-semibold`) above a `footNote` (`text-caption text-ink-faint max-w-[42ch]`). Right = `<Button variant="solid" disabled={buildDisabled} onClick={onBuild}>Confirm & build →</Button>`.

**Presentational contract:** NO `useState` in `ScopePlate`. Everything derived (the preview segments, the foot state text, the gate `buildDisabled`, the count) is a PROP computed by the consumer/page.

**Files:** Create `src/print/components/ScopePlate.tsx` + `ScopingAssembly.stories.tsx`; Test `test/print-components/ScopePlate.test.tsx`.

**Props:**
```ts
import type { ReactNode } from "react";
import type { PhaseSegment } from "../ui"; // re-exported from ui/index.ts (verified)
export interface ScopePlateProps {
  seriesName: string;
  onEditTitle?: () => void;
  /** AutoGroup body (embed <strong> for emphasis). */
  grouping: ReactNode;
  orderedBy: string;
  orderNote: ReactNode;
  onChangeOrder?: () => void;
  /** Count line, e.g. "6 samples · low to high". */
  count: string;
  onUndo?: () => void;
  undoLabel?: string;
  /** ScopeSampleRow×N (children-slotting). */
  rows: ReactNode;
  /** ScopeCandidateRow×N, or the page's empty-state node. */
  candidates: ReactNode;
  /** Preview strip segments (low→high). */
  preview: PhaseSegment[];
  footState: { kind: "warn" | "ready"; text: string };
  footNote: ReactNode;
  buildDisabled?: boolean;
  onBuild?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}
```

- [ ] **Step 1 — failing test** (`test/print-components/ScopePlate.test.tsx`) asserting via data/role/text:
  - renders the series name, the "New series" kicker text, the `orderedBy` value + the order note, the count string;
  - renders the `rows` and `candidates` slot content (pass simple `<div data-testid="row-slot" />` / `<div data-testid="cand-slot" />` sentinels and assert they're present);
  - renders a `PhaseStrip` (`ps-seg` present) from a 2-segment `preview`;
  - the **gate**: `buildDisabled` → the Confirm button has the `disabled` attribute; `footState.kind="warn"` shows the warn text;
  - `onBuild` fires on click when not disabled; `onUndo` button only present when `onUndo` is passed.

- [ ] **Step 2 — implement** `ScopePlate.tsx` per the section list. Use conditional spreads for optional handlers. Copy `cx`. Snap off-scale serif/px sizes to the nearest named role (no `text-[Npx]` in `components/`).

- [ ] **Step 3 — `ScopingAssembly.stories.tsx`** — the page simulation (the Batch 7/9/10/11 deferral idiom). It owns the real state:
  - a `SERIES` array (6 members: id, name, value, key, flagged, phase + `coexistWith` for the 2-phase members) and a `CANDIDATES` array (1 member: the `smp_21` 1-1 lipid line with the accent `<strong>` why), each mapped to a `realTraces` fixture trace (reuse the 5 available ids cyclically);
  - `useState` for the SERIES flagged-map + an `undo` history stack (mirror the mockup's flag-toggle + add-candidate + undo);
  - derives `preview: PhaseSegment[]` from each member's phase(+coexist), the foot state (`flags>0 → warn "{n} values to check before you can build"`, else `ready "All N values confirmed — ready to build"`), `buildDisabled = flags>0`, and the count;
  - maps SERIES→`<ScopeSampleRow … onToggleFlag>` (sorted by `key`) into `rows`, CANDIDATES→`<ScopeCandidateRow … onAdd>` (or the `.cand-empty` italic node) into `candidates`;
  - wraps in the centred `bg-paper p-10` stage so it reads as the page. Title `"components/ScopingAssembly"`, default export `meta`, one `Page` story.

- [ ] **Step 4 — gate + commit** (named files only, trailer):
`git commit -m "feat(print/components): ScopePlate + ScopingAssembly — series-scoping worksheet"`

---

## Batch verification (after all 4 tasks)

```
npx vitest run test/print-ui/FlagButton.test.tsx test/print-components/ScopeSampleRow.test.tsx \
  test/print-components/ScopeCandidateRow.test.tsx test/print-components/ScopePlate.test.tsx
npm run lint:design          # clean — proves placement-only
npx tsc --noEmit -p tsconfig.build.json
npm run build-storybook      # exit 0
```
Manual fidelity: `npm run storybook` → `components/ScopingAssembly` side-by-side with `docs/redesign-mockups/series-scoping.html`. Then update `docs/greenfield-component-ledger.md` (a new deferred-registry row + bump the Layer-3 / modals counts; the `ScopePlate`+`SampleRow`(s) entries move from ⬜ to ✅). Ledger update is committed; the plan file is NOT.
```
