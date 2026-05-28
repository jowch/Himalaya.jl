# R3-D — Focus chrome alignment Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Bring the focus workspace's rail / tools / notes / reflections chrome to mockup parity (issue #257, round-3 findings R3-F02/F03/F04/F06/F07/F09/F10 + U-4).

**Architecture:** Surgical, component-local edits across five existing focus components plus a one-line DESIGN.md doc note. No new state fields, no schema work, no styles.css role additions (reuse `text-meta`). All visual findings are asserted via `data-*` attributes and DOM-structure queries, never Tailwind class-string equality. Builds cleanly on top of wave-1 #254's accent-rationing edits already merged on `main` (auto-peak triangles → phase color, `baseColor` precedence, retired "auto peak" legend swatch) — none of those are touched.

**Tech Stack:** React + TypeScript, Vitest + Testing Library + JSDOM, Tailwind (Print token vocabulary), Zustand store.

---

## Branch & verify gate

- Branch: `r3-d-focus-workspace` off `main` (`6d73d64`).
- Worktree dir: `/Users/me/projects/Himalaya.jl/.claude/worktrees/r3-e-series-finish` (directory name is cosmetic; the checked-out branch is `r3-d-focus-workspace` at `6d73d64`).
- Verify gate (run from `packages/HimalayaUI/frontend/`): `npm test` (Vitest, one-shot) AND `npm run build` (tsc --noEmit + vite build) — both must be green before the PR.

## File map

| File | Responsibility / change |
|---|---|
| `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx` (~:352-378) | R3-F02 — Speculative behind `<details>` disclosure; CTA gated on open/has-speculative. |
| `packages/HimalayaUI/frontend/src/components/PlotCard.tsx` (~:474-496, TitleStrip) | R3-F03 — QRange numeric inputs → `Auto-fit` + `+ Peak` ghost buttons. U-4 — zoom indicator in TitleStrip when `xDomain !== null` (off full range). |
| `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx` | R3-F03 — add-peak "armed" mode prop so `+ Peak` arms placement + carries visible state. |
| `packages/HimalayaUI/frontend/src/components/FocusNotesMargin.tsx` (:24-51) | R3-F04 — Print note region: `Notes` head + count badge, quiet body text with mono `q ≈ N.NNN` refs, dashed `✎ Add a note…` affordance. |
| `packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx` (:108-126, :142-174, :224-226) | R3-F09 — sticky `<th>` uses `text-meta` (reuse, no new token). R3-F10 — unindexed `<tr>` dims uniformly via `data-indexed` + opacity. R3-F07 — right-side header N-of-M summary for panel-head symmetry. |
| `DESIGN.md` §5 | R3-F06 — document the `order` (shipped) vs `hkl` (mockup) column-header rename. |

## Two flags for team-lead (resolve at plan-review)

These two are design judgment calls inside the issue's scope; resolve before implementing Task 7.

**Flag 1 — R3-F03 "+ Peak" semantics + the fate of `QNumInput`.**
The issue says `+ Peak` "enters add-peak mode in TraceViewer." But TraceViewer has **no** add-peak mode today — a click on empty plot area *already* adds a manual peak unconditionally (`TraceViewer.tsx:281-283`). To make `+ Peak` meaningful (and match the mockup's `.armed` terracotta state, `focus-workspace.html:177-178`), I propose adding an `addArmed` boolean + `onToggleAddArmed` to TraceViewer: `+ Peak` toggles armed; armed shows the terracotta `.armed` look + a `data-add-armed` flag on the viewer root. Click-to-add still works whether or not armed (NO regression to existing behavior — armed only adds an affordance + visible state, it does not gate placement). Separately, `QNumInput` is **exported** and has 5 unit tests (`PlotCard.test.tsx`). The issue offers a fallback: "fold numeric editing into a popover off the segmented control." I propose **keeping `QNumInput` exported but removing the always-visible `<QRange>` cluster from TitleStrip** — the tested unit survives untouched, the numeric pair stops reading as toolbar-soup, and a future popover can reuse it. **Need team-lead's OK on: (a) the armed-mode approach (vs a true placement-gating mode, which is a larger change), (b) keeping QNumInput exported-but-unmounted vs deleting it + its 5 tests.**

**Flag 2 — R3-F09 `text-meta` is 11.5px, not 10.5px.**
The issue/findings claim `text-meta` is "10.5px — 1px delta invisible." Verified in `styles.css:87,268`: `text-meta` is `--text-sm` = **11.5px**. Current `<th>` is `text-[9.5px]`. The real jump is 9.5 → 11.5 (**2px**), not 1px. Per team-lead's explicit directive I am **taking the reuse path (use `text-meta`), NOT adding a `text-th-caps` role to styles.css.** Flagging only so the 2px header bump is a known, accepted consequence rather than a review surprise. (Confirmed: styles.css will not be touched.)

---

## Task 1: R3-F09 — Reflections sticky `<th>` uses the `text-meta` scale token

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx:108-126`
- Test: `packages/HimalayaUI/frontend/test/FocusReflectionsTable.test.tsx`

- [ ] **Step 1: Write the failing test** — assert the sticky header cells carry the scale-token role, not an inline pixel size. Add inside the existing `describe("FocusReflectionsTable (#209)")`:

```tsx
it("R3-F09: column headers use the text-meta scale token, not an inline px size", async () => {
  renderTable();
  await screen.findByTestId("focus-reflections");
  const ths = document.querySelectorAll('[data-testid="focus-reflections"] thead th');
  expect(ths.length).toBe(4);
  for (const th of ths) {
    // Fixed-Scale Rule: no inline text-[Npx] on the sticky header.
    expect(th.className).not.toMatch(/text-\[\d/);
    expect(th.className).toContain("text-meta");
  }
});
```

- [ ] **Step 2: Run test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`): `npm test -- FocusReflectionsTable`
Expected: FAIL — current `<th>` className contains `text-[9.5px]` and lacks `text-meta`.

- [ ] **Step 3: Implement** — in each of the four `<th>` classNames (`:108-126`), replace `text-[9.5px] font-bold uppercase tracking-[0.07em]` with `text-meta uppercase tracking-[0.07em] font-bold`. (`text-meta` carries `font-weight:500`; keep the explicit `font-bold` after it so bold wins. Drop only the `text-[9.5px]` literal.) Apply to all four header cells; the two numeric cells keep their trailing `font-mono`.

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- FocusReflectionsTable`
Expected: PASS (all existing tests + the new one).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx packages/HimalayaUI/frontend/test/FocusReflectionsTable.test.tsx
git commit -m "R3-F09: reflections th uses text-meta scale token (no inline px)"
```

## Task 2: R3-F10 — Unindexed reflections row dims uniformly

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx:142-174`
- Test: `packages/HimalayaUI/frontend/test/FocusReflectionsTable.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
it("R3-F10: unindexed row is flagged data-indexed='false' and dimmed as a whole", async () => {
  renderTable();
  const row3 = await screen.findByTestId("reflection-row-3"); // unclaimed peak
  const row1 = screen.getByTestId("reflection-row-1");        // Pn3m-claimed
  expect(row3).toHaveAttribute("data-indexed", "false");
  expect(row1).toHaveAttribute("data-indexed", "true");
  // Whole-row dim (mockup `.refl tr.unindexed`), not just the phase cell.
  expect(row3.className).toMatch(/opacity-55/);
  expect(row1.className).not.toMatch(/opacity-55/);
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- FocusReflectionsTable`
Expected: FAIL — rows have no `data-indexed` attribute and no `opacity-55`.

- [ ] **Step 3: Implement** — on the `<tr>` (`:142-161`), add `data-indexed={indexed ? "true" : "false"}` and append `opacity-55` to the className when `!indexed`. `indexed` is already in scope at `:133`. Compute the className:

```tsx
className={[
  "cursor-pointer border-b border-hair last:border-b-0 transition-colors",
  indexed ? "" : "opacity-55",
].join(" ")}
data-indexed={indexed ? "true" : "false"}
```

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- FocusReflectionsTable`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx packages/HimalayaUI/frontend/test/FocusReflectionsTable.test.tsx
git commit -m "R3-F10: unindexed reflections row dims as a whole (opacity-55)"
```

## Task 3: R3-F07 — Reflections panel-head right-side cluster

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx:224-226`
- Test: `packages/HimalayaUI/frontend/test/FocusReflectionsTable.test.tsx`

The mockup detector panel-head is `[label | exposure cluster]`; the reflections panel-head ships label-only. Add a right-side `N / M` summary (lifted from the footer's `covered` / `peaks.length`, already computed at `:84`) so the two lower-row panels read symmetric.

- [ ] **Step 1: Write the failing test**

```tsx
it("R3-F07: reflections panel-head carries a right-side N-of-M summary cluster", async () => {
  renderTable();
  await screen.findByTestId("focus-reflections");
  const cluster = screen.getByTestId("focus-reflections-head-summary");
  // 2 of 3 peaks claimed in the fixture.
  expect(cluster).toHaveTextContent("2 / 3");
});
```

- [ ] **Step 2: Run test to verify it fails**

Run: `npm test -- FocusReflectionsTable`
Expected: FAIL — no `focus-reflections-head-summary` testid.

- [ ] **Step 3: Implement** — the panel-head `<div>` at `:224-226` already uses `flex items-center justify-between gap-3`. Add the summary span after the `Reflections` label:

```tsx
<div className="card-header flex items-center justify-between gap-3">
  <span className="text-meta uppercase tracking-wider">Reflections</span>
  <span
    data-testid="focus-reflections-head-summary"
    className="text-meta tabular-nums text-ink-faint"
  >
    {covered} / {peaks.length}
  </span>
</div>
```

Keep it unconditional for symmetry; at 0/0 it reads fine and mirrors the detector head always showing its cluster (the body still shows the empty hint).

- [ ] **Step 4: Run test to verify it passes**

Run: `npm test -- FocusReflectionsTable`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/FocusReflectionsTable.tsx packages/HimalayaUI/frontend/test/FocusReflectionsTable.test.tsx
git commit -m "R3-F07: reflections panel-head N-of-M summary (detector-panel symmetry)"
```

## Task 4: R3-F06 — Document the `order` vs `hkl` column rename in DESIGN.md §5

**Files:**
- Modify: `DESIGN.md` §5 (Signature components, after the `Phase chip` bullet ~:226)

No test (doc-only). The shipped header is `order` (peaks are ratio-series positions, not Miller triples); the mockup says `hkl`. Document the deliberate rename.

- [ ] **Step 1: Add the doc note** — append a bullet under §5 Signature components:

```markdown
- **Reflections card:** the lower-row table beside the detector window. Column headers are `phase · order · q · d`. The mockup labels the second column `hkl`; the shipped header is **`order`** by deliberate rename — SAXS peaks in our data model carry a 1-based ratio-series position, not an (h,k,l) Miller triple (`lib/seriesRatio.ts`, `src/phase.jl phaseratios`). "order" reads as "1st reflection, 2nd reflection…", which is what the value is. Treat the live `order` header as correct; the mockup `hkl` is the divergence.
```

- [ ] **Step 2: Verify it renders**

Run: `grep -n "order.*hkl\|Reflections card" DESIGN.md`
Expected: the new bullet appears in §5.

- [ ] **Step 3: Commit**

```bash
git add DESIGN.md
git commit -m "R3-F06: document order vs hkl reflections column rename (DESIGN §5)"
```

## Task 5: R3-F04 — Notes margin Print treatment

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/FocusNotesMargin.tsx:24-51`
- Test: `packages/HimalayaUI/frontend/test/FocusNotesMargin.test.tsx`

Stage 1 only (no schema). Render a `Notes` head with a count badge (`1` when `sample.notes` non-empty, hidden when empty); the saved note as quiet body text with `q ≈ N.NNN` substrings wrapped in a mono span; a dashed `✎ Add a note…` affordance wrapping the editor. The focus-gate sync (`focused` ref) and `onSaveNotes`-on-blur contract MUST be preserved — the two existing tests depend on `focus-notes-input` staying the edit target with the same draft-clobber guard, and `getByDisplayValue` targeting it.

- [ ] **Step 1: Write the failing tests** — add to `describe("FocusNotesMargin")`:

```tsx
it("R3-F04: shows a count badge when notes are present", () => {
  render(<FocusNotesMargin sample={SAMPLE} onSaveNotes={() => {}} />);
  expect(screen.getByTestId("focus-notes-count")).toHaveTextContent("1");
});

it("R3-F04: hides the count badge when notes are empty", () => {
  render(<FocusNotesMargin sample={{ ...SAMPLE, notes: "" }} onSaveNotes={() => {}} />);
  expect(screen.queryByTestId("focus-notes-count")).toBeNull();
});

it("R3-F04: mono-formats q ≈ N.NNN references in the note body", () => {
  render(<FocusNotesMargin sample={{ ...SAMPLE, notes: "weak shoulder at q ≈ 0.064 here" }} onSaveNotes={() => {}} />);
  const refs = screen.getAllByTestId("focus-notes-qref");
  expect(refs.length).toBeGreaterThanOrEqual(1);
  expect(refs[0]).toHaveTextContent("q ≈ 0.064");
});

it("R3-F04: renders the dashed add-a-note affordance", () => {
  render(<FocusNotesMargin sample={{ ...SAMPLE, notes: "" }} onSaveNotes={() => {}} />);
  expect(screen.getByTestId("focus-notes-add")).toBeInTheDocument();
});
```

(NOTE the existing fixture `SAMPLE.notes` is `"weak shoulder at q≈0.064"` — no spaces around `≈`. The qref test uses its OWN sample with spaced `q ≈ 0.064` to match the mockup form; the regex below accepts optional whitespace so both render mono.)

- [ ] **Step 2: Run tests to verify they fail**

Run: `npm test -- FocusNotesMargin`
Expected: FAIL — no count / qref / add testids exist.

- [ ] **Step 3: Implement** — rewrite the component. Add a `renderNoteBody` helper splitting on a `q ≈ N.NNN` regex (optional whitespace around `≈`) and wrapping matches in mono spans. Keep the focus-gated textarea editor.

```tsx
import { useEffect, useState } from "react";
import type { Sample } from "../api";

interface Props {
  sample: Sample;
  /** Called on blur with the edited notes value. */
  onSaveNotes: (notes: string) => void;
}

// Matches the mockup's mono `q ≈ 0.064` reference (whitespace around ≈ optional,
// so the existing "q≈0.064" fixture also formats).
const Q_REF = /q\s*≈\s*\d+(?:\.\d+)?/g;

function renderNoteBody(text: string): JSX.Element[] {
  const out: JSX.Element[] = [];
  let last = 0;
  let key = 0;
  let m: RegExpExecArray | null;
  Q_REF.lastIndex = 0;
  while ((m = Q_REF.exec(text)) !== null) {
    if (m.index > last) out.push(<span key={key++}>{text.slice(last, m.index)}</span>);
    out.push(
      <span key={key++} data-testid="focus-notes-qref"
            className="font-mono font-semibold text-ink-soft">{m[0]}</span>,
    );
    last = m.index + m[0].length;
  }
  if (last < text.length) out.push(<span key={key++}>{text.slice(last)}</span>);
  return out;
}

/**
 * FocusNotesMargin — the Notes margin on the focus workspace (mockup `.notes-margin`:
 * a quiet ruled column, deliberately not a card). Stage-1 Print treatment (R3-F04):
 * a `Notes` head + count badge, the saved note rendered as quiet body text with
 * mono-formatted `q ≈ N.NNN` references, and a dashed `✎ Add a note…` affordance.
 * Threaded per-note metadata (author + timestamp) is stage 2, deferred until the
 * schema supports it.
 *
 * Focus-gated edit (mirrors SampleMetadataCard / QNumInput): the external
 * `sample.notes` is synced into the draft only while the textarea is NOT focused,
 * so a background refetch or another user's save can't clobber a mid-edit draft.
 */
export function FocusNotesMargin({ sample, onSaveNotes }: Props): JSX.Element {
  const [draft, setDraft] = useState(sample.notes ?? "");
  const [focused, setFocused] = useState(false);

  useEffect(() => {
    if (!focused) setDraft(sample.notes ?? "");
  }, [sample.id, sample.notes, focused]);

  const saved = (sample.notes ?? "").trim();
  const hasNote = saved.length > 0;

  return (
    <aside data-testid="focus-notes-margin"
           className="flex flex-col gap-[15px] overflow-y-auto border-l border-hair bg-paper px-[19px] py-[22px]">
      <div className="flex items-baseline gap-[7px]">
        <span className="text-meta uppercase tracking-wider text-ink-faint">Notes</span>
        {hasNote && (
          <span data-testid="focus-notes-count"
                className="font-mono text-[10px] font-bold text-ink-faint">1</span>
        )}
      </div>

      {hasNote && (
        <div data-testid="focus-notes-body"
             className="text-xs leading-[1.6] text-ink-soft">
          {renderNoteBody(saved)}
        </div>
      )}

      <div data-testid="focus-notes-add"
           className="mt-auto border-b border-dashed border-hair-strong pb-[7px] pt-[6px]
                      text-xs text-ink-faint before:text-[var(--color-accent)] before:content-['✎_']">
        <textarea
          data-testid="focus-notes-input"
          value={draft}
          onChange={(e) => setDraft(e.target.value)}
          onFocus={() => setFocused(true)}
          onBlur={() => { setFocused(false); onSaveNotes(draft); }}
          placeholder="Add a note…"
          className="min-h-[44px] w-full resize-none bg-transparent text-xs
                     text-ink-soft outline-none placeholder:text-ink-faint"
        />
      </div>
    </aside>
  );
}
```

(`✎ ` prefix via Tailwind `before:content` arbitrary value matches mockup `.nm-add::before`. The textarea stays the editor so the existing focus-gate + blur + `getByDisplayValue` tests pass unchanged.)

- [ ] **Step 4: Run tests to verify all pass**

Run: `npm test -- FocusNotesMargin`
Expected: PASS — the three original tests (`renders the sample's notes` via `getByDisplayValue`, draft-clobber guard, `onSaveNotes`-on-blur) PLUS the four new ones.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/FocusNotesMargin.tsx packages/HimalayaUI/frontend/test/FocusNotesMargin.test.tsx
git commit -m "R3-F04: Print notes margin — count badge, mono q-refs, dashed add affordance"
```

## Task 6: R3-F02 — Speculative section behind a disclosure

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PhasePanel.tsx:352-378`
- Test: `packages/HimalayaUI/frontend/test/PhasePanel.test.tsx`

Pull Speculative below-the-fold behind a `<details>` disclosure so the rail reads as two calm sections (Phase call + Candidates). `<details open={speculatives.length > 0}>` auto-opens when the user has built speculatives (so they don't vanish) but is collapsed by default. The CTA lives inside the `<details>` body, so it only shows when open or when speculatives force-open it — satisfying "render the CTA only on hover/focus or when at least one speculative exists."

- [ ] **Step 1: Confirm scope** — `speculatives` (`:251`), `memberIds`, `toggle`, `setHoveredIndex`, `deleteIndex` (`:231`), `openBuilder` (`:233`), and the `exposureId` PROP (`:222`) are all in scope at `:352` (used in the current block). PhasePanel takes `exposureId` as a prop; tests render `renderWithProviders(<PhasePanel exposureId={42} />)` with a per-test `vi.stubGlobal("fetch", …)` returning `[]` for `/indices` + `/groups` (mirror `PhasePanel.test.tsx:90-103`). There is NO shared `renderPhasePanel` helper.

- [ ] **Step 2: Write the failing test** — add to `PhasePanel.test.tsx`, mirroring the sibling candidate-multiselect tests' fetch-stub + `renderWithProviders` setup with a fixture that has ZERO speculative indices (i.e. no index with `kind: "speculative"`):

```tsx
it("R3-F02: Speculative lives behind a collapsed disclosure when none exist", async () => {
  // stub fetch so /indices and /groups settle to empty (no speculative indices),
  // exactly like the sibling tests in this file, then:
  renderWithProviders(<PhasePanel exposureId={42} />);
  const disclosure = await screen.findByTestId("speculative-disclosure");
  // Contract: the Speculative section is a NATIVE <details> that is collapsed
  // by default (the browser hides the body, incl. the CTA, until opened).
  expect(disclosure.tagName.toLowerCase()).toBe("details");
  expect(disclosure).not.toHaveAttribute("open");
  // The CTA is a DESCENDANT of the (closed) disclosure — i.e. behind-the-fold,
  // not a sibling that renders unconditionally. NOTE: JSDOM does not honor the
  // native display:none of closed <details> content, so `queryByTestId` would
  // still FIND the CTA — assert containment, not absence.
  const cta = screen.getByTestId("add-speculative-button");
  expect(disclosure).toContainElement(cta);
});

it("R3-F02: disclosure is open when speculatives exist (user builds stay visible)", async () => {
  // stub fetch so /indices returns one index with kind: "speculative" and a
  // /groups payload that does NOT include it; mirror the sibling fixture shape.
  renderWithProviders(<PhasePanel exposureId={42} />);
  const disclosure = await screen.findByTestId("speculative-disclosure");
  expect(disclosure).toHaveAttribute("open");
});
```

(Inline the same `vi.stubGlobal("fetch", …)` wiring a neighboring test uses; the file does per-test fetch stubs, not a shared helper. The second test pins the auto-open-when-non-empty behavior — `open={speculatives.length > 0}`.)

- [ ] **Step 3: Run test to verify it fails**

Run: `npm test -- PhasePanel`
Expected: FAIL — current Speculative block is a plain `<div>` and the CTA always renders.

- [ ] **Step 4: Implement** — replace the Speculative `<div className="flex flex-col gap-2.5">…</div>` block (`:353-378`) with a `<details>` disclosure, keeping `RailHead` inside the `<summary>`:

```tsx
{/* Speculative — user-built sub-minpeaks indices. R3-F02: below-the-fold
    behind a disclosure so the rail reads as the calm two-section output
    (Phase call + Candidates). Auto-open when speculatives exist so the
    user's own builds stay visible. */}
<details
  data-testid="speculative-disclosure"
  open={speculatives.length > 0}
  className="flex flex-col gap-2.5"
>
  <summary className="cursor-pointer list-none">
    <RailHead>Speculative</RailHead>
  </summary>
  {speculatives.length > 0 && (
    <div className="mt-2.5 flex flex-col gap-[7px]">
      {speculatives.map((ix) => (
        <CandidateRow
          key={ix.id}
          index={ix}
          inCall={memberIds.has(ix.id)}
          onToggle={() => toggle(ix)}
          onHover={() => setHoveredIndex(ix.id)}
          onLeave={() => setHoveredIndex(undefined)}
          onDelete={() => deleteIndex.mutate(ix.id)}
        />
      ))}
    </div>
  )}
  <button
    type="button"
    data-testid="add-speculative-button"
    className="mt-2.5 w-full text-xs text-ink-faint border border-dashed border-hair rounded-md py-1.5 hover:text-ink hover:bg-paper-sunk transition-colors"
    onClick={() => openBuilder(exposureId)}
  >
    + Add speculative
  </button>
</details>
```

(`<details>` natively hides body content when closed, so the closed-no-speculative state shows only the `summary` head — the two-section mockup. `list-none` drops the default disclosure triangle so the rail head reads as before.)

- [ ] **Step 5: Run test to verify it passes + no regressions**

Run: `npm test -- PhasePanel`
Expected: PASS — new test + all existing PhasePanel tests.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PhasePanel.tsx packages/HimalayaUI/frontend/test/PhasePanel.test.tsx
git commit -m "R3-F02: Speculative section behind a collapsed disclosure (calm two-section rail)"
```

## Task 7: R3-F03 — Tools cluster: Auto-fit + `+ Peak` ghost buttons; retire QRange inputs

**BLOCKED on team-lead's ruling (Flag 1).** Plan below assumes the proposed approach: add an armed mode; remove the visible QRange cluster; keep `QNumInput` exported.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/TraceViewer.tsx` (add `addArmed` / `onToggleAddArmed`; viewer-root flag)
- Modify: `packages/HimalayaUI/frontend/src/components/PlotCard.tsx` (TitleStrip tools cluster `:474-496`; `addArmed` state; pass-through)
- Test: `packages/HimalayaUI/frontend/test/PlotCard.headerSlot.test.tsx`, `packages/HimalayaUI/frontend/test/TraceViewer.test.tsx`

- [ ] **Step 1: TraceViewer — add armed-mode props (test first).** The TraceViewer test renders the component directly with a `defaultProps` spread (`TraceViewer.test.tsx:110-135`) and asserts on the existing `data-testid="trace-viewer"` root (`TraceViewer.tsx:754`) via `container.querySelector`. Mirror that existing render block:

```tsx
it("R3-F03: armed mode is reflected on the viewer root via data-add-armed", () => {
  const trace = { q: [0.1, 0.2, 0.3], I: [10, 20, 30], sigma: [1, 1, 1] };
  const { container } = render(
    <TraceViewer
      trace={trace}
      peaks={[]}
      activeGroupIndices={[]}
      hoveredIndex={undefined}
      addArmed={true}
      {...defaultProps}
    />,
  );
  expect(
    container.querySelector('[data-testid="trace-viewer"]'),
  ).toHaveAttribute("data-add-armed", "true");
});
```

- [ ] **Step 2: Run — fails** (`npm test -- TraceViewer`): `addArmed` is not a valid prop and the root has no `data-add-armed`.

- [ ] **Step 3: Implement TraceViewer** — add optional props to `TraceViewerProps`:

```tsx
  /** R3-F03: when true, the trace plate is in "place a peak" mode — the tools
   *  cluster's `+ Peak` ghost is armed. Click-to-add empty-area placement works
   *  regardless; armed only carries the visible state + a cursor cue. */
  addArmed?: boolean;
  onToggleAddArmed?: () => void;
```

Destructure `addArmed` in the component signature. The host root at `TraceViewer.tsx:751-756` ALREADY carries `data-testid="trace-viewer"`; add `data-add-armed={addArmed ? "true" : "false"}` to that same `<div>`. No change to the click handler — empty-area click still adds a peak; arming is presentational only. Keeps the change surgical and #254-compatible. (`onToggleAddArmed` is accepted for symmetry but the toggle owner is PlotCard; TraceViewer needs only `addArmed` to render the state — leave `onToggleAddArmed` optional/unused here unless a future click-cue uses it.)

- [ ] **Step 4: Run — passes** (`npm test -- TraceViewer`).

- [ ] **Step 5: PlotCard — tools cluster (test first).** Add to `PlotCard.headerSlot.test.tsx`:

```tsx
it("R3-F03: tools cluster shows Auto-fit + Add Peak ghost buttons, no numeric q-range inputs", () => {
  renderWithProviders(
    <PlotCard headerSlot={<div data-testid="custom-header">Custom</div>} />,
  );
  expect(screen.getByTestId("tool-autofit")).toBeInTheDocument();
  expect(screen.getByTestId("tool-add-peak")).toBeInTheDocument();
  // The numeric q-range editing pair is retired from the visible toolbar.
  expect(screen.queryByTestId("q-range-controls")).toBeNull();
  expect(screen.queryByTestId("q-range-min")).toBeNull();
});
```

- [ ] **Step 6: Run — fails** (`npm test -- PlotCard.headerSlot`).

- [ ] **Step 7: Implement PlotCard.** Add `const [addArmed, setAddArmed] = useState(false);`; reset it in the existing `useEffect([activeExposureId])` reset alongside `setXDomain(null)`. Pass `addArmed` + `onToggleAddArmed={() => setAddArmed((v) => !v)}` into the focus-body `<TraceViewer>`. In `TitleStrip`:
  - Rename the "fit features" button copy to `Auto-fit`, set `data-testid="tool-autofit"`, ghost `tool-btn` styling.
  - Add a `+ Peak` ghost button `data-testid="tool-add-peak"` calling a new `onToggleAddPeak` prop; armed adds `data-armed` + the terracotta look.
  - Remove the `<QRange .../>` element. Keep `XScaleToggle` + the export cluster + divider.
  - Thread `addArmed: boolean` + `onToggleAddPeak: () => void` through `TitleStripProps`; wire them in PlotCard's `<TitleStrip>` render.

  Tools cluster becomes:

```tsx
<div className="shrink-0 flex items-center gap-2">
  <ZoomIndicator zoomed={isZoomed} onReset={onFitFeatures} />
  <XScaleToggle xType={xType} onSetXType={onSetXType} />
  <button
    type="button"
    onClick={onFitFeatures}
    disabled={!canFit}
    data-testid="tool-autofit"
    title="Auto-zoom to peaks (or post-beam region)"
    className="rounded-md border border-hair-strong bg-plate px-2.5 py-1
               text-[11px] font-semibold text-ink hover:bg-paper-sunk
               disabled:opacity-40 disabled:cursor-default whitespace-nowrap"
  >
    Auto-fit
  </button>
  <button
    type="button"
    onClick={onToggleAddPeak}
    data-testid="tool-add-peak"
    data-armed={addArmed ? "true" : "false"}
    aria-pressed={addArmed}
    title="Click the trace to place a peak"
    className={[
      "rounded-md border px-2.5 py-1 text-[11px] font-semibold whitespace-nowrap transition-colors",
      addArmed
        ? "bg-print-accent border-print-accent text-paper"
        : "border-hair-strong bg-plate text-ink hover:bg-paper-sunk",
    ].join(" ")}
  >
    + Peak
  </button>
  <span className="w-px h-4 bg-hair-strong" aria-hidden="true" />
  <FigureExportControls
    spec={exportSpec}
    filenameStem={exportFilenameStem}
    ariaContext="trace plot"
    disabled={exportDisabled}
  />
</div>
```

(`ZoomIndicator` arrives in Task 8 — if Task 8 is done first, it is already present; if Task 7 is done first, render it then. The inline `text-[11px]` MATCHES the mockup's `.tool-btn`/`.seg button` size (`focus-workspace.html:165,172`); it is the mockup-specified control scale, not a Fixed-Scale §3 one-off. Fallback if the reviewer disagrees: `text-xs` to match `XScaleToggle`. Call out in the PR description.)

- [ ] **Step 8: Run — passes** (`npm test -- PlotCard`). Confirm the 5 `QNumInput` tests in `PlotCard.test.tsx` still pass (QNumInput stays exported; only its mounting in TitleStrip is removed).

- [ ] **Step 9: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/TraceViewer.tsx packages/HimalayaUI/frontend/src/components/PlotCard.tsx packages/HimalayaUI/frontend/test/PlotCard.headerSlot.test.tsx packages/HimalayaUI/frontend/test/TraceViewer.test.tsx
git commit -m "R3-F03: tools cluster Auto-fit + Peak ghosts; retire numeric q-range from toolbar"
```

## Task 8: U-4 — Zoom-state indicator in the TitleStrip

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PlotCard.tsx` (export `ZoomIndicator`; TitleStrip wiring; `isZoomed` prop)
- Test: `packages/HimalayaUI/frontend/test/PlotCard.headerSlot.test.tsx`

Show a small `text-meta` terracotta ghost-button ("zoomed · reset") when the visible x-domain is off the full range; clicking it calls `onFitFeatures` (Auto-fit), pairing with Task 7. "Off-default" = PlotCard's own `xDomain` state is non-null (null = full trace). Extract a pure `ZoomIndicator` so it can be unit-tested without driving a wheel event through JSDOM.

- [ ] **Step 1: Write the failing tests** — add to `PlotCard.headerSlot.test.tsx`:

```tsx
it("U-4: no zoom indicator when the trace is at full range (default)", () => {
  renderWithProviders(
    <PlotCard headerSlot={<div data-testid="custom-header">Custom</div>} />,
  );
  expect(screen.queryByTestId("zoom-indicator")).toBeNull();
});
```

Plus a direct unit on the extracted component (import `ZoomIndicator` from `../src/components/PlotCard`):

```tsx
import { ZoomIndicator } from "../src/components/PlotCard";

it("U-4: ZoomIndicator renders a reset ghost button only when zoomed", async () => {
  const onReset = vi.fn();
  const { rerender } = render(<ZoomIndicator zoomed={false} onReset={onReset} />);
  expect(screen.queryByTestId("zoom-indicator")).toBeNull();
  rerender(<ZoomIndicator zoomed={true} onReset={onReset} />);
  const btn = screen.getByTestId("zoom-indicator");
  expect(btn).toHaveTextContent(/zoomed/i);
  await userEvent.setup().click(btn);
  expect(onReset).toHaveBeenCalled();
});
```

(Add `import { render } from "@testing-library/react"` / `userEvent` if not already imported in this file.)

- [ ] **Step 2: Run — fails** (`npm test -- PlotCard.headerSlot`): no `ZoomIndicator` export / `zoom-indicator` testid.

- [ ] **Step 3: Implement** — add an exported pure component near the other PlotCard sub-components:

```tsx
export function ZoomIndicator(
  { zoomed, onReset }: { zoomed: boolean; onReset: () => void },
): JSX.Element | null {
  if (!zoomed) return null;
  return (
    <button
      type="button"
      data-testid="zoom-indicator"
      onClick={onReset}
      title="Reset to full q-range"
      className="text-meta text-print-accent hover:underline whitespace-nowrap"
    >
      zoomed · reset
    </button>
  );
}
```

In PlotCard, pass `isZoomed={xDomain !== null}` to `<TitleStrip>`; add `isZoomed: boolean` to `TitleStripProps`; render `<ZoomIndicator zoomed={isZoomed} onReset={onFitFeatures} />` at the start of the tools `<div>` (as shown in Task 7's cluster). Reset → `onFitFeatures` per the issue ("ghost button that calls `fitFeatures`"; pairs with Auto-fit).

- [ ] **Step 4: Run — passes** (`npm test -- PlotCard`).

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PlotCard.tsx packages/HimalayaUI/frontend/test/PlotCard.headerSlot.test.tsx
git commit -m "U-4: zoom-state indicator in TitleStrip (zoomed · reset → Auto-fit)"
```

## Task 9: Full verify gate

- [ ] **Step 1: Full Vitest run**

Run (from `packages/HimalayaUI/frontend/`): `npm test`
Expected: all suites green (no regressions in PlotCard / TraceViewer / FocusNotesMargin / FocusReflectionsTable / PhasePanel / headerSlot).

- [ ] **Step 2: Build**

Run: `npm run build`
Expected: `tsc --noEmit` clean + vite build succeeds.

- [ ] **Step 3: No NEW dark-era tokens** (R3-F05 is out of scope; just don't ADD relics)

Run: `git diff main -- packages/HimalayaUI/frontend/src/components | grep -E '^\+' | grep -E 'bg-bg|text-fg|border-border|bg-bg-subtle|bg-bg-hover'`
Expected: empty (used Print tokens: `text-ink`, `text-ink-faint`, `text-ink-soft`, `bg-plate`, `bg-paper-sunk`, `border-hair`, `border-hair-strong`, `bg-print-accent`, `text-print-accent`, `text-paper`, `text-meta`).

- [ ] **Step 4: styles.css untouched** (reuse path held)

Run: `git diff main --stat -- packages/HimalayaUI/frontend/src/styles.css`
Expected: no output.

## Self-review — spec coverage

| Done-when (issue #257) | Task |
|---|---|
| Rail two sections; Speculative below-fold/disclosure; CTA gated | Task 6 (R3-F02) |
| Tools cluster: Auto-fit + `+ Peak`; numeric q-range gone | Task 7 (R3-F03) |
| Notes margin Print-styled: count badge + mono q-refs + dashed add; textarea-first gone | Task 5 (R3-F04) |
| Zoom indicator in TitleStrip when off-default; click auto-fits | Task 8 (U-4) |
| Reflections sticky `<th>` uses a scale token | Task 1 (R3-F09) |
| Unindexed rows dim uniformly | Task 2 (R3-F10) |
| Reflections panel header right-side cluster (detector symmetry) | Task 3 (R3-F07) |
| DESIGN.md §5 documents the `order` column rename | Task 4 (R3-F06) |

All eight findings mapped. No styles.css role added (reuse `text-meta`). All visual assertions use `data-*` / DOM structure / scale-token-contract regexes, never Tailwind class-string equality on styling. Build + test gate in Task 9. Tasks are independent except Task 7 ↔ Task 8 share the TitleStrip tools `<div>` (do Task 8 first or merge their cluster edit); Task 7 is gated on Flag 1.
