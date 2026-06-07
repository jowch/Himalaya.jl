# Greenfield Batch 11 — Contact-sheet slice (`SheetTable` + `CullBar`)

> Execute with superpowers:subagent-driven-development (fresh implementer per task + two-stage review).
> Commit trailer (every commit, exact last line): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Worktree frontend dir: `packages/HimalayaUI/frontend`. Typecheck: `npx tsc --noEmit -p tsconfig.build.json`.
> Commit ONLY named files. NEVER `git add -A`/`.`. NEVER stage `src/bones/*` or anything under `docs/superpowers/plans/`.

## Context

The greenfield "The Print" rebuild (branch `worktree-greenfield-ui-rebuild`, NOT merged). Batch 10 (Loupe) is done. This slice builds the **contact-sheet** surface's two remaining Layer-3 composites, derived from `docs/redesign-mockups/sample-table.html`:

- **`SheetTable`** — the `.sheet` plate: aligned column header + a rows region that slots `SampleTableRow`×N.
- **`CullBar`** — the floating dark batch-reject bar (`.cullbar`) that appears on multi-frame select.

Every row leaf already exists (`SampleTableRow` ✅, which exports `SAMPLE_TABLE_COLS`). `CullBar` trips **one refactor-on-contact**: the `Button` primitive has no light-on-dark ghost variant, so the bar's Restore/Clear actions need a new `ghostInverse` variant added in the exempt `ui/` layer (Task 1). The Reject action maps to the existing `variant="accent"`.

**Layer-4 deferral (consistent with Batch 7/9/10):** the `.sheet-shell` page chrome — the `.head` (Contact-sheet `Kicker` + sub + the `6/9 samples screened` `ProgressBar`), the `.kb-legend`, the loupe view, and the page's selection-Set state — is page assembly. A combined `ContactSheetAssembly` story (Task 4) simulates the page for fidelity without building the page component.

### Verified leaf APIs (checked against live source — do NOT re-derive)

- **`SampleTableRow`** (`./SampleTableRow`) — props `{ name; sampleId; screened?; exposures: GalleryExposure[]; selectedExposureId?; onSelectExposure?; kept; total; dropped?; tags: Tag[]; onAddTag?; onRemoveTag?; phase?: string|null; className? }`. `data-testid="sample-table-row"`, `data-screened`. **Exports `SAMPLE_TABLE_COLS`** (the 5-track grid template string) — the header MUST reuse this exact constant so columns align.
- **`Card`** (`../ui`) — `elevated` → the `.card` rule (bg-plate + hairline + radius + soft warm shadow). Default no padding. Pass `className="overflow-hidden"` so the header border + rows clip to the rounded corners.
- **`Kicker`** (`../ui`) — `{ tone?: "accent"|"faint"; as?: "div"|"span"|"h2"|"h3" }`. `tone="faint"` = `text-kicker text-kicker-faint` (uppercase/700/tracked/ink-faint) — exactly the `.sheet-head` column-label role.
- **`Button`** (`../ui`) — `variant?: "solid"|"accent"|"ghost"|"danger"|"outline"` (+ `ghostInverse` after Task 1); `armed?`. `data-variant` reflects the variant. Renders a `<button>`; forwards `...props` (incl. `onClick`, `aria-label`).
- **`ProgressBar`** (`../ui`) — `{ value; total; label?; className? }`, `data-testid="progressbar"`, `role="progressbar"`. (assembly only)
- **`cx` helper** — there is NO shared export; each component defines a tiny local `function cx(...parts: Array<string|false|null|undefined>): string { return parts.filter(Boolean).join(" "); }` (the house idiom — copy it into each new file).
- **Fixtures** — `src/print/fixtures/thumbs/{37,65,66,67,93}.png` (import `?url`). Reuse `SampleTableRow.stories.tsx`'s EXPOSURES/TAGS shapes.

### Design-guard contract (`npm run lint:design` MUST stay green)

`src/print/components/**` is **NOT** guard-exempt; `src/print/ui/**` **IS**. In `components/`:
- **ALLOWED:** token color classes (`bg-ink`, `text-paper`, `bg-plate`, `bg-accent`, `border-hair`/`border-hair-strong`, `text-ink-faint`), named type roles (`text-xs`/`text-sm`/`text-data`/`text-kicker`), `font-mono`/`font-semibold`/`font-bold`/`uppercase`/`tracking-*`, `opacity-*`, arbitrary SPACING/geometry (`bottom-7`, `min-w-[220px]`, `gap-1`, `px-2.5`), positioning (`fixed`, `left-1/2`, `-translate-x-1/2`, `z-50`), `rounded-md`/`rounded-full`, `shadow-lg`.
- **BANNED:** arbitrary-value APPEARANCE (`text-[…]`, `rounded-[10px]`, `bg-[…]`), raw `oklch(`/hex literals, side-stripe borders. If you need appearance a primitive doesn't expose → extend the primitive in `ui/` (Task 1 is the only such case here).

---

## Task 1: `Button` `ghostInverse` variant (refactor-on-contact)

**Why:** `CullBar` is a dark `bg-ink` floating pill. Its Restore/Clear actions are quiet text buttons that must read as *light muted → full paper on hover*. Every existing `Button` variant assumes the light paper surface (dark text), so on the dark bar they'd be invisible. Add one dark-surface ghost variant. The mockup distinguishes Restore (`oklch 0.80`) from Clear (`oklch 0.66`); we deliberately collapse both to one `ghostInverse` (consistent with the system's other collapses — radius→1 step, SymmetrySelector→SegmentedControl). The Reject button stays `variant="accent"`.

**Files:**
- Modify: `src/print/ui/Button.tsx`
- Modify: `src/print/ui/Button.stories.tsx`
- Test: `test/print-ui/Button.test.tsx`

- [ ] **Step 1 — failing test.** In `test/print-ui/Button.test.tsx`, add a case asserting the new variant renders and is tagged:
```tsx
it("renders the ghostInverse variant for dark surfaces", () => {
  render(<Button variant="ghostInverse">Restore</Button>);
  const btn = screen.getByRole("button", { name: "Restore" });
  expect(btn.getAttribute("data-variant")).toBe("ghostInverse");
});
```
Run: `npx vitest run test/print-ui/Button.test.tsx` → FAIL (TS: `"ghostInverse"` not assignable to `ButtonVariant`).

- [ ] **Step 2 — implement.** In `Button.tsx`: add `"ghostInverse"` to the `ButtonVariant` union, and a `variantClass` entry. `ui/` is guard-exempt, but prefer tokens + `opacity`:
```ts
ghostInverse:
  "bg-transparent text-paper opacity-70 hover:opacity-100 border border-transparent " +
  "focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent",
```
(Light paper text dimmed at rest, full on hover; transparent fill; sits on the dark bar.)

- [ ] **Step 3 — verify test passes.** Run: `npx vitest run test/print-ui/Button.test.tsx` → PASS.

- [ ] **Step 4 — story.** In `Button.stories.tsx`, add a `GhostInverse` story rendered on a dark wrapper so it's visible:
```tsx
export const GhostInverse: Story = {
  args: { variant: "ghostInverse", children: "Restore" },
  render: (args) => (
    <div className="bg-ink p-4 rounded-md inline-flex">
      <Button {...args} />
    </div>
  ),
};
```
(Match the file's existing meta/Story idiom — if it uses `satisfies`, keep that; if base `args` in meta, follow suit.)

- [ ] **Step 5 — gate.** Run: `npx vitest run test/print-ui/Button.test.tsx` (PASS) · `npm run lint:design` (clean) · `npx tsc --noEmit -p tsconfig.build.json` (clean).

- [ ] **Step 6 — commit** (named files only):
```bash
git add src/print/ui/Button.tsx src/print/ui/Button.stories.tsx test/print-ui/Button.test.tsx
git commit -m "feat(print/ui): Button ghostInverse variant — light-muted ghost for dark surfaces

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `CullBar` composite

**Mockup:** `docs/redesign-mockups/sample-table.html` — `.cullbar` CSS (lines 297–326) + markup (lines 570–574):
```html
<div class="cullbar"> <!-- .show toggles visible -->
  <span class="cb-count"><b id="cb-n">0</b> frames selected</span>
  <button class="cb-reject">Drop<span class="cb-k">X</span></button>
  <button class="cb-restore">Restore</button>
  <button class="cb-clear">Clear<span class="cb-k">Esc</span></button>
</div>
```
`.cullbar`: fixed, `left:50%; bottom:28px; translateX(-50%)`, flex row gap 4px, `bg-ink`/`text-paper`, radius 10px, big drop shadow, hidden (opacity 0 + `pointer-events:none` + `translateY(20px)`) until `.show`. `.cb-count` 12.5px/600 with a mono `<b>` number. `.cb-k` mono 9.5px dim — the keyboard hint inside reject/clear.

**Radius note:** the system has one radius step (`rounded-md` = 5px); use `rounded-md` (NOT `rounded-[10px]`, which the guard bans and which would reopen the collapsed-radius decision). `Toast` does the same. Note this deviation in the component doc-comment.

**Presentational** (Batch 7 contract): no `useState`; count + handlers are props.

**Files:**
- Create: `src/print/components/CullBar.tsx`
- Create: `src/print/components/CullBar.stories.tsx`
- Test: `test/print-components/CullBar.test.tsx`

### API
```tsx
export interface CullBarProps {
  /** Number of currently-selected frames. Shown as the mono count. */
  count: number;
  /** When true the bar is visible (slid up + interactive); false → hidden
   *  (faded, translated down, pointer-events:none) so it can animate in/out
   *  while staying mounted. Default false. */
  show?: boolean;
  onReject?: () => void;   // "Drop" (accent) — reject the selected frames
  onRestore?: () => void;  // "Restore" (ghostInverse) — un-reject
  onClear?: () => void;    // "Clear" (ghostInverse) — clear the selection
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}
```

### Structure
- Root `<div data-testid="cull-bar" data-show={show ? "true" : "false"}>` with classes:
  - base: `fixed left-1/2 bottom-7 -translate-x-1/2 z-50 flex items-center gap-1 bg-ink text-paper rounded-md shadow-lg pl-4 pr-2 py-[7px] transition-opacity`
  - visibility (via `cx`): `show ? "opacity-100" : "opacity-0 pointer-events-none"`
  - (`bottom-7` = 28px; `gap-1` = 4px. The translateY-in animation isn't load-bearing — `transition-opacity` + the opacity swap is enough; the page owns mount.)
- Count: `<span className="text-sm font-semibold mr-2.5"><b className="font-mono">{count}</b> frames selected</span>`. (`text-sm`≈14px nearest the 12.5px role; `mr-2.5`=10px.)
- Reject: `<Button variant="accent" onClick={onReject}>Drop<span className="font-mono text-xs opacity-60 ml-1">X</span></Button>`
- Restore: `<Button variant="ghostInverse" onClick={onRestore}>Restore</Button>`
- Clear: `<Button variant="ghostInverse" onClick={onClear}>Clear<span className="font-mono text-xs opacity-60 ml-1">Esc</span></Button>`
- Use `cx` for the root class composition (local helper). Forward `onClick` only when the handler is defined is NOT required — `Button` forwards `onClick={undefined}` harmlessly; pass directly.

- [ ] **Step 1 — failing test.** `test/print-components/CullBar.test.tsx`:
```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { CullBar } from "../../src/print/components/CullBar";

describe("CullBar", () => {
  it("shows the selected-frame count", () => {
    render(<CullBar count={3} show />);
    const bar = screen.getByTestId("cull-bar");
    expect(bar).toHaveAttribute("data-show", "true");
    expect(bar.textContent).toContain("3");
    expect(bar.textContent).toContain("frames selected");
  });

  it("is hidden (data-show=false) when show is omitted", () => {
    render(<CullBar count={0} />);
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
  });

  it("wires reject, restore, and clear", async () => {
    const onReject = vi.fn(), onRestore = vi.fn(), onClear = vi.fn();
    render(<CullBar count={2} show onReject={onReject} onRestore={onRestore} onClear={onClear} />);
    await userEvent.click(screen.getByRole("button", { name: /drop/i }));
    await userEvent.click(screen.getByRole("button", { name: /restore/i }));
    await userEvent.click(screen.getByRole("button", { name: /clear/i }));
    expect(onReject).toHaveBeenCalledOnce();
    expect(onRestore).toHaveBeenCalledOnce();
    expect(onClear).toHaveBeenCalledOnce();
  });

  it("uses the accent reject button and ghostInverse for restore/clear", () => {
    render(<CullBar count={1} show />);
    expect(screen.getByRole("button", { name: /drop/i }).getAttribute("data-variant")).toBe("accent");
    expect(screen.getByRole("button", { name: /restore/i }).getAttribute("data-variant")).toBe("ghostInverse");
    expect(screen.getByRole("button", { name: /clear/i }).getAttribute("data-variant")).toBe("ghostInverse");
  });
});
```
Run: `npx vitest run test/print-components/CullBar.test.tsx` → FAIL (module missing).

- [ ] **Step 2 — implement** `CullBar.tsx` per the API + structure above (local `cx`; doc-comment noting the presentational contract + the `rounded-md` radius deviation from the mockup's 10px).

- [ ] **Step 3 — verify test passes.** Run: `npx vitest run test/print-components/CullBar.test.tsx` → PASS.

- [ ] **Step 4 — stories.** `CullBar.stories.tsx`, `title: "components/CullBar"`. Render on a tall light wrapper (e.g. `<div className="bg-paper h-[200px] relative">`) so the fixed bar is visible in the frame. Stories: `Selected` (count 3, show), `Hidden` (count 0, show false), `Interactive` (a `useState` count harness with the three handlers nudging/clearing the count — mirrors LoupeSidePanel.Interactive).

- [ ] **Step 5 — gate.** `npx vitest run test/print-components/CullBar.test.tsx` (PASS) · `npm run lint:design` (clean) · `npx tsc --noEmit -p tsconfig.build.json` (clean).

- [ ] **Step 6 — commit:**
```bash
git add src/print/components/CullBar.tsx src/print/components/CullBar.stories.tsx test/print-components/CullBar.test.tsx
git commit -m "feat(print/components): CullBar — floating batch-reject action bar

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `SheetTable` composite

**Mockup:** `sample-table.html` — `.sheet` (lines 152–167) + `.sheet-head`/`.COLS` (lines 159–168, 487–498). The `.sheet` is a `Card`-elevated plate, `overflow-hidden`, containing a `.sheet-head` (paper-sunk header row, bottom hairline-strong, 5 uppercase column labels) then the rows region. Column labels: **Sample · Exposures · Kept · Tags · Status**.

**Children-slotting** (the `Gallery` pattern): the page maps samples → `SampleTableRow` (owning per-row handlers + selection); `SheetTable` lays them out under the aligned header. Column alignment holds because BOTH the header and every `SampleTableRow` use the exported `SAMPLE_TABLE_COLS`.

**Files:**
- Create: `src/print/components/SheetTable.tsx`
- Create: `src/print/components/SheetTable.stories.tsx`
- Test: `test/print-components/SheetTable.test.tsx`

### API
```tsx
export interface SheetTableProps {
  /** SampleTableRow elements, already mapped/sorted by the page. */
  children?: ReactNode;
  /** Rendered (instead of the rows region) when there are zero children. */
  empty?: ReactNode;
  /** PLACEMENT-ONLY, appended last. */
  className?: string;
}
```

### Structure
- Import `Card`, `Kicker` from `../ui`; `SAMPLE_TABLE_COLS` from `./SampleTableRow`; `React`/`ReactNode`. Local `cx`.
- Root: `<Card elevated className={cx("overflow-hidden", className)} data-testid="sheet-table">`. (Card forwards `data-testid` via `...rest`.)
- Header: `<div data-testid="sheet-head" className="bg-paper-sunk border-b border-hair-strong"><div className="grid" style={{ gridTemplateColumns: SAMPLE_TABLE_COLS }}>` then five `<Kicker tone="faint" className="px-4 py-2.5">Sample</Kicker>` … (`px-4`=16px, `py-2.5`=10px — matches the `.sheet-head` cell padding). Labels in order: Sample, Exposures, Kept, Tags, Status.
- Rows region: `<div data-testid="sheet-rows">{children}</div>`. If `React.Children.count(children) === 0 && empty != null`, render `<div data-testid="sheet-empty">{empty}</div>` inside the rows region instead of the children.

### Doc-comment must note
The grid template is the **shared exported `SAMPLE_TABLE_COLS`** (NOT a re-derived `grid-cols-[…]`) — that identity is what proves header/row column alignment. Children-slotting keeps `SheetTable` presentational; the page owns row data, per-row handlers, and selection.

- [ ] **Step 1 — failing test.** `test/print-components/SheetTable.test.tsx`:
```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { SheetTable } from "../../src/print/components/SheetTable";
import { SampleTableRow } from "../../src/print/components/SampleTableRow";

const row = (id: string) => (
  <SampleTableRow key={id} name={`Sample ${id}`} sampleId={id} exposures={[]} kept={1} total={1} tags={[]} />
);

describe("SheetTable", () => {
  it("renders the five aligned column headers", () => {
    render(<SheetTable>{[row("a")]}</SheetTable>);
    const head = screen.getByTestId("sheet-head");
    ["Sample", "Exposures", "Kept", "Tags", "Status"].forEach((label) =>
      expect(head.textContent).toContain(label),
    );
  });

  it("slots the provided rows", () => {
    render(<SheetTable>{[row("a"), row("b"), row("c")]}</SheetTable>);
    expect(screen.getAllByTestId("sample-table-row")).toHaveLength(3);
  });

  it("renders the empty node when there are no rows", () => {
    render(<SheetTable empty={<span>No samples</span>}>{[]}</SheetTable>);
    expect(screen.getByTestId("sheet-empty")).toBeInTheDocument();
    expect(screen.queryByTestId("sample-table-row")).toBeNull();
  });
});
```
Run: `npx vitest run test/print-components/SheetTable.test.tsx` → FAIL (module missing).

- [ ] **Step 2 — implement** `SheetTable.tsx` per the API + structure.

- [ ] **Step 3 — verify test passes.** Run: `npx vitest run test/print-components/SheetTable.test.tsx` → PASS.

- [ ] **Step 4 — stories.** `SheetTable.stories.tsx`, `title: "components/SheetTable"`. Build a small EXPOSURES/TAGS fixture (copy the shape from `SampleTableRow.stories.tsx`, importing the same `thumbs/*.png?url`). Render 3–4 `SampleTableRow` children (mix screened/unscreened, indexed/not-indexed, one with a dropped count). Stories: `Sheet` (full rows) and `Empty` (no children + an `empty` node). Wrap in `<div className="bg-paper p-6">`.

- [ ] **Step 5 — gate.** `npx vitest run test/print-components/SheetTable.test.tsx` (PASS) · `npm run lint:design` (clean) · `npx tsc --noEmit -p tsconfig.build.json` (clean).

- [ ] **Step 6 — commit:**
```bash
git add src/print/components/SheetTable.tsx src/print/components/SheetTable.stories.tsx test/print-components/SheetTable.test.tsx
git commit -m "feat(print/components): SheetTable — contact-sheet plate with aligned column header

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `ContactSheetAssembly` page-sim story + ledger update

A page-simulation story (NOT a component — no new `.tsx` component file) proving the full contact-sheet view composes: the `.head` (Contact-sheet `Kicker` + sub line + `6/9 samples screened` `ProgressBar`), the `SheetTable` with mapped `SampleTableRow` children, and the `CullBar` wired to a real selection-Set. Mirrors `LoupeAssembly.stories.tsx`.

**Files:**
- Create: `src/print/components/ContactSheetAssembly.stories.tsx`
- Modify: `docs/greenfield-component-ledger.md`

- [ ] **Step 1 — assembly story.** `ContactSheetAssembly.stories.tsx`, `title: "components/ContactSheetAssembly"`, `component: SheetTable` (for the Storybook docs anchor). A `ContactSheetView` function component with:
  - `useState<Set<number>>` for selected exposure ids across rows (the selection that drives the cullbar).
  - A `SAMPLES` fixture array (4–5 samples, real `thumbs/*.png?url` exposures, mixed screened/phase/dropped/tags).
  - Layout: `<div className="max-w-[1240px] mx-auto">` → a `.head` row (`<Kicker tone="accent">Contact sheet</Kicker>` + a sub `<p className="text-body text-ink-soft">` + a right-aligned `ProgressBar` block showing screened count) → `<SheetTable>` mapping `SAMPLES` to `SampleTableRow` (each row's `onSelectExposure` toggles the id in the Set; `selectedExposureId` reflects the most-recent or membership as suits the row API) → `<CullBar count={selected.size} show={selected.size > 0} onReject={…clear} onRestore={…} onClear={() => setSelected(new Set())} />`.
  - Keep it a faithful page sketch — the goal is fidelity verification, not production wiring. A single-select-per-row is fine if multi-frame selection is awkward with the row API; the cullbar count can track the Set size.
  - One default export; a `Page` story: `export const Page: Story = { render: () => <ContactSheetView /> };` (use `const meta: Meta<typeof SheetTable> =` — render-only).

- [ ] **Step 2 — gate the story.** Run: `npm run lint:design` (clean) · `npx tsc --noEmit -p tsconfig.build.json` (clean). (No unit test for the assembly — it's a visual sim, like `LoupeAssembly`.)

- [ ] **Step 3 — ledger.** In `docs/greenfield-component-ledger.md`: flip the `SheetTable` (line ~129) and `CullBar` (line ~170) rows to ✅ with file paths; add a **Batch 11** row to the *Out-of-scope & deferred registry* table summarizing the slice (composites built, the `Button.ghostInverse` refactor-on-contact, the `rounded-md` radius deviation, the children-slotting + presentational decisions, the Layer-4 deferral of `.sheet-shell`/`.kb-legend`/loupe-view/selection-state); update the L3 coverage tally (`SheetTable` + `CullBar` move from ⬜ to ✅; remaining ⬜ = `ScopePlate`+scoping `SampleRow`, `BuilderRail`); strike contact-sheet from the next-up bullets.

- [ ] **Step 4 — commit:**
```bash
git add src/print/components/ContactSheetAssembly.stories.tsx docs/greenfield-component-ledger.md
git commit -m "feat(print/components): contact-sheet assembly story + Batch 11 ledger

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Batch verification (after all tasks)

- `npm run build-storybook` → exit 0.
- Serve `storybook-static` on a port; Playwright-MCP visual check of `components/SheetTable` (header aligned over rows), `components/CullBar` (dark bar, accent Drop + muted Restore/Clear), and `components/ContactSheetAssembly` (full view; cullbar appears on selection) vs `sample-table.html`.
- Final: `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json` clean across the whole tree.
- Clean up: kill the server, remove screenshots + `.playwright-mcp` dirs.

## Out of scope (this slice)

`.sheet-shell` plate-shell page component, `.kb-legend`, the loupe view (built in Batch 10 — assembly only here), real selection/keyboard wiring, the `Beamtime` facet filter, the `view-seg` Contact-sheet/Loupe toggle. All Layer-4 page assembly.
