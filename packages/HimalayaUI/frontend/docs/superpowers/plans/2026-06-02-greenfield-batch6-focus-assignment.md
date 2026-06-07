# Greenfield Batch 6 — Focus assignment slice (`AssignmentCart`)

> Execute with superpowers:subagent-driven-development (fresh implementer + review).
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Branch `worktree-greenfield-ui-rebuild` — NOT merged, NOT pushed. Do NOT offer merge/finish/PR.
> Commit ONLY the named files per task (never `git add -A`). Never touch/stage `src/bones/registry.ts`,
> `src/bones/contact-sheet.bones.json`, or any pre-existing untracked plan doc.
> Typecheck with `npx tsc --noEmit -p tsconfig.build.json` (root tsc has unrelated errors).

## Context

The greenfield "The Print" rebuild's Layers 0–2 are complete. The build frontier is the
**Layer-3 tier-2 panels** (`src/print/components/`, design-guard NON-exempt → placement-only:
named-role type classes + token utilities + layout only; NO `text-[…]`/`rounded-[…]`/raw
`oklch(`/hex/side-stripes/inline style colors). The user chose the **Focus assignment** slice.

Verified against live source this session:
- `PhaseBlock` (`src/print/components/PhaseBlock.tsx`) EXISTS: `{ phase; score:0..1; meta:ReactNode;
  onRemove?; className? }` — serif name (`text-headline`) + score + dismiss `IconButton` + `meta`
  line (`text-data text-ink-soft`) + `ScoreBar size="bar"`.
- `CandidateRow` + `CandidateList` (`src/print/components/CandidateRow.tsx`) EXIST. `CandidateList`
  is the thin `flex flex-col gap-2` wrapper. **No new work needed on the candidate side.**
- So the one genuinely-new panel is **`AssignmentCart`** — the `.phasecall` plate from
  `docs/redesign-mockups/focus-workspace.html`.

Mockup source of truth — `docs/redesign-mockups/focus-workspace.html`:
- `.rail` aside → `.rail-sec` "Phase call" → `.phasecall` plate (`bg plate`, `border hair`, `radius`,
  `overflow hidden`) containing EITHER:
  - `.pc-empty` (`padding 15px; font-size 12px; color ink-faint`) — copy: "No phase assigned —
    every peak is unindexed. Check a candidate below." (mockup JS lines 890-891), OR
  - an optional `.pc-tag` (`padding 9px 15px 0`; 10px uppercase 0.08em bold ink-faint) reading
    "Coexistence · N phases" — shown ONLY when >1 phase (mockup JS line 893-894) — followed by one
    `.pc-block` per phase, separated by `.pc-block + .pc-block { border-top: 1px hair }`
    (`padding: 11px 15px 13px`).
- Each `.pc-block` = `PhaseBlock` (name/score/meta/bar). The mockup block ALSO carries a
  `.pcb-series` line ("series <name>", mono ink-soft) below the bar — PhaseBlock currently lacks it.
  Task A adds it as an OPTIONAL content-only prop (guard-clean: text only).

Verified token/utility facts (cite these; they resolve):
- `.pc-tag` → `text-kicker text-kicker-faint` (styles.css:270-273; `text-kicker` is uppercase/700/
  tracked, color comes from the `-faint`/`-accent` modifier).
- `.pc-empty` copy → `text-caption text-ink-faint` (styles.css:246,268 — the named placeholder role).
- plate chrome → `bg-plate border border-hair rounded overflow-hidden` (all `@theme` tokens).
- full-bleed dividers via the Gallery-proven `React.Children.map` + positional `border-t border-hair`
  on a full-width child wrapper (NOT `divide-y` on a padded wrapper — that would inset the hairline).
- `cx(...)` is a per-file local helper (copy the one in `Gallery.tsx`/`SeriesCard.tsx`).
- horizontal block padding `px-4` (≈15px) on each block wrapper + tag + empty; PhaseBlock's own
  `py-3` gives the vertical.

## Task A — `PhaseBlock` gains an optional `series` line (refactor-on-contact)

**Files:**
- Modify: `src/print/components/PhaseBlock.tsx`
- Modify test: `test/print-components/PhaseBlock.test.tsx`

- [ ] **Step 1 — failing test.** Add to `PhaseBlock.test.tsx`:

```tsx
it("renders an optional series line when `series` is provided", () => {
  render(<PhaseBlock phase="Im3m" score={0.8} meta="a = 14.2 nm · 5 reflections" series="QII alkane" />);
  expect(screen.getByText(/QII alkane/)).toBeInTheDocument();
  expect(screen.getByText(/series/i)).toBeInTheDocument();
});

it("omits the series line when `series` is absent", () => {
  render(<PhaseBlock phase="Im3m" score={0.8} meta="a = 14.2 nm" />);
  expect(screen.queryByText(/^series/i)).not.toBeInTheDocument();
});
```

- [ ] **Step 2 — run, expect fail** (`series` not a prop yet). `npm test -- PhaseBlock`.

- [ ] **Step 3 — implement.** Add `series?: ReactNode;` to `PhaseBlockProps`. After the `ScoreBar`,
  when `series != null` render a line mirroring `.pcb-series` (label faint, value ink-soft, mono):

```tsx
{series != null && (
  <div className="text-caption text-ink-faint mt-1.5">
    series <b className="text-ink-soft font-semibold font-mono">{series}</b>
  </div>
)}
```

  (`font-mono` is a layout/family utility, not appearance-color — guard-clean. `text-caption`/
  `text-ink-*` are named roles/tokens.) Keep everything else unchanged.

- [ ] **Step 4 — run, expect pass.** `npm test -- PhaseBlock`.
- [ ] **Step 5 — guard + tsc.** `npm run lint:design` clean; `npx tsc --noEmit -p tsconfig.build.json` clean.
- [ ] **Step 6 — commit** ONLY `src/print/components/PhaseBlock.tsx test/print-components/PhaseBlock.test.tsx`:

```
feat(print): PhaseBlock optional series line (.pcb-series from focus-workspace mockup)
```

## Task B — `AssignmentCart` panel

**Files:**
- Create: `src/print/components/AssignmentCart.tsx`
- Create: `src/print/components/AssignmentCart.stories.tsx`
- Create: `test/print-components/AssignmentCart.test.tsx`

- [ ] **Step 1 — failing test.** `test/print-components/AssignmentCart.test.tsx`:

```tsx
import { render, screen } from "@testing-library/react";
import { AssignmentCart } from "../../src/print/components/AssignmentCart";
import { PhaseBlock } from "../../src/print/components/PhaseBlock";

const block = (phase: string) => (
  <PhaseBlock key={phase} phase={phase} score={0.8} meta="a = 14.2 nm · 5 reflections" />
);

describe("AssignmentCart", () => {
  it("shows the empty state when there are no phase blocks", () => {
    render(<AssignmentCart />);
    expect(screen.getByTestId("assignment-empty")).toBeInTheDocument();
    expect(screen.getByText(/no phase assigned/i)).toBeInTheDocument();
    expect(screen.queryByTestId("coexistence-tag")).not.toBeInTheDocument();
  });

  it("lets the caller override the empty copy", () => {
    render(<AssignmentCart empty="Nothing here yet" />);
    expect(screen.getByText("Nothing here yet")).toBeInTheDocument();
  });

  it("renders one block and NO coexistence tag for a single phase", () => {
    render(<AssignmentCart>{block("Im3m")}</AssignmentCart>);
    expect(screen.getAllByTestId("phase-block")).toHaveLength(1);
    expect(screen.queryByTestId("coexistence-tag")).not.toBeInTheDocument();
    expect(screen.getByTestId("assignment-cart").dataset.phaseCount).toBe("1");
  });

  it("shows a 'Coexistence · N phases' tag when more than one phase is present", () => {
    render(<AssignmentCart>{[block("Ia3d"), block("Im3m")]}</AssignmentCart>);
    expect(screen.getAllByTestId("phase-block")).toHaveLength(2);
    const tag = screen.getByTestId("coexistence-tag");
    expect(tag).toHaveTextContent(/coexistence/i);
    expect(tag).toHaveTextContent(/2 phases/);
    expect(screen.getByTestId("assignment-cart").dataset.phaseCount).toBe("2");
  });

  it("forwards a placement-only className", () => {
    render(<AssignmentCart className="mt-4">{block("Im3m")}</AssignmentCart>);
    expect(screen.getByTestId("assignment-cart").className).toContain("mt-4");
  });
});
```

- [ ] **Step 2 — run, expect fail** (module missing). `npm test -- AssignmentCart`.

- [ ] **Step 3 — implement** `src/print/components/AssignmentCart.tsx`:

```tsx
import { Children, type ReactNode } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface AssignmentCartProps {
  /** PhaseBlock children — one per assigned phase. */
  children?: ReactNode;
  /** Override the empty-state copy. */
  empty?: ReactNode;
  /** PLACEMENT-ONLY. Appended last. */
  className?: string;
}

const DEFAULT_EMPTY =
  "No phase assigned — every peak is unindexed. Check a candidate below.";

export function AssignmentCart({
  children,
  empty,
  className,
}: AssignmentCartProps): JSX.Element {
  const count = Children.count(children);
  return (
    <div
      data-testid="assignment-cart"
      data-phase-count={count}
      className={cx(
        "bg-plate border border-hair rounded overflow-hidden",
        className,
      )}
    >
      {count === 0 ? (
        <div
          data-testid="assignment-empty"
          className="text-caption text-ink-faint p-4"
        >
          {empty ?? DEFAULT_EMPTY}
        </div>
      ) : (
        <>
          {count > 1 && (
            <div
              data-testid="coexistence-tag"
              className="text-kicker text-kicker-faint px-4 pt-2.5"
            >
              Coexistence · {count} phases
            </div>
          )}
          {Children.map(children, (child, i) => (
            <div className={cx("px-4", i > 0 && "border-t border-hair")}>
              {child}
            </div>
          ))}
        </>
      )}
    </div>
  );
}
```

- [ ] **Step 4 — run, expect pass.** `npm test -- AssignmentCart`.

- [ ] **Step 5 — stories** `src/print/components/AssignmentCart.stories.tsx`. Idiom: `import type
  { Meta, StoryObj } from "@storybook/react-vite";` + loose `const meta: Meta<typeof AssignmentCart>
  = {...}` (NOT `satisfies`); `@storybook/test` NOT installed. Stories: `Empty` (no children),
  `SinglePhase` (one PhaseBlock, score 0.84, meta "a = 14.2 nm · 5 reflections", series "QII alkane"),
  `Coexistence` (two PhaseBlocks: Ia3d 0.71 + Im3m 0.66 with `onRemove={() => {}}` so the dismiss
  affordance shows). Wrap each in a `~300px`-wide `bg-paper-sunk p-4` div to mimic the rail.

- [ ] **Step 6 — guard + tsc + full suite.** `npm run lint:design` clean; `npx tsc --noEmit -p
  tsconfig.build.json` clean; `npm test -- AssignmentCart PhaseBlock` green.

- [ ] **Step 7 — commit** ONLY the three new files:

```
feat(print): AssignmentCart panel (phase-call plate + coexistence tag + empty state)
```

## Verification (whole batch, run by controller)

- `npx tsc --noEmit -p tsconfig.build.json` → exit 0
- `npm run lint:design` → clean
- `npm test` → full suite green (expect +~7 tests)
- `npm run build-storybook` → exit 0
- Visually confirm `components/AssignmentCart` Empty / SinglePhase / Coexistence stories against the
  `.phasecall` region of `focus-workspace.html`.

## Task C — ledger + memory (controller, after B is reviewed-green)

- Flip the `AssignmentCart` / `CandidateList` Layer-3 rows to ✅ (with paths) in
  `docs/greenfield-component-ledger.md`; bump the coverage summary (Layer-3 panels 3 → 5 ✅);
  add a decisions-registry entry (AssignmentCart owns count→tag + positional-divider + empty;
  CandidateList was already built in Batch 1's CandidateRow file). Note the optional `series` prop
  added to PhaseBlock.
- Update `project_greenfield_composite_layer.md` + `MEMORY.md` index line (new NEXT = Series reading
  or Focus plates).
```