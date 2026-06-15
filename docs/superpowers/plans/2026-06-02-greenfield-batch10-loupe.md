# Greenfield Batch 10 — Loupe slice (`BigFrame` + `LoupeSidePanel`)

> Execute with superpowers:subagent-driven-development (fresh implementer per task + frontend-reviewer + fix loop).
> Commit trailer (every commit, exact last line): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Commit ONLY the named files per task (never `git add -A`/`.`). Never stage `src/bones/*`, `docs/superpowers/plans/*`.
> Typecheck: `npx tsc --noEmit -p tsconfig.build.json`. Work from the worktree frontend dir.

## Context

The greenfield "The Print" rebuild (branch `worktree-greenfield-ui-rebuild`, NOT merged). Layer-3 frontier. This slice builds the **Loupe** surface's two composites, derived from `../../../docs/redesign-mockups/sample-table.html` (the `.loupe-shell` region, CSS lines 325–442, markup lines 510–564 — REPO ROOT relative to frontend dir).

**All leaf deps already exist** — this is pure composition, no `ui/` refactor-on-contact:
- `BigFrame` ← `DetectorImage` (`src/print/detector/`) + `RejectOverlay` (`src/print/ui/`)
- `LoupeSidePanel` ← `MetaList` + `Kicker` + `KbLegend` (`src/print/ui/`) + `Verdict` + `RepresentativeBox` (`src/print/components/`) + `TagList` (`src/print/ui/`)

**Verified leaf APIs (do not re-derive — these are checked against live source):**
- `DetectorImage`: `{ src: string | null; size: "thumb" | "full"; lutVariant?: DetectorLutVariant; className? }`. Placeholder `data-testid="detector-image-placeholder"` when `src===null`. Forwards `className` to its wrapper `<div>`. Mirror `DetectorPanel`'s usage: `size="full"`, `className="h-full w-full"`, conditional-spread `lutVariant`.
- `DetectorLutVariant = "neutral" | "warm"` — `import type { DetectorLutVariant } from "../detector/detectorLut"`.
- `RejectOverlay`: `{ className? }`. `data-testid="reject-overlay"`, `aria-hidden`. SVG terracotta ✕, scales to fill a `relative` parent (give it `absolute inset-0`).
- `MetaList`: `{ entries: MetaEntry[]; className? }`. `MetaEntry = { key: string; value: ReactNode }`. `data-testid="meta-list"`. `import { MetaList, type MetaEntry } from "../ui"` (barrel).
- `Verdict`: `{ dropped: boolean; hint?: ReactNode; onToggle?: () => void; className? }`. Renders a `Dot` + bold "Dropped"/"Kept" + a `Button` labelled "Restore"/"Drop". `import { Verdict } from "./Verdict"`.
- `RepresentativeBox`: `{ isRepresentative: boolean; onSetRepresentative?: () => void; children?: ReactNode; className? }`. `Button` labelled "Set as representative". `import { RepresentativeBox } from "./RepresentativeBox"`.
- `TagList`: `{ tags: Tag[]; onAdd?; onRemove?; editable?; size?; maxVisible?; className? }`. `Tag = { key: string; value?: string }` from `../ui/tag`. `data-testid="tag-list"`. `import { TagList } from "../ui"`.
- `KbLegend`: `{ shortcuts: Shortcut[]; className? }`. `Shortcut = { keyLabel: string; description: string }`. `data-testid="kb-legend"`. `import { KbLegend, type Shortcut } from "../ui"`.
- `Kicker`: `{ tone?: "accent" | "faint"; as?: "div"|"span"|"h2"|"h3"; ...HTMLAttributes }`. `.side-h` == `Kicker tone="faint"` (uppercase 700 tracked ink-faint), use `as="h3"`. `import { Kicker } from "../ui"`.
- Tokens available (guard-clean, generated from `--color-*` in styles.css `@theme`): `bg-frame-edge`, `text-frame-tag` (the light caption tint over the dark frame), `bg-accent`, `text-plate` (near-white, for on-accent text).

**Design-guard reminder:** `src/print/components/**` is NOT exempt. Compose primitives + layout/token utilities only. ALLOWED: named type roles (`text-display`), token color classes (`bg-frame-edge`/`text-frame-tag`/`bg-accent`/`text-plate`/`text-ink-faint`), arbitrary SPACING/geometry (`px-2`, `py-[3px]`, `max-w-[500px]`), `opacity-40`, `uppercase`/`font-bold`/`tracking-wide`. BANNED: `text-[…]`, `rounded-[…]`, raw `oklch(`/hex, side-stripe borders. `npm run lint:design` must stay green.

**Scope boundary (do NOT build):** the `.loupe-plate` shell (bordered plate + `loupe-back` button + serif `loupe-head` + 2-col `loupe-body` grid) and the `loupe-strip` filmstrip are **Layer-4 page assembly** — deferred, exactly as Batch 7/9 deferred their pages. A combined "loupe assembly" *story* (simulating the page with a `ThumbnailGallery` filmstrip + 2-col grid, inline in the story) gives visual verification without building the page component.

---

## Task 1: `BigFrame`

The large square detector frame: the loupe's hero. Composes `DetectorImage` (filling a square frame) + a `RejectOverlay` (the big ✕) + a "Dropped" accent pill + a mono caption — the last three shown/styled per `rejected`.

**Files:**
- Create: `src/print/components/BigFrame.tsx`
- Create: `src/print/components/BigFrame.stories.tsx`
- Test: `test/print-components/BigFrame.test.tsx`

### Spec (from mockup `.big-frame`, CSS 354–374, markup 519–522)

- Root `<div data-testid="big-frame" data-rejected={rejected ? "true" : undefined}>`: `relative`, `aspect-square`, `rounded`, `overflow-hidden`, `bg-frame-edge border border-frame-edge`, `max-w-[500px]`, `mx-auto`. Forward placement `className`.
- `DetectorImage` fills it: `size="full"`, `className={cx("h-full w-full", rejected && "opacity-40")}`, conditional-spread `lutVariant`.
- When `rejected`:
  - "Dropped" pill — `data-role="dropped-tag"`, absolute top-left (`absolute left-3 top-3`), `bg-accent text-plate`, `uppercase font-bold tracking-wide`, `text-xs`, `rounded-sm`, `px-2 py-[3px]`. Text: `Dropped`.
  - `<RejectOverlay className="absolute inset-0" />` (the big ✕).
- Caption (always when `caption` given) — `data-role="frame-caption"`, absolute bottom-left (`absolute bottom-2.5 left-3`), `font-mono text-frame-tag text-sm tracking-wide`. Renders `{caption}`.

### Props

```tsx
import type { ReactNode } from "react";
import type { DetectorLutVariant } from "../detector/detectorLut";

export interface BigFrameProps {
  /** Detector image URL; null → DetectorImage frame-window placeholder. */
  src: string | null;
  /** Mono caption set over the dark frame (the `.frame-tag`), e.g. "frame 65 · 0.40 s". */
  caption?: ReactNode;
  /** true → dims the image, shows the "Dropped" pill + the grease-pencil ✕. */
  rejected?: boolean;
  /** Forwarded to DetectorImage (colormap). */
  lutVariant?: DetectorLutVariant;
  /** PLACEMENT ONLY. */
  className?: string;
}
```

- [ ] **Step 1 — failing test** `test/print-components/BigFrame.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import { BigFrame } from "../../src/print/components/BigFrame";

describe("BigFrame", () => {
  it("renders the frame with a detector-image placeholder when src is null", () => {
    render(<BigFrame src={null} />);
    expect(screen.getByTestId("big-frame")).toBeInTheDocument();
    expect(screen.getByTestId("detector-image-placeholder")).toBeInTheDocument();
  });

  it("renders the mono caption when given", () => {
    render(<BigFrame src={null} caption="frame 65 · 0.40 s" />);
    const cap = screen.getByText("frame 65 · 0.40 s");
    expect(cap).toHaveAttribute("data-role", "frame-caption");
  });

  it("kept (default): no reject overlay, no dropped tag, data-rejected absent", () => {
    render(<BigFrame src={null} caption="f1" />);
    expect(screen.queryByTestId("reject-overlay")).not.toBeInTheDocument();
    expect(screen.queryByText("Dropped")).not.toBeInTheDocument();
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected", "true");
  });

  it("rejected: shows the dropped tag and the reject overlay, flags data-rejected", () => {
    render(<BigFrame src={null} caption="f1" rejected />);
    expect(screen.getByTestId("reject-overlay")).toBeInTheDocument();
    const tag = screen.getByText("Dropped");
    expect(tag).toHaveAttribute("data-role", "dropped-tag");
    expect(screen.getByTestId("big-frame")).toHaveAttribute("data-rejected", "true");
  });
});
```

- [ ] **Step 2** — run `npm test -- BigFrame` → FAIL (module not found).
- [ ] **Step 3** — implement `src/print/components/BigFrame.tsx` per spec (local `cx` helper as in sibling files; `DetectorImage` from `../detector`, `RejectOverlay` from `../ui`).
- [ ] **Step 4** — write `BigFrame.stories.tsx`: `title: "components/BigFrame"`, `component: BigFrame`. Stories: `Kept` (real fixture `import thumb37 from "../fixtures/thumbs/37.png?url"`, `src={thumb37}`, `caption="frame 37 · 0.40 s"`), `Dropped` (same src, `rejected`, caption), `Empty` (`src={null}`, caption). Wrap each in a sizing container (e.g. `<div className="bg-paper p-6 w-[520px]">`).
- [ ] **Step 5** — `npm test -- BigFrame` (green) + `npm run lint:design` (clean) + `npx tsc --noEmit -p tsconfig.build.json` (clean).
- [ ] **Step 6** — commit ONLY `src/print/components/BigFrame.tsx src/print/components/BigFrame.stories.tsx test/print-components/BigFrame.test.tsx`.

---

## Task 2: `LoupeSidePanel`

The loupe aside (`.loupe-side`): a vertical stack (gap 18px) of five blocks — "This exposure" (`MetaList`), the `Verdict`, the `RepresentativeBox`, "Sample tags" (`TagList`), and "Keys" (`KbLegend`). **Presentational** (Batch 7 contract): holds NO interaction `useState`; every datum is a prop with a page-owned handler.

**Files:**
- Create: `src/print/components/LoupeSidePanel.tsx`
- Create: `src/print/components/LoupeSidePanel.stories.tsx`
- Test: `test/print-components/LoupeSidePanel.test.tsx`

### Spec (from mockup `.loupe-side`, CSS 379–442, markup 525–561)

- Root `<aside data-testid="loupe-side-panel">`: `flex flex-col gap-[18px]`. Forward placement `className`.
- Section header helper = `<Kicker tone="faint" as="h3" className="mb-2">{label}</Kicker>` (the `.side-h`).
- Block 1 — `<div>` (section): header "This exposure" + `<MetaList entries={meta} />`.
- Block 2 — `<Verdict dropped={dropped} {...(onToggleDrop ? { onToggle: onToggleDrop } : {})} />` (no header).
- Block 3 — `<RepresentativeBox isRepresentative={isRepresentative} {...(onSetRepresentative ? { onSetRepresentative } : {})} />` (no header).
- Block 4 — section: header "Sample tags" + `<TagList tags={tags} editable {...(onAddTag ? { onAdd: onAddTag } : {})} {...(onRemoveTag ? { onRemove: onRemoveTag } : {})} />`.
- Block 5 — section: header "Keys" + `<KbLegend shortcuts={shortcuts} className="flex-col gap-1" />` (the loupe lists keys vertically; `KbLegend` is a wrapping flex row — `flex-col` placement override stacks them).

`exactOptionalPropertyTypes: true` is on — forward optional callbacks to children via conditional spread (children's matching props are not `| undefined`).

### Props

```tsx
import type { MetaEntry, Shortcut, Tag } from "../ui"; // re-exported types; verify barrel — else import Tag from "../ui/tag"

export interface LoupeSidePanelProps {
  /** "This exposure" metadata rows. */
  meta: MetaEntry[];
  /** Verdict state + toggle. */
  dropped: boolean;
  onToggleDrop?: () => void;
  /** Representative state + setter. */
  isRepresentative: boolean;
  onSetRepresentative?: () => void;
  /** Sample tags + edit handlers. */
  tags: Tag[];
  onAddTag?: (t: Tag) => void;
  onRemoveTag?: (t: Tag) => void;
  /** Key legend — defaults to the canonical loupe keys. */
  shortcuts?: Shortcut[];
  /** PLACEMENT ONLY. */
  className?: string;
}

export const LOUPE_KEYS: Shortcut[] = [
  { keyLabel: "← →", description: "flip frames" },
  { keyLabel: "X", description: "drop / restore" },
  { keyLabel: "R", description: "set representative" },
  { keyLabel: "Esc", description: "back to the sheet" },
];
```

(In the signature use `shortcuts = LOUPE_KEYS` as the default.)

> NOTE for implementer: confirm `MetaEntry`, `Shortcut`, `Tag` are re-exported from `../ui` (barrel `src/print/ui/index.ts`). If any is not, import it from its source (`../ui/tag` for `Tag`, `../ui/MetaList` for `MetaEntry`, `../ui/KbLegend` for `Shortcut`). Verify before writing the import — do not assume.

- [ ] **Step 1 — failing test** `test/print-components/LoupeSidePanel.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { LoupeSidePanel } from "../../src/print/components/LoupeSidePanel";

const meta = [
  { key: "frame", value: "65" },
  { key: "exposure", value: "0.40 s" },
];
const tags = [{ key: "LL37" }, { key: "temp", value: "37C" }];

function setup(over = {}) {
  const props = {
    meta,
    dropped: false,
    isRepresentative: false,
    tags,
    ...over,
  };
  render(<LoupeSidePanel {...props} />);
  return props;
}

describe("LoupeSidePanel", () => {
  it("renders the four section headers + the composed leaves", () => {
    setup();
    expect(screen.getByTestId("loupe-side-panel")).toBeInTheDocument();
    expect(screen.getByText("This exposure")).toBeInTheDocument();
    expect(screen.getByText("Sample tags")).toBeInTheDocument();
    expect(screen.getByText("Keys")).toBeInTheDocument();
    expect(screen.getByTestId("meta-list")).toBeInTheDocument();
    expect(screen.getByTestId("tag-list")).toBeInTheDocument();
    expect(screen.getByTestId("kb-legend")).toBeInTheDocument();
  });

  it("reflects kept verdict and wires the drop toggle", async () => {
    const onToggleDrop = vi.fn();
    setup({ onToggleDrop });
    // Verdict renders a Button labelled "Drop" when kept.
    await userEvent.click(screen.getByRole("button", { name: /drop/i }));
    expect(onToggleDrop).toHaveBeenCalledOnce();
  });

  it("wires set-representative", async () => {
    const onSetRepresentative = vi.fn();
    setup({ onSetRepresentative });
    await userEvent.click(screen.getByRole("button", { name: /set as representative/i }));
    expect(onSetRepresentative).toHaveBeenCalledOnce();
  });

  it("shows the default loupe keys", () => {
    setup();
    expect(screen.getByText("flip frames")).toBeInTheDocument();
    expect(screen.getByText("set representative")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2** — run `npm test -- LoupeSidePanel` → FAIL.
- [ ] **Step 3** — implement `src/print/components/LoupeSidePanel.tsx` per spec. Verify the type-import barrel claim first (see NOTE).
- [ ] **Step 4** — write `LoupeSidePanel.stories.tsx`: `title: "components/LoupeSidePanel"`. Stories: `Kept` (static, noop handlers, `isRepresentative: false`), `DroppedRepresentative` (`dropped`, `isRepresentative`), and an `Interactive` story with a `useState` wrapper (page simulation: local `dropped`/`isRep` state toggled by the handlers). Wrap in `<div className="bg-paper p-5 w-[300px]">`.
- [ ] **Step 5** — `npm test -- LoupeSidePanel` + `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json` (all clean).
- [ ] **Step 6** — commit ONLY the three named files.

---

## Task 3: Combined assembly story + ledger update

- [ ] **Step 1** — Add a `LoupeAssembly` story (in `BigFrame.stories.tsx`, or a new `src/print/components/LoupeSidePanel.stories.tsx` export) that simulates the Layer-4 loupe page: a 2-col grid (`grid grid-cols-[1fr_286px] gap-7`) with `BigFrame` + a `ThumbnailGallery align="center" size="lg"` filmstrip on the left and `LoupeSidePanel` on the right, using inline mockup data + a `useState` page-sim harness (selected frame, dropped, rep). This is for visual fidelity verification only — NOT a new component. Use real fixture thumbs.
- [ ] **Step 2** — `npm run build-storybook` (exit 0). Visually verify `BigFrame` (kept + dropped), `LoupeSidePanel`, and the assembly against `sample-table.html` loupe region via Playwright MCP against `storybook-static`.
- [ ] **Step 3** — Update `docs/greenfield-component-ledger.md`: flip `BigFrame` + `LoupeSidePanel` ⬜→✅ with paths; add a Batch 10 decision-registry row (loupe slice = pure composition, no new primitive; plate-shell/filmstrip/back/head deferred to Layer-4 page); update L3 coverage counts; strike loupe from the next-up bullets.
- [ ] **Step 4** — Update memory (`project_greenfield_composite_layer.md` + `MEMORY.md` index line) with the Batch 10 paragraph + revised NEXT (remaining: Contact sheet, Series scoping, Series builder).
- [ ] **Step 5** — commit ONLY the story file(s) touched + `docs/greenfield-component-ledger.md` (memory files are outside the repo, written separately — never staged).

## Verification (whole batch)

`npm test -- print-components` green · `npm run lint:design` clean · `npx tsc --noEmit -p tsconfig.build.json` clean · `npm run build-storybook` exit 0 · frontend-reviewer clean · all surfaces visually verified vs `sample-table.html`.
