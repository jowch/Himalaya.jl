# Greenfield Batch 7 — Focus plates (TracePlate · DetectorPanel · CombsPanel) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Commit ONLY the named files per task (never `git add -A`). Do NOT touch `src/bones/registry.ts`, `src/bones/contact-sheet.bones.json`, or untracked plan docs.

**Goal:** Build the Focus workspace's three hero panels as Layer-3 composites — `TracePlate` (hero trace), `DetectorPanel` (detector frame + exposure strip), `CombsPanel` (comb/residual toggle) — wiring the already-built renderers (`TracePlot`, `DetectorImage`, `CombChart`/`ResidualChart`) into placement-only panel chrome derived from `docs/redesign-mockups/2026-05-29-focus-plot.html`.

**Architecture:** Each panel composes existing primitives + renderers + a new tiny shared `PanelHeader`. Closed-look / placement-only: chrome via the `Card` primitive (`elevated` for the hero plate, flat for the two lower panels), labels via `Kicker tone="faint"`, no inline appearance literals. Cross-child / interaction state (scale mode, exposure selection, comb view, hovered-q) is **lifted to the consumer** — panels are presentational and take state + callbacks (the Batch 6 "container owns cross-child state" rule, applied one level up: the *page* owns these, the panel renders them).

**Tech Stack:** React + TypeScript (strict), Tailwind v4 token utilities, Storybook (`@storybook/react-vite`), Vitest + RTL.

---

## Verified facts (checked against live source 2026-06-02 — do NOT re-derive, but DO re-read the cited file before editing)

- **`Card`** (`src/print/ui/Card.tsx`): `as`, `elevated` (→ `.card` lifted shadow; the hero plate), `border` ("hair"|"strong"), `draft`, `padding` ("sm"|"md"|"lg" → p-3/p-4/p-7), `selected`, `className`. Flat default = `rounded-md bg-plate border border-hair`, NO shadow (the two lower panels). Pads asymmetrically? Pass NO `padding` and use a placement padding utility.
- **`PlateHeader`** (`src/print/components/PlateHeader.tsx`): `kicker?`, `title`, `subtitle?`, `as` ("h1"|"h2"|"h3", default "h2"), `children` (right tools slot), `className`. Renders accent Kicker + `text-display` serif title + `text-data text-ink-faint` mono subtitle + right slot. This IS the hero `.plate-head`.
- **`ToolBar`** (`src/print/components/ToolBar.tsx`): `children`, `className`. `role="toolbar"`, `flex items-center gap-2 flex-shrink-0`.
- **`Button`** (`src/print/ui/Button.tsx`): `variant` ("solid"|"accent"|"ghost"|"danger"|"outline", default "ghost") + **`armed`** boolean (terracotta fill, `aria-pressed`, `data-armed="true"`). Auto-fit = `<Button variant="ghost">`; "+ Peak" = `<Button armed={addPeakArmed}>`.
- **`SegmentedControl<T>`** (`src/print/ui/SegmentedControl.tsx`): `options: {value,label,...}[]`, `value`, `onChange`, `variant` ("bordered"|"plain"), `size` ("xs"|"sm"|"md"; **xs is the "focus-plot mini comb switch"**), `aria-label` (REQUIRED), `testId`, `className`. Active segment = `bg-ink text-paper`.
- **`Kicker`** (`src/print/ui/Kicker.tsx`): `tone` ("accent"|"faint", default "faint"), `as` ("div"|"span"|"h2"|"h3"). `tone="faint"` = the `.panel-h` label (uppercase 700 tracked `text-ink-faint`).
- **`TracePlot`** (`src/print/plot/TracePlot.tsx`): `trace: TraceModel` (**required**, `{ trace: Trace; peaks: PlotPeak[]; phase: string|null }` — peaks & phase NON-optional), `xType?`/`yType?: ScaleType`, `axes?`, `xLabel?`/`yLabel?`, `interaction?: TracePlotInteraction|false`, `height?`, `width?`, `paperColor?`, `className?`, `"data-testid"?`, `show?`, `highlightPeakIds?`, `yHeadroom?`. `TracePlotInteraction = { onXDomain; onAddPeak?; onClickPeak?; onReset?; hitTolerancePx? }`. Build a model via the `modelFor(member)` pattern in `TracePlot.stories.tsx` (imports `realTraces`, `realMembers`, `SeriesMember`).
- **`PlotPeak`** (`src/print/plot/marks/PlotPeaks.tsx`): `{ id; q; intensity?; source: "auto"|"manual"; excluded?; predictedAbsent?; hot?; color?; label? }`.
- **`DetectorImage`** (`src/print/detector/DetectorImage.tsx`): `src: string|null`, `size: "thumb"|"full"`, `lutVariant?: DetectorLutVariant` ("neutral"|"warm"), `className`. Placeholder state → `data-testid="detector-image-placeholder"`; image state → `<canvas role="img" aria-label="Detector image">`.
- **`ThumbnailGallery`** (`src/print/components/ThumbnailGallery.tsx`): `exposures: GalleryExposure[]`, `selectedId?`, `onSelect?`, `size?` ("sm"|"lg"), `align?` ("start"|"center"), `className`. `GalleryExposure = { id; src: string|null; frameNo?; representative?; rejected? }`. Root `data-testid="thumbnail-gallery"`.
- **`Thumbnail`** (`src/print/ui/Thumbnail.tsx`): `src`, `frameNo?`, `representative?`, `rejected?`, `selected?`, `size?` ("sm"=62px|"lg"=70px), `onClick?`, `title?`, `className`. Root `data-testid="thumbnail"` (button).
- **`CombChart`** (`src/print/comb/CombChart.tsx`): `assigned: CombSeries[]`, `hovered?`, `leftover: number[]` (**required**), `hoveredQ?`, `onHoverQ?`, `maxWidth?`, `className?`. Wrapping div, no testid.
- **`ResidualChart`** (`src/print/comb/ResidualChart.tsx`): `assigned: CombSeries[]`, `hovered?`, `hoveredQ?`, `onHoverQ?`, `maxWidth?`, `className?` (NO `leftover`).
- **`CombLegend`** (`src/print/components/CombLegend.tsx`): `items?: ReadonlyArray<"auto"|"manual"|"predicted-absent"|"excluded">`, `className?`. **NOTE:** this is the peak-glyph vocab, NOT the mockup `.combs-legend` tooth vocab (predicted&observed / predicted-absent / leftover). Use as-is for now; the tooth-specific legend is a DEFERRED `CombLegend` refinement (record in ledger, do not block).
- **Fixtures:** `comb.fixtures.ts` exports `PN3M`, `IM3M` (`CombSeries`), `LEFTOVER` (`number[]`). `realTraces.ts` exports `realTraces: Record<number, Trace>` (keys 37/65/66/67/93). `realSeriesMembers.ts` exports `realMembers: SeriesMember[]`. Detector thumbs: `import thumb65 from "../fixtures/thumbs/65.png?url"` (37/65/66/67/93).
- **Token `--color-frame-edge`** exists (`styles.css:47`) → `bg-frame-edge` / `border-frame-edge` are guard-clean for the detector dark box. Also `--color-frame-tag` (`styles.css:51`).
- **Guard** (`npm run lint:design`): `src/print/components/**` is NOT exempt. Allowed: named-role type classes (`text-display`, `text-data`, `text-caption`, `text-kicker*`), token utilities (`bg-plate`, `bg-frame-edge`, `border-hair`, `text-ink-faint`, `text-print-accent`…), font utilities (`font-mono`, `font-semibold`), layout/placement (flex/grid/gap/p-*/mt-*/aspect-square/min-h-0/flex-1…). BANNED: `text-[…]`, `rounded-[…]`, raw `oklch(`/`rgba(`/quoted hex, side-stripe borders, inline `style` colour. `src/print/ui/**` IS exempt (Task 1's Thumbnail xs edit may use any appearance).

## Conventions (per-component recipe)

- Files: `src/print/components/<Name>.tsx` + `src/print/components/<Name>.stories.tsx`; tests in `test/print-components/<Name>.test.tsx`.
- Local `cx` helper per file: `function cx(...parts: Array<string | false | null | undefined>): string { return parts.filter(Boolean).join(" "); }`.
- Storybook idiom: `import type { Meta, StoryObj } from "@storybook/react-vite";` + loose `const meta: Meta<typeof X> = {…}`. `@storybook/test` NOT installed → `const noop = () => {};`. One default export per file. Interactive stories use a small `function XDemo()` wrapper with `useState` (mirror `ToolBar.stories.tsx`).
- Tests assert via `data-*` / roles / text / attributes, NEVER Tailwind/SVG class strings (single `className`-forwarding `.toContain()` is the one allowed exception).
- TDD per component: write failing test → run (fail) → implement → write story → `npm test -- print-components/<Name>` + `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json` green → commit named files.

---

## Task 1: `PanelHeader` (shared `.panel-head`)

The small uppercase header shared by `DetectorPanel` and `CombsPanel`: a faint kicker label on the left + an optional right tools slot. Derived from `.panel-head` / `.panel-h` (focus-plot mockup lines 190–197).

**Files:**
- Create: `src/print/components/PanelHeader.tsx`
- Create: `src/print/components/PanelHeader.stories.tsx`
- Test: `test/print-components/PanelHeader.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
// test/print-components/PanelHeader.test.tsx
import { render, screen } from "@testing-library/react";
import { PanelHeader } from "../../src/print/components/PanelHeader";

describe("PanelHeader", () => {
  it("renders the label text", () => {
    render(<PanelHeader label="Detector image" />);
    expect(screen.getByText("Detector image")).toBeInTheDocument();
  });
  it("renders a right-side tools slot when children are given", () => {
    render(
      <PanelHeader label="Reflections — comb">
        <button data-testid="tool">x</button>
      </PanelHeader>,
    );
    expect(screen.getByTestId("panel-header")).toBeInTheDocument();
    expect(screen.getByTestId("tool")).toBeInTheDocument();
  });
  it("omits the tools slot when no children", () => {
    render(<PanelHeader label="Detector image" />);
    expect(screen.queryByTestId("panel-header-tools")).not.toBeInTheDocument();
  });
  it("forwards a placement-only className", () => {
    render(<PanelHeader label="X" className="mb-5" />);
    expect(screen.getByTestId("panel-header").className).toContain("mb-5");
  });
});
```

- [ ] **Step 2: Run → fail.** `npm test -- print-components/PanelHeader` → FAIL (module not found).

- [ ] **Step 3: Implement**

```tsx
// src/print/components/PanelHeader.tsx
import type { ReactNode } from "react";
import { Kicker } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface PanelHeaderProps {
  /** The uppercase section label (.panel-h), e.g. "Detector image". */
  label: ReactNode;
  /** Optional right-side tools slot (exposure strip / comb-view toggle). */
  children?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function PanelHeader({ label, children, className }: PanelHeaderProps): JSX.Element {
  return (
    <div
      data-testid="panel-header"
      className={cx("flex items-center justify-between gap-2.5 mb-3", className)}
    >
      <Kicker tone="faint">{label}</Kicker>
      {children != null && <div data-testid="panel-header-tools">{children}</div>}
    </div>
  );
}
```

- [ ] **Step 4: Story**

```tsx
// src/print/components/PanelHeader.stories.tsx
import type { Meta, StoryObj } from "@storybook/react-vite";
import { PanelHeader } from "./PanelHeader";
import { SegmentedControl } from "../ui";

const meta: Meta<typeof PanelHeader> = {
  title: "components/PanelHeader",
  component: PanelHeader,
};
export default meta;
type Story = StoryObj<typeof meta>;

export const LabelOnly: Story = {
  args: { label: "Detector image" },
};

export const WithTools: Story = {
  render: () => (
    <PanelHeader label="Reflections — comb">
      <SegmentedControl
        size="xs"
        options={[
          { value: "comb", label: "comb" },
          { value: "resid", label: "indexing space" },
        ]}
        value="comb"
        onChange={() => {}}
        aria-label="comb view"
      />
    </PanelHeader>
  ),
};
```

- [ ] **Step 5: Verify.** `npm test -- print-components/PanelHeader` PASS · `npm run lint:design` clean · `npx tsc --noEmit -p tsconfig.build.json` clean.

- [ ] **Step 6: Commit**

```bash
git add src/print/components/PanelHeader.tsx src/print/components/PanelHeader.stories.tsx test/print-components/PanelHeader.test.tsx
git commit -m "feat(print): PanelHeader shared .panel-head (faint label + tools slot)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 2: `Thumbnail` / `ThumbnailGallery` `xs` size (refactor-on-contact)

The Focus DetectorPanel exposure strip is denser (≈30px, mockup `.e-thumb` lines 200–212) than the contact-sheet (`sm` 62px) / loupe (`lg` 70px) strips. Add an `xs` size to `Thumbnail` and thread it through `ThumbnailGallery`. Both files are in `src/print/ui/`? — **`Thumbnail` is in `src/print/ui/`, `ThumbnailGallery` is in `src/print/components/`.** `Thumbnail.tsx` (ui) is guard-exempt; `ThumbnailGallery.tsx` (components) must stay placement-only (it only forwards `size`, so no appearance added there).

**Files:**
- Modify: `src/print/ui/Thumbnail.tsx` (add `xs` to the `size` union + its size class)
- Modify: `src/print/components/ThumbnailGallery.tsx` (widen `size` prop type to include `xs`)
- Test: `test/print-ui/Thumbnail.test.tsx` (add an xs case) — VERIFY this path exists first (`ls test/print-ui/Thumbnail.test.tsx`); if the test lives elsewhere, adapt.

- [ ] **Step 1: Read** `src/print/ui/Thumbnail.tsx` fully. Note the exact `size` union, the size→dimension class map (find the `sm`/`lg` entry, e.g. a `Record<"sm"|"lg", string>` of width/height classes), and the existing `data-size` attribute if any.

- [ ] **Step 2: Write the failing test** (append to the existing Thumbnail test file; adapt to its current style)

```tsx
it("renders the xs (dense exposure-strip) size", () => {
  render(<Thumbnail src={null} size="xs" />);
  const t = screen.getByTestId("thumbnail");
  expect(t.dataset.size).toBe("xs"); // requires Thumbnail to set data-size={size}
});
```

If `Thumbnail` does not currently expose `data-size`, add `data-size={size}` to its root in Step 4 (it is a non-appearance test hook, allowed). If it already encodes size differently, assert that instead — do NOT assert class strings.

- [ ] **Step 3: Run → fail.** `npm test -- Thumbnail` → FAIL (`"xs"` not assignable / data-size undefined).

- [ ] **Step 4: Implement.** In `Thumbnail.tsx`: widen `size?: "sm" | "lg"` → `size?: "xs" | "sm" | "lg"`; add an `xs` entry to the size-class record (≈30px box; mirror the existing entries' shape, e.g. `xs: "w-[30px] h-[30px]"` — `w-[30px]` arbitrary sizing is permitted here because `ui/` is guard-exempt, but prefer a spacing-scale class if one matches, e.g. `size-7.5` is 30px if available; otherwise the arbitrary value is fine in `ui/`). Ensure `data-size={size}` is on the root button. In `ThumbnailGallery.tsx`: widen the `size?` prop type to `"xs" | "sm" | "lg"` (it just forwards `size` to each `Thumbnail`; no appearance change — stays placement-only).

- [ ] **Step 5: Verify.** `npm test -- Thumbnail` PASS · `npm run lint:design` clean (confirms `ThumbnailGallery` added no banned utility) · `npx tsc --noEmit -p tsconfig.build.json` clean.

- [ ] **Step 6: Commit**

```bash
git add src/print/ui/Thumbnail.tsx src/print/components/ThumbnailGallery.tsx test/print-ui/Thumbnail.test.tsx
git commit -m "feat(print): Thumbnail xs size for the dense Focus exposure strip

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 3: `TracePlate` (the hero trace plate)

The lifted hero plate: `PlateHeader` (kicker/title/subtitle + a `ToolBar` of scale toggle + Auto-fit + armed "+ Peak") over the `TracePlot`. Derived from `.plate` / `.plate-head` / `.tools` (focus-plot mockup lines 116–176, 430–448). Presentational: scale mode + add-peak armed state are props.

**Files:**
- Create: `src/print/components/TracePlate.tsx`
- Create: `src/print/components/TracePlate.stories.tsx`
- Test: `test/print-components/TracePlate.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
// test/print-components/TracePlate.test.tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { TracePlate } from "../../src/print/components/TracePlate";
import type { TraceModel } from "../../src/print/plot/TracePlot";

const model: TraceModel = {
  trace: { q: [0.02, 0.05, 0.1, 0.2], I: [10, 40, 20, 5], sigma: [1, 1, 1, 1] },
  peaks: [{ id: 0, q: 0.05, intensity: 40, source: "auto" }],
  phase: "Pn3m",
};

const base = {
  title: "Lipid 1-2 + LL37",
  trace: model,
  scale: "log" as const,
  onScaleChange: () => {},
};

describe("TracePlate", () => {
  it("renders the title and the trace plot region", () => {
    render(<TracePlate {...base} kicker="Integration" subtitle="smp_09" />);
    expect(screen.getByTestId("trace-plate")).toBeInTheDocument();
    expect(screen.getByText("Lipid 1-2 + LL37")).toBeInTheDocument();
    expect(screen.getByText("Integration")).toBeInTheDocument();
    expect(screen.getByText("smp_09")).toBeInTheDocument();
  });
  it("calls onScaleChange when the scale toggle is used", () => {
    const onScaleChange = vi.fn();
    render(<TracePlate {...base} onScaleChange={onScaleChange} />);
    fireEvent.click(screen.getByText("linear q"));
    expect(onScaleChange).toHaveBeenCalledWith("lin");
  });
  it("shows Auto-fit and fires onAutoFit", () => {
    const onAutoFit = vi.fn();
    render(<TracePlate {...base} onAutoFit={onAutoFit} />);
    fireEvent.click(screen.getByText("Auto-fit"));
    expect(onAutoFit).toHaveBeenCalled();
  });
  it("reflects the armed '+ Peak' toggle and fires onToggleAddPeak", () => {
    const onToggleAddPeak = vi.fn();
    render(<TracePlate {...base} addPeakArmed onToggleAddPeak={onToggleAddPeak} />);
    const peak = screen.getByText("+ Peak");
    expect(peak).toHaveAttribute("aria-pressed", "true");
    fireEvent.click(peak);
    expect(onToggleAddPeak).toHaveBeenCalled();
  });
  it("forwards a placement-only className", () => {
    render(<TracePlate {...base} className="mt-6" />);
    expect(screen.getByTestId("trace-plate").className).toContain("mt-6");
  });
});
```

- [ ] **Step 2: Run → fail.** `npm test -- print-components/TracePlate` → FAIL.

- [ ] **Step 3: Implement**

```tsx
// src/print/components/TracePlate.tsx
import type { ReactNode } from "react";
import { Card, Button, SegmentedControl } from "../ui";
import { PlateHeader } from "./PlateHeader";
import { ToolBar } from "./ToolBar";
import { TracePlot, type TraceModel, type TracePlotInteraction } from "../plot/TracePlot";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export type TraceScale = "log" | "lin";

export interface TracePlateProps {
  /** Accent eyebrow, e.g. "Integration". */
  kicker?: ReactNode;
  /** Serif plate title (the sample name). */
  title: ReactNode;
  /** Mono subtitle line (ids · facility · representative exposure). */
  subtitle?: ReactNode;
  /** The fully-built trace model (peaks + phase resolved by the caller). */
  trace: TraceModel;
  /** Current x-scale mode; the consumer owns it. */
  scale: TraceScale;
  onScaleChange: (next: TraceScale) => void;
  /** Auto-fit action. Omit → button hidden. */
  onAutoFit?: () => void;
  /** Add-peak armed (terracotta) toggle state + handler. */
  addPeakArmed?: boolean;
  onToggleAddPeak?: () => void;
  /** Forwarded TracePlot interaction (zoom / add / select). */
  interaction?: TracePlotInteraction | false;
  /** Plot height in px. Default 360. */
  plotHeight?: number;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function TracePlate({
  kicker,
  title,
  subtitle,
  trace,
  scale,
  onScaleChange,
  onAutoFit,
  addPeakArmed = false,
  onToggleAddPeak,
  interaction,
  plotHeight = 360,
  className,
}: TracePlateProps): JSX.Element {
  return (
    <Card
      as="section"
      elevated
      data-testid="trace-plate"
      className={cx("px-6 pt-5 pb-4", className)}
    >
      <PlateHeader kicker={kicker} title={title} subtitle={subtitle} as="h1">
        <ToolBar>
          <SegmentedControl
            options={[
              { value: "log", label: "log q" },
              { value: "lin", label: "linear q" },
            ]}
            value={scale}
            onChange={onScaleChange}
            aria-label="q scale"
          />
          {onAutoFit && (
            <Button variant="ghost" onClick={onAutoFit}>
              Auto-fit
            </Button>
          )}
          {onToggleAddPeak && (
            <Button armed={addPeakArmed} onClick={onToggleAddPeak}>
              + Peak
            </Button>
          )}
        </ToolBar>
      </PlateHeader>
      <TracePlot
        trace={trace}
        height={plotHeight}
        xType={scale === "log" ? "log" : "linear"}
        axes
        interaction={interaction}
        paperColor="var(--color-plate)"
        data-testid="trace-plate-plot"
        className="mt-2"
      />
    </Card>
  );
}
```

> If `Card`'s `as="section"` + extra `data-testid`/`className` typing complains, pass `data-testid` through `...rest` (Card spreads rest). If `xType`'s `ScaleType` literal is not `"log"`/`"linear"`, read `src/print/plot/projection.ts` for the exact union and map accordingly. `paperColor="var(--color-plate)"` is a CSS var STRING value (not a banned literal) — mirrors `CardFigure`.

- [ ] **Step 4: Run → pass.** `npm test -- print-components/TracePlate`.

- [ ] **Step 5: Story** (build the model with the `modelFor` pattern from `TracePlot.stories.tsx`)

```tsx
// src/print/components/TracePlate.stories.tsx
import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { TracePlate, type TraceScale } from "./TracePlate";
import type { TraceModel } from "../plot/TracePlot";
import type { PlotPeak } from "../plot/marks/PlotPeaks";
import { realTraces } from "../fixtures/realTraces";
import { realMembers } from "../fixtures/realSeriesMembers";
import type { SeriesMember } from "../../api";

function modelFor(member: SeriesMember): TraceModel {
  const trace = realTraces[member.exposure_id as number]!;
  const peaks: PlotPeak[] = (member.snapshot?.effective_peaks ?? []).map((p) => ({
    id: p.id,
    q: p.q,
    intensity: p.intensity,
    source: p.source,
  }));
  return { trace, peaks, phase: member.snapshot?.confirmed_index?.phase ?? null };
}

const heroModel = modelFor(realMembers[0]!); // exp 65

const meta: Meta<typeof TracePlate> = {
  title: "components/TracePlate",
  component: TracePlate,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function HeroDemo() {
  const [scale, setScale] = useState<TraceScale>("log");
  const [armed, setArmed] = useState(false);
  return (
    <div style={{ maxWidth: 1180 }}>
      <TracePlate
        kicker="Integration"
        title="Lipid 1-2 + LL37 1:0.5"
        subtitle="smp_09 · SSRL Apr 2026 · representative exposure smp_09_e03"
        trace={heroModel}
        scale={scale}
        onScaleChange={setScale}
        onAutoFit={() => {}}
        addPeakArmed={armed}
        onToggleAddPeak={() => setArmed((p) => !p)}
        interaction={{ onXDomain: () => {} }}
      />
    </div>
  );
}

export const Hero: Story = { render: () => <HeroDemo /> };
```

- [ ] **Step 6: Verify.** `npm test -- print-components/TracePlate` PASS · `npm run lint:design` clean · `npx tsc --noEmit -p tsconfig.build.json` clean.

- [ ] **Step 7: Commit**

```bash
git add src/print/components/TracePlate.tsx src/print/components/TracePlate.stories.tsx test/print-components/TracePlate.test.tsx
git commit -m "feat(print): TracePlate hero plate (PlateHeader+ToolBar+TracePlot)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 4: `DetectorPanel` (detector frame + exposure strip)

Flat panel: `PanelHeader` ("Detector image" + a tools slot for the exposure strip) over a dark `.det-box` framing `DetectorImage size="full"`, plus a hint line. Derived from the detector `.panel` (focus-plot mockup lines 184–227, 452–461). Presentational: the exposure strip + current `src` are supplied by the consumer (page owns `curExp`).

**Files:**
- Create: `src/print/components/DetectorPanel.tsx`
- Create: `src/print/components/DetectorPanel.stories.tsx`
- Test: `test/print-components/DetectorPanel.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
// test/print-components/DetectorPanel.test.tsx
import { render, screen } from "@testing-library/react";
import { DetectorPanel } from "../../src/print/components/DetectorPanel";

describe("DetectorPanel", () => {
  it("renders the default label and the detector frame", () => {
    render(<DetectorPanel src={null} />);
    expect(screen.getByTestId("detector-panel")).toBeInTheDocument();
    expect(screen.getByText("Detector image")).toBeInTheDocument();
    expect(screen.getByTestId("detector-frame")).toBeInTheDocument();
  });
  it("renders a custom label", () => {
    render(<DetectorPanel src={null} label="Real source" />);
    expect(screen.getByText("Real source")).toBeInTheDocument();
  });
  it("renders the header tools slot (exposure strip)", () => {
    render(<DetectorPanel src={null} tools={<div data-testid="expo">e</div>} />);
    expect(screen.getByTestId("expo")).toBeInTheDocument();
  });
  it("renders a hint line when provided", () => {
    render(<DetectorPanel src={null} hint="The real source." />);
    expect(screen.getByText("The real source.")).toBeInTheDocument();
  });
  it("forwards a placement-only className", () => {
    render(<DetectorPanel src={null} className="h-full" />);
    expect(screen.getByTestId("detector-panel").className).toContain("h-full");
  });
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement**

```tsx
// src/print/components/DetectorPanel.tsx
import type { ReactNode } from "react";
import { Card } from "../ui";
import { DetectorImage } from "../detector";
import type { DetectorLutVariant } from "../detector/detectorLut";
import { PanelHeader } from "./PanelHeader";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface DetectorPanelProps {
  /** Current exposure image URL; null → DetectorImage placeholder. */
  src: string | null;
  /** Section label. Default "Detector image". */
  label?: ReactNode;
  /** Header right slot — the exposure switcher (ThumbnailGallery). */
  tools?: ReactNode;
  /** Detector colormap; "neutral" lets the ring overlay own colour. */
  lutVariant?: DetectorLutVariant;
  /** Caption under the frame. */
  hint?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function DetectorPanel({
  src,
  label = "Detector image",
  tools,
  lutVariant,
  hint,
  className,
}: DetectorPanelProps): JSX.Element {
  return (
    <Card
      as="section"
      data-testid="detector-panel"
      className={cx("flex flex-col p-4", className)}
    >
      <PanelHeader label={label}>{tools}</PanelHeader>
      <div
        data-testid="detector-frame"
        className="bg-frame-edge border border-frame-edge rounded overflow-hidden aspect-square mx-auto w-full"
      >
        <DetectorImage src={src} size="full" lutVariant={lutVariant} />
      </div>
      {hint != null && (
        <div className="text-caption text-ink-faint mt-2.5">{hint}</div>
      )}
    </Card>
  );
}
```

> If `DetectorLutVariant` is not exported from `../detector/detectorLut`, import it from `../detector` (check `detector/index.ts`) or inline the `"neutral" | "warm"` union. If `DetectorImage size="full"` over-rotates inside a square frame, drop `aspect-square` and let the image's own intrinsic ratio drive height (verify visually in Storybook).

- [ ] **Step 4: Run → pass.**

- [ ] **Step 5: Story** (exposure strip via `ThumbnailGallery size="xs"` from Task 2)

```tsx
// src/print/components/DetectorPanel.stories.tsx
import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { DetectorPanel } from "./DetectorPanel";
import { ThumbnailGallery, type GalleryExposure } from "./ThumbnailGallery";
import thumb37 from "../fixtures/thumbs/37.png?url";
import thumb65 from "../fixtures/thumbs/65.png?url";
import thumb66 from "../fixtures/thumbs/66.png?url";
import thumb67 from "../fixtures/thumbs/67.png?url";

const SRCS = [thumb65, thumb66, thumb67, thumb37];
const EXPOSURES: GalleryExposure[] = SRCS.map((src, i) => ({
  id: i,
  src,
  frameNo: 65 + i,
  representative: i === 0,
}));

const meta: Meta<typeof DetectorPanel> = {
  title: "components/DetectorPanel",
  component: DetectorPanel,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function DetectorDemo() {
  const [cur, setCur] = useState(0);
  return (
    <div style={{ width: 372 }}>
      <DetectorPanel
        src={SRCS[cur]!}
        hint="The real source. Rings in phase colour; hover a peak, ring, or comb tooth — the triple lights up."
        tools={
          <ThumbnailGallery
            exposures={EXPOSURES}
            selectedId={cur}
            onSelect={setCur}
            size="xs"
          />
        }
      />
    </div>
  );
}

export const Default: Story = { render: () => <DetectorDemo /> };

export const Empty: Story = {
  args: { src: null, hint: "No exposure selected." },
};
```

- [ ] **Step 6: Verify.** `npm test -- print-components/DetectorPanel` PASS · `npm run lint:design` clean · `npx tsc --noEmit -p tsconfig.build.json` clean.

- [ ] **Step 7: Commit**

```bash
git add src/print/components/DetectorPanel.tsx src/print/components/DetectorPanel.stories.tsx test/print-components/DetectorPanel.test.tsx
git commit -m "feat(print): DetectorPanel (frame + exposure strip + hint)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 5: `CombsPanel` (comb / indexing-space toggle)

Flat panel: `PanelHeader` ("Reflections — comb" + a `SegmentedControl size="xs"` comb/indexing-space toggle) over a body that renders `CombChart` (with `leftover`) or `ResidualChart` per `view`, plus a `CombLegend` footer. Derived from the combs `.panel` (focus-plot mockup lines 229–260, 463–482). Presentational: `view` + `hoveredQ` are props.

**Files:**
- Create: `src/print/components/CombsPanel.tsx`
- Create: `src/print/components/CombsPanel.stories.tsx`
- Test: `test/print-components/CombsPanel.test.tsx`

- [ ] **Step 1: Write the failing test**

```tsx
// test/print-components/CombsPanel.test.tsx
import { render, screen, fireEvent } from "@testing-library/react";
import { CombsPanel } from "../../src/print/components/CombsPanel";
import { PN3M, IM3M, LEFTOVER } from "../../src/print/comb/comb.fixtures";

const base = {
  assigned: [PN3M, IM3M],
  leftover: LEFTOVER,
  view: "comb" as const,
  onViewChange: () => {},
};

describe("CombsPanel", () => {
  it("renders the label and the comb-view toggle", () => {
    render(<CombsPanel {...base} />);
    expect(screen.getByTestId("combs-panel")).toBeInTheDocument();
    expect(screen.getByText("Reflections — comb")).toBeInTheDocument();
    expect(screen.getByText("comb")).toBeInTheDocument();
    expect(screen.getByText("indexing space")).toBeInTheDocument();
  });
  it("renders the comb body and the legend in comb view", () => {
    render(<CombsPanel {...base} />);
    expect(screen.getByTestId("combs-body")).toBeInTheDocument();
  });
  it("calls onViewChange when the toggle is used", () => {
    const onViewChange = vi.fn();
    render(<CombsPanel {...base} onViewChange={onViewChange} />);
    fireEvent.click(screen.getByText("indexing space"));
    expect(onViewChange).toHaveBeenCalledWith("resid");
  });
  it("reflects the resid view via data-view", () => {
    render(<CombsPanel {...base} view="resid" />);
    expect(screen.getByTestId("combs-panel").dataset.view).toBe("resid");
  });
  it("forwards a placement-only className", () => {
    render(<CombsPanel {...base} className="h-full" />);
    expect(screen.getByTestId("combs-panel").className).toContain("h-full");
  });
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement**

```tsx
// src/print/components/CombsPanel.tsx
import type { ReactNode } from "react";
import { Card, SegmentedControl } from "../ui";
import { CombChart, ResidualChart, type CombSeries } from "../comb";
import { PanelHeader } from "./PanelHeader";
import { CombLegend } from "./CombLegend";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export type CombView = "comb" | "resid";

export interface CombsPanelProps {
  /** Assigned phases → comb rows. */
  assigned: CombSeries[];
  /** Leftover (unindexed) observed-peak q-values; comb view only. */
  leftover?: number[];
  /** Which view is showing; the consumer owns it. */
  view: CombView;
  onViewChange: (next: CombView) => void;
  /** Incoming q-link from the trace/detector. */
  hoveredQ?: number;
  onHoverQ?: (q?: number) => void;
  /** Section label. Default "Reflections — comb". */
  label?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function CombsPanel({
  assigned,
  leftover = [],
  view,
  onViewChange,
  hoveredQ,
  onHoverQ,
  label = "Reflections — comb",
  className,
}: CombsPanelProps): JSX.Element {
  return (
    <Card
      as="section"
      data-testid="combs-panel"
      data-view={view}
      className={cx("flex flex-col p-4", className)}
    >
      <PanelHeader label={label}>
        <SegmentedControl
          size="xs"
          options={[
            { value: "comb", label: "comb" },
            { value: "resid", label: "indexing space" },
          ]}
          value={view}
          onChange={onViewChange}
          aria-label="comb view"
        />
      </PanelHeader>
      <div data-testid="combs-body" className="flex-1 min-h-0 overflow-hidden">
        {view === "comb" ? (
          <CombChart
            assigned={assigned}
            leftover={leftover}
            hoveredQ={hoveredQ}
            onHoverQ={onHoverQ}
          />
        ) : (
          <ResidualChart assigned={assigned} hoveredQ={hoveredQ} onHoverQ={onHoverQ} />
        )}
      </div>
      <CombLegend className="mt-2.5" />
    </Card>
  );
}
```

> `CombView` ("comb"|"resid") must match the `SegmentedControl` option `value`s exactly so `onChange`'s `T` infers to `CombView` (if TS widens to `string`, annotate `options` as `SegmentOption<CombView>[]` — import the type from `../ui`). `CombLegend`'s default item set is the peak-glyph vocab — acceptable for this slice; the exact tooth-vocab legend (predicted&observed / predicted-absent / leftover) is a DEFERRED `CombLegend` refinement (record in Task 6 ledger).

- [ ] **Step 4: Run → pass.**

- [ ] **Step 5: Story**

```tsx
// src/print/components/CombsPanel.stories.tsx
import { useState } from "react";
import type { Meta, StoryObj } from "@storybook/react-vite";
import { CombsPanel, type CombView } from "./CombsPanel";
import { PN3M, IM3M, LEFTOVER } from "../comb/comb.fixtures";

const meta: Meta<typeof CombsPanel> = {
  title: "components/CombsPanel",
  component: CombsPanel,
  parameters: { layout: "padded" },
};
export default meta;
type Story = StoryObj<typeof meta>;

function CombsDemo() {
  const [view, setView] = useState<CombView>("comb");
  const [hoveredQ, setHoveredQ] = useState<number | undefined>(undefined);
  return (
    <div style={{ width: 720, height: 416 }}>
      <CombsPanel
        assigned={[PN3M, IM3M]}
        leftover={LEFTOVER}
        view={view}
        onViewChange={setView}
        hoveredQ={hoveredQ}
        onHoverQ={setHoveredQ}
        className="h-full"
      />
    </div>
  );
}

export const Default: Story = { render: () => <CombsDemo /> };
```

- [ ] **Step 6: Verify.** `npm test -- print-components/CombsPanel` PASS · `npm run lint:design` clean · `npx tsc --noEmit -p tsconfig.build.json` clean.

- [ ] **Step 7: Commit**

```bash
git add src/print/components/CombsPanel.tsx src/print/components/CombsPanel.stories.tsx test/print-components/CombsPanel.test.tsx
git commit -m "feat(print): CombsPanel (comb/indexing-space toggle + legend)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```

---

## Task 6: Ledger + memory + batch gate

**Files:**
- Modify: `docs/greenfield-component-ledger.md`
- Modify: `/Users/me/.claude/projects/-Users-me-projects-Himalaya-jl/memory/project_greenfield_composite_layer.md`
- Modify: `/Users/me/.claude/projects/-Users-me-projects-Himalaya-jl/memory/MEMORY.md`

- [ ] **Step 1: Full batch gate.** Run (capture once): `npx tsc --noEmit -p tsconfig.build.json` (exit 0) · `npm run lint:design` (clean) · `npm test -- print` (or full `npm test`) all green · `npm run build-storybook` (exit 0).

- [ ] **Step 2: Visual fidelity.** Serve `storybook-static` and screenshot `components/TracePlate` (Hero), `components/DetectorPanel` (Default), `components/CombsPanel` (Default, both views) against `2026-05-29-focus-plot.html`.

- [ ] **Step 3: Ledger.** Flip `TracePlate`, `DetectorPanel`, `CombsPanel` rows to ✅ (Batch 7); add `PanelHeader` as a new tier-1/shared composite ✅; bump the Layer-3 coverage line (5 → 8 ✅). Add decisions-registry entries: (a) "Focus lower panels share `PanelHeader`; panels are presentational — scale/exposure/view/hovered-q state lifted to the consumer (Batch 7)"; (b) "`CombsPanel` footer uses the existing peak-glyph `CombLegend`; the mockup's comb-tooth legend vocab (predicted&observed / leftover) is a DEFERRED `CombLegend` refinement"; (c) "`Thumbnail` gained an `xs` (≈30px) size for the dense Focus exposure strip (refactor-on-contact, Batch 7)". Strike the Focus-plates items from the frontier; note the Focus page (Layer 4) is now assemblable.

- [ ] **Step 4: Memory.** Update `project_greenfield_composite_layer.md` (frontmatter description + a BATCH 7 paragraph + rewrite NEXT) and the `MEMORY.md` index line (L3 8✅, new gotchas: presentational-panels-lift-state + CombLegend-tooth-vocab-deferred + Thumbnail-xs).

- [ ] **Step 5: Commit**

```bash
git add docs/greenfield-component-ledger.md
git commit -m "docs(print): ledger — Batch 7 Focus-plates slice (TracePlate/DetectorPanel/CombsPanel/PanelHeader)

Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>"
```
(Memory files are outside the repo — no commit.)

---

## Verification (whole batch)

- `npx tsc --noEmit -p tsconfig.build.json` exit 0.
- `npm run lint:design` clean (proves all four new `components/` files are placement-only).
- `npm test` green (new suites: PanelHeader, TracePlate, DetectorPanel, CombsPanel + the Thumbnail xs case).
- `npm run build-storybook` exit 0; `components/{TracePlate,DetectorPanel,CombsPanel,PanelHeader}` stories render and match the focus-plot mockup.
- Working tree: only the named files committed; `src/bones/*` and untracked plan docs untouched. Branch `worktree-greenfield-ui-rebuild` remains unmerged.
