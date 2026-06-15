# Greenfield Batch 9 — Series reading slice

> **For agentic workers:** REQUIRED SUB-SKILL: superpowers:subagent-driven-development — fresh implementer per task + frontend-reviewer + fix loop. Tasks run SEQUENTIALLY (each commits; one git actor).
> Commit trailer (every commit, exact last line): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Commit ONLY the named files (never `git add -A`/`.`). NEVER stage `src/bones/registry.ts`, `src/bones/contact-sheet.bones.json`, or anything under `docs/superpowers/plans/`.
> Typecheck: `npx tsc --noEmit -p tsconfig.build.json`. Guard: `npm run lint:design`. Work from `packages/HimalayaUI/frontend`.

**Goal:** Build the series-plot *reading surface* composites — `SeriesMemberRow`, `MemberList`, `ReadingPanel`, `SeriesPlate` — plus the `Swatch` coexistence/empty refactor-on-contact.

**Architecture:** Bottom-up, placement-only composites over existing primitives + the `WaterfallChart` renderer. Panels are presentational: every datum is a prop, cross-component hover state is lifted to the page (the story simulates the page). Derived from `docs/redesign-mockups/2026-05-29-series-plot.html`.

**Out of scope (recorded deviations):** waterfall/heatmap rep toggle (HeatmapChart out-of-scope); vertical phase-strip companion (needs an unbuilt L1 renderer aligned to WaterfallChart's internal row geometry).

---

## Mockup reference (`2026-05-29-series-plot.html`)

- `.plate` (lines 91–96, 264–297) → **SeriesPlate**: `.plate-head` (kicker+serif h1+mono h-sub) + `.tools` (two segs) + `.figure-wrap` (#waterfall + #phasestrip) + `.plate-foot` (legend + mono offset note).
- `.reading` / `.rd-row` / `.rd-coex` (lines 178–192, 646–675) → **ReadingPanel** + existing **ReadingRow**.
- `.member` / `.m-*` (lines 164–175, 677–715) → **MemberList** + new **SeriesMemberRow**.

Type-scale snaps (named roles only, no arbitrary `text-[…]`): `.m-name` 12px/600 → `text-meta font-semibold`; `.m-var`/`.m-data` mono → `text-data`; `.plate-foot` 11px → `text-meta`; mono offset note → `text-data`. Plate h1 27px = `text-display` (via PlateHeader). `.reading` box = `Card`.

---

## Task 1: Extend `Swatch` (coexistence + empty modes) — refactor-on-contact (ui/, guard-exempt)

**Files:**
- Modify: `src/print/ui/Swatch.tsx`
- Test: `test/print-ui/Swatch.test.tsx` (extend existing)

The series-plot member swatch is solid (1 phase), a 135° two-phase **gradient** (coexistence), or a transparent **dashed** cell (form factor), at 11px. `Swatch` is the only place these appearance modes can live (placement-only consumers can't author a gradient). The stub comment at the file foot already anticipates this.

- [ ] **Step 1: Extend the Swatch test** (`test/print-ui/Swatch.test.tsx`) — read the existing file first, then add:

```tsx
import { phaseColor } from "../../src/phases";

it("renders a coexistence gradient blending both phases", () => {
  const { getByTestId } = render(<div data-testid="w"><Swatch phase="Pn3m" coexistWith="Lamellar" shape="circle" /></div>);
  const sw = getByTestId("w").firstChild as HTMLElement;
  expect(sw.getAttribute("data-coexist")).toBe("Lamellar");
  // The DOM canonicalizes the OKLCH substrings on write, so compare against an
  // identically-built reference (round-trip), as the solid-color test does.
  const ref = document.createElement("span");
  ref.style.background = `linear-gradient(135deg, ${phaseColor("Pn3m")} 48%, ${phaseColor("Lamellar")} 52%)`;
  expect(sw.style.background).toBe(ref.style.background);
});

it("renders an empty (form-factor) swatch with no phase fill", () => {
  const { getByTestId } = render(<div data-testid="w"><Swatch phase="Pn3m" empty /></div>);
  const sw = getByTestId("w").firstChild as HTMLElement;
  expect(sw.getAttribute("data-empty")).toBe("true");
  // empty wins over phase fill — no solid phase color as the background
  expect(sw.style.background).not.toContain(phaseColor("Pn3m"));
});

it("size md is 11px", () => {
  const { getByTestId } = render(<div data-testid="w"><Swatch phase="Pn3m" size="md" /></div>);
  const sw = getByTestId("w").firstChild as HTMLElement;
  expect(sw.className).toContain("w-[11px]");
});
```

(Keep/adjust any existing assertions — defaults must stay `sm`/9px/solid so ReadingRow + builder MemberRow are unchanged.)

- [ ] **Step 2: Run → fail** — `npm test -- print-ui/Swatch`

- [ ] **Step 3: Rewrite `Swatch.tsx`:**

```tsx
import { phaseColor } from "../../phases";

export type SwatchShape = "square" | "circle";
export type SwatchSize = "sm" | "md";

interface SwatchProps {
  /** Phase whose color fills the swatch (via phaseColor). Ignored when `empty`. */
  phase: string;
  /** "square" (rounded color chip, default) | "circle" (reading-row dot). */
  shape?: SwatchShape;
  /** "sm" 9px (default) | "md" 11px (series-plot member). */
  size?: SwatchSize;
  /** Second phase → a 135° gradient blending both (coexistence). */
  coexistWith?: string;
  /** Form factor: transparent fill + dashed hairline, no phase color. */
  empty?: boolean;
  /** PLACEMENT-ONLY. */
  className?: string;
}

const shapeClass: Record<SwatchShape, string> = {
  square: "rounded-sm",
  circle: "rounded-full",
};
const sizeClass: Record<SwatchSize, string> = {
  sm: "h-[9px] w-[9px]",
  md: "h-[11px] w-[11px]",
};

export function Swatch({
  phase,
  shape = "square",
  size = "sm",
  coexistWith,
  empty = false,
  className = "",
}: SwatchProps): JSX.Element {
  const background = empty
    ? "transparent"
    : coexistWith
      ? `linear-gradient(135deg, ${phaseColor(phase)} 48%, ${phaseColor(coexistWith)} 52%)`
      : phaseColor(phase);
  return (
    <span
      data-swatch
      data-phase={phase}
      data-shape={shape}
      {...(coexistWith ? { "data-coexist": coexistWith } : {})}
      {...(empty ? { "data-empty": "true" } : {})}
      aria-hidden="true"
      className={`inline-block shrink-0 ${sizeClass[size]} ${shapeClass[shape]} ${empty ? "border-[1.4px] border-dashed border-hair-strong" : ""} ${className}`.trim()}
      style={{ background }}
    />
  );
}
```

(`border-[1.4px]` is an arbitrary **border-width** — geometry, not color; allowed. The color is the token `border-hair-strong`. ui/ is guard-exempt regardless, but keep it clean.)

- [ ] **Step 4: export the new type** — add `SwatchSize` to the `Swatch` re-export in `src/print/ui/index.ts` (alongside `SwatchShape`).

- [ ] **Step 5: Gate** — `npm test -- print-ui/Swatch` green; `npm run lint:design` clean; `npx tsc --noEmit -p tsconfig.build.json` clean.

- [ ] **Step 6: Commit** — `git add src/print/ui/Swatch.tsx src/print/ui/index.ts test/print-ui/Swatch.test.tsx` then commit (`feat(print/ui): Swatch coexistence gradient + empty + size modes`).

---

## Task 2: `SeriesMemberRow` leaf (components/, placement-only)

**Files:**
- Create: `src/print/components/SeriesMemberRow.tsx`, `src/print/components/SeriesMemberRow.stories.tsx`
- Test: `test/print-components/SeriesMemberRow.test.tsx`

The series-plot `.member` row (NOT the builder `.trow` — that's the existing `MemberRow`). Swatch decodes the phase(s); names are colored per phase and joined by " + "; the variable value sits at the row's right; a mono data line carries lattice + q₁. Fully presentational.

- [ ] **Step 1: Test** (`test/print-components/SeriesMemberRow.test.tsx`):

```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SeriesMemberRow } from "../../src/print/components/SeriesMemberRow";

describe("SeriesMemberRow", () => {
  it("renders a single-phase member: solid swatch, one colored name, variable + data line", () => {
    const { getByTestId, getByText, container } = render(
      <SeriesMemberRow phases={["Pn3m"]} variableValue="1:0" dataLine="a = 205 Å · q₁ 0.061 Å⁻¹" />,
    );
    expect(getByTestId("series-member-row")).toBeTruthy();
    expect(getByText("Pn3m")).toBeTruthy();
    expect(getByText("1:0")).toBeTruthy();
    expect(getByText("a = 205 Å · q₁ 0.061 Å⁻¹")).toBeTruthy();
    const sw = container.querySelector("[data-swatch]")!;
    expect(sw.getAttribute("data-coexist")).toBeNull();
    expect(sw.getAttribute("data-empty")).toBeNull();
  });

  it("renders a coexistence member: gradient swatch + both names", () => {
    const { getByText, container } = render(
      <SeriesMemberRow phases={["Pn3m", "Lamellar"]} variableValue="1:0.5" dataLine="a 195 · d 60 Å · q₁ 0.057 Å⁻¹" />,
    );
    expect(getByText("Pn3m")).toBeTruthy();
    expect(getByText("Lamellar")).toBeTruthy();
    expect(container.querySelector("[data-swatch]")!.getAttribute("data-coexist")).toBe("Lamellar");
  });

  it("renders a form-factor member: empty swatch + 'Form factor'", () => {
    const { getByText, container } = render(
      <SeriesMemberRow phases={[]} variableValue="1:1.5" dataLine="no Bragg peaks · q₁ —" />,
    );
    expect(getByText("Form factor")).toBeTruthy();
    expect(container.querySelector("[data-swatch]")!.getAttribute("data-empty")).toBe("true");
  });

  it("marks hot state and fires hover handlers", () => {
    const onHover = vi.fn(), onLeave = vi.fn();
    const { getByTestId } = render(
      <SeriesMemberRow phases={["Pn3m"]} variableValue="1:0" dataLine="x" hot onHover={onHover} onLeave={onLeave} />,
    );
    const row = getByTestId("series-member-row");
    expect(row.getAttribute("data-hot")).toBe("true");
    fireEvent.mouseEnter(row); fireEvent.mouseLeave(row);
    expect(onHover).toHaveBeenCalledTimes(1);
    expect(onLeave).toHaveBeenCalledTimes(1);
  });
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement** (`src/print/components/SeriesMemberRow.tsx`):

```tsx
import type { ReactNode } from "react";
import { Swatch, PhaseLabel } from "../ui";

export interface SeriesMemberRowProps {
  /** Dominant + any coexisting phase names, display order. `[]` = form factor. */
  phases: string[];
  /** The variable value (e.g. "1:0.5"), right-aligned, mono. */
  variableValue: ReactNode;
  /** Page-derived lattice + first-peak line (e.g. "a = 195 Å · q₁ 0.057 Å⁻¹"). */
  dataLine: ReactNode;
  /** Page-owned highlight (synced with the waterfall hot row). */
  hot?: boolean;
  onHover?: () => void;
  onLeave?: () => void;
  /** PLACEMENT-ONLY. Appended last. */
  className?: string;
}

// The series-plot reading "member" (mockup `.member`) — a self-decoding row:
// swatch shows phase(s), names are colored per phase, the variable sits right,
// and a mono data line carries lattice + q₁. Presentational: `hot` + handlers
// are page-owned so hovering here / the waterfall / nothing stays in sync.
export function SeriesMemberRow({
  phases,
  variableValue,
  dataLine,
  hot = false,
  onHover,
  onLeave,
  className,
}: SeriesMemberRowProps): JSX.Element {
  const formFactor = phases.length === 0;
  return (
    <div
      data-testid="series-member-row"
      {...(hot ? { "data-hot": "true" } : {})}
      onMouseEnter={onHover}
      onMouseLeave={onLeave}
      className={`flex items-center gap-[9px] px-[9px] py-[7px] rounded cursor-pointer border ${hot ? "bg-plate border-hair" : "border-transparent hover:bg-plate"}${className ? ` ${className}` : ""}`}
    >
      {formFactor ? (
        <Swatch phase="" empty size="md" shape="circle" />
      ) : phases.length >= 2 ? (
        <Swatch phase={phases[0]!} coexistWith={phases[1]!} size="md" shape="circle" />
      ) : (
        <Swatch phase={phases[0]!} size="md" shape="circle" />
      )}
      <div className="flex-1 min-w-0 flex flex-col gap-px">
        <div className="flex items-baseline justify-between gap-2">
          <span className="inline-flex items-baseline min-w-0">
            {formFactor ? (
              <span className="text-meta font-semibold text-ink-faint">Form factor</span>
            ) : (
              phases.map((p, i) => (
                <span key={p} className="inline-flex items-baseline">
                  {i > 0 && <span className="text-meta text-ink-faint mx-0.5"> + </span>}
                  <PhaseLabel phase={p} className="text-meta font-semibold">{p}</PhaseLabel>
                </span>
              ))
            )}
          </span>
          <span className="text-data text-ink-soft shrink-0">{variableValue}</span>
        </div>
        <div className="text-data text-ink-faint truncate">{dataLine}</div>
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Story** (`SeriesMemberRow.stories.tsx`) — `title: "components/SeriesMemberRow"`, one default export. A `Gallery` story rendering a solid, a coexistence, and a form-factor row stacked (gap-0.5) so all three swatch modes are visible.

- [ ] **Step 5: Gate** — `npm test -- print-components/SeriesMemberRow`; lint:design; tsc build config.

- [ ] **Step 6: Commit** — `git add src/print/components/SeriesMemberRow.tsx src/print/components/SeriesMemberRow.stories.tsx test/print-components/SeriesMemberRow.test.tsx`.

---

## Task 3: `MemberList` panel (components/)

**Files:**
- Create: `src/print/components/MemberList.tsx`, `src/print/components/MemberList.stories.tsx`
- Test: `test/print-components/MemberList.test.tsx`

Data-driven list of `SeriesMemberRow`s with controlled hover. Renders members in the order given (the page reverses so list-top = waterfall-top).

- [ ] **Step 1: Test:**

```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { MemberList, type MemberDatum } from "../../src/print/components/MemberList";

const members: MemberDatum[] = [
  { key: "a", phases: ["Pn3m"], variableValue: "1:0", dataLine: "a = 205 Å · q₁ 0.061 Å⁻¹" },
  { key: "b", phases: ["Pn3m", "Lamellar"], variableValue: "1:0.5", dataLine: "a 195 · d 60 Å" },
  { key: "c", phases: [], variableValue: "1:1.5", dataLine: "no Bragg peaks · q₁ —" },
];

describe("MemberList", () => {
  it("renders one row per member in order", () => {
    const { getAllByTestId } = render(<MemberList members={members} />);
    expect(getAllByTestId("series-member-row")).toHaveLength(3);
  });
  it("marks the hovered key hot and reports hover/leave", () => {
    const onHoverMember = vi.fn();
    const { getAllByTestId } = render(
      <MemberList members={members} hoveredKey="b" onHoverMember={onHoverMember} />,
    );
    const rows = getAllByTestId("series-member-row");
    expect(rows[1]!.getAttribute("data-hot")).toBe("true");
    expect(rows[0]!.getAttribute("data-hot")).toBeNull();
    fireEvent.mouseEnter(rows[0]!);
    expect(onHoverMember).toHaveBeenCalledWith("a");
    fireEvent.mouseLeave(rows[0]!);
    expect(onHoverMember).toHaveBeenCalledWith(undefined);
  });
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement:**

```tsx
import type { ReactNode } from "react";
import { SeriesMemberRow } from "./SeriesMemberRow";

export interface MemberDatum {
  key: string;
  phases: string[];
  variableValue: ReactNode;
  dataLine: ReactNode;
}

export interface MemberListProps {
  /** Members in display order (page reverses so top = high variable). */
  members: MemberDatum[];
  /** Controlled hot key, synced with the waterfall hot row. */
  hoveredKey?: string;
  onHoverMember?: (key?: string) => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function MemberList({ members, hoveredKey, onHoverMember, className }: MemberListProps): JSX.Element {
  return (
    <div data-testid="member-list" className={`flex flex-col gap-0.5${className ? ` ${className}` : ""}`}>
      {members.map((m) => (
        <SeriesMemberRow
          key={m.key}
          phases={m.phases}
          variableValue={m.variableValue}
          dataLine={m.dataLine}
          hot={hoveredKey === m.key}
          {...(onHoverMember ? { onHover: () => onHoverMember(m.key), onLeave: () => onHoverMember(undefined) } : {})}
        />
      ))}
    </div>
  );
}
```

- [ ] **Step 4: Story** — `title: "components/MemberList"`. A stateful `OpenDemo` holding `hoveredKey` so hovering a row highlights it; data = the mockup LL37 titration members (reversed: 1:1.5 form-factor at top → 1:0 Pn3m at bottom). Keep data inline in the story (page-owned, mirrors CustomIndexModal.stories pattern).

- [ ] **Step 5: Gate.**

- [ ] **Step 6: Commit** — `git add src/print/components/MemberList.tsx src/print/components/MemberList.stories.tsx test/print-components/MemberList.test.tsx`.

---

## Task 4: `ReadingPanel` panel (components/)

**Files:**
- Create: `src/print/components/ReadingPanel.tsx`, `src/print/components/ReadingPanel.stories.tsx`
- Test: `test/print-components/ReadingPanel.test.tsx`

The derived "phases present" box: `Card` wrapping `ReadingRow`×N + optional coexistence / form-factor footer lines (`.rd-coex`: mono ink-soft, hairline top).

- [ ] **Step 1: Test:**

```tsx
import { render } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ReadingPanel, type ReadingDatum } from "../../src/print/components/ReadingPanel";

const readings: ReadingDatum[] = [
  { phase: "Pn3m", span: "1:0 → 1:0.5", lattice: "a 205 → 195 Å" },
  { phase: "Lamellar", span: "1:0.5 → 1:1", lattice: "d 60 → 55 Å" },
];

describe("ReadingPanel", () => {
  it("renders one ReadingRow per reading", () => {
    const { getAllByTestId } = render(<ReadingPanel readings={readings} />);
    expect(getAllByTestId("reading-row")).toHaveLength(2);
  });
  it("renders coexistence and form-factor notes when provided", () => {
    const { getByText } = render(
      <ReadingPanel readings={readings} coexistenceNote="coexistence at 1:0.5" formFactorNote="form factor only at 1:1.5" />,
    );
    expect(getByText("coexistence at 1:0.5")).toBeTruthy();
    expect(getByText("form factor only at 1:1.5")).toBeTruthy();
  });
  it("omits notes when not provided", () => {
    const { queryByTestId } = render(<ReadingPanel readings={readings} />);
    expect(queryByTestId("reading-coex")).toBeNull();
  });
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement:**

```tsx
import type { ReactNode } from "react";
import { Card } from "../ui";
import { ReadingRow } from "./ReadingRow";

export interface ReadingDatum {
  phase: string;
  span: ReactNode;
  lattice: ReactNode;
}

export interface ReadingPanelProps {
  readings: ReadingDatum[];
  /** e.g. "coexistence at 1:0.5". */
  coexistenceNote?: ReactNode;
  /** e.g. "form factor only at 1:1.5". */
  formFactorNote?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function ReadingPanel({ readings, coexistenceNote, formFactorNote, className }: ReadingPanelProps): JSX.Element {
  return (
    <Card data-testid="reading-panel" className={`px-[13px] py-[11px] flex flex-col gap-[9px]${className ? ` ${className}` : ""}`}>
      {readings.map((r, i) => (
        <ReadingRow key={i} phase={r.phase} span={r.span} lattice={r.lattice} />
      ))}
      {coexistenceNote != null && (
        <div data-testid="reading-coex" className="text-data text-ink-soft mt-0.5 pt-2 border-t border-hair">{coexistenceNote}</div>
      )}
      {formFactorNote != null && (
        <div data-testid="reading-ff" className="text-data text-ink-soft mt-0.5 pt-2 border-t border-hair">{formFactorNote}</div>
      )}
    </Card>
  );
}
```

(Verify `Card` forwards `data-testid` + arbitrary className and renders a plate/hair box without forcing elevation — TracePlate uses `<Card elevated …>`; here omit `elevated`. If Card's base padding conflicts, the placement padding above overrides; confirm by reading `Card.tsx`.)

- [ ] **Step 4: Story** — `title: "components/ReadingPanel"`. Data = mockup-derived Pn3m + Lamellar readings + both notes.

- [ ] **Step 5: Gate.**

- [ ] **Step 6: Commit** — `git add src/print/components/ReadingPanel.tsx src/print/components/ReadingPanel.stories.tsx test/print-components/ReadingPanel.test.tsx`.

---

## Task 5: `SeriesPlate` panel (components/)

**Files:**
- Create: `src/print/components/SeriesPlate.tsx`, `src/print/components/SeriesPlate.stories.tsx`
- Test: `test/print-components/SeriesPlate.test.tsx`

Mirror `TracePlate`: `Card elevated` + `PlateHeader`(as h1) + `ToolBar`(scale SegmentedControl only) + `WaterfallChart` + optional foot. Heatmap toggle omitted; phase-strip companion deferred.

- [ ] **Step 1: Test:**

```tsx
import { render, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { SeriesPlate } from "../../src/print/components/SeriesPlate";
import { TRANSITION } from "../../src/print/waterfall/waterfall.fixtures";

describe("SeriesPlate", () => {
  it("renders header, waterfall, scale toggle, and the foot", () => {
    const { getByTestId, getByText, getAllByRole } = render(
      <SeriesPlate kicker="Series" title="LL37 titration" subtitle="7 exposures" rows={TRANSITION}
        scale="log" onScaleChange={() => {}} legendPhases={["Pn3m", "Lamellar"]}
        footHint="peaks are light anchors" footNote="offset ×1.0 · log I" />,
    );
    expect(getByTestId("series-plate")).toBeTruthy();
    expect(getByText("LL37 titration")).toBeTruthy();
    expect(getByTestId("waterfall")).toBeTruthy();
    expect(getAllByRole("button").some((b) => b.textContent === "log q")).toBe(true);
    expect(getByText("offset ×1.0 · log I")).toBeTruthy();
  });
  it("does not render a heatmap toggle", () => {
    const { queryByText } = render(
      <SeriesPlate title="x" rows={TRANSITION} scale="log" onScaleChange={() => {}} />,
    );
    expect(queryByText(/heatmap/i)).toBeNull();
  });
  it("wires the scale toggle", () => {
    const onScaleChange = vi.fn();
    const { getByText } = render(
      <SeriesPlate title="x" rows={TRANSITION} scale="log" onScaleChange={onScaleChange} />,
    );
    fireEvent.click(getByText("linear q"));
    expect(onScaleChange).toHaveBeenCalledWith("lin");
  });
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement:**

```tsx
import type { ReactNode } from "react";
import { Card, SegmentedControl, Swatch } from "../ui";
import { PlateHeader } from "./PlateHeader";
import { ToolBar } from "./ToolBar";
import { WaterfallChart } from "../waterfall";
import type { WaterfallRow } from "../waterfall/waterfallModel";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export type SeriesScale = "log" | "lin";

export interface SeriesPlateProps {
  kicker?: ReactNode;
  title: ReactNode;
  subtitle?: ReactNode;
  /** Member rows low→high (rendered bottom-up by WaterfallChart). */
  rows: WaterfallRow[];
  scale: SeriesScale;
  onScaleChange: (next: SeriesScale) => void;
  /** Controlled hot row + cursor q, lifted to the page (synced with MemberList). */
  hoveredKey?: string;
  onHoverRow?: (key?: string) => void;
  hoveredQ?: number;
  onHoverQ?: (q?: number) => void;
  /** Foot legend phase dots (e.g. ["Pn3m","Lamellar"]). */
  legendPhases?: string[];
  /** Foot hint (e.g. "peaks are light anchors — hover a trace to read its indexing"). */
  footHint?: ReactNode;
  /** Mono right-aligned foot note (e.g. "offset ×1.0 · log I"). */
  footNote?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function SeriesPlate({
  kicker, title, subtitle, rows, scale, onScaleChange,
  hoveredKey, onHoverRow, hoveredQ, onHoverQ,
  legendPhases, footHint, footNote, className,
}: SeriesPlateProps): JSX.Element {
  const hasFoot = (legendPhases && legendPhases.length > 0) || footHint != null || footNote != null;
  return (
    <Card as="section" elevated data-testid="series-plate" className={cx("px-6 pt-5 pb-[18px]", className)}>
      <PlateHeader kicker={kicker} title={title} subtitle={subtitle} as="h1">
        <ToolBar>
          <SegmentedControl
            options={[{ value: "log", label: "log q" }, { value: "lin", label: "linear q" }]}
            value={scale}
            onChange={onScaleChange}
            aria-label="q scale"
          />
        </ToolBar>
      </PlateHeader>
      <WaterfallChart
        rows={rows}
        xType={scale === "log" ? "log" : "linear"}
        className="mt-2"
        {...(hoveredKey !== undefined ? { hoveredKey } : {})}
        {...(onHoverRow ? { onHoverRow } : {})}
        {...(hoveredQ !== undefined ? { hoveredQ } : {})}
        {...(onHoverQ ? { onHoverQ } : {})}
      />
      {hasFoot && (
        <div data-testid="series-plate-foot" className="mt-3 pt-[11px] border-t border-hair flex items-center justify-between text-meta text-ink-faint">
          <div className="flex items-center gap-[14px]">
            {legendPhases?.map((p) => (
              <span key={p} className="inline-flex items-center gap-[5px]">
                <Swatch phase={p} shape="circle" />
                <span>{p}</span>
              </span>
            ))}
            {footHint != null && <span>{footHint}</span>}
          </div>
          {footNote != null && <span className="text-data">{footNote}</span>}
        </div>
      )}
    </Card>
  );
}
```

(Verify `WaterfallChart` is re-exported from `../waterfall` (index.ts); if not, import from `../waterfall/WaterfallChart`. Confirm `SegmentedControl` infers the `"log"|"lin"` union from `options` — `TracePlate` proves the pattern.)

- [ ] **Step 4: Story** (`SeriesPlate.stories.tsx`) — `title: "components/SeriesPlate"`. A stateful demo holding `scale` so the toggle works; `rows={TRANSITION}` (from `waterfall.fixtures`); `legendPhases`, `footHint`, `footNote` from the mockup foot. Wrap in a max-w container so the plate reads at plate width.

- [ ] **Step 5: Gate.**

- [ ] **Step 6: Commit** — `git add src/print/components/SeriesPlate.tsx src/print/components/SeriesPlate.stories.tsx test/print-components/SeriesPlate.test.tsx`.

---

## Task 6: Final gate + ledger/memory reconciliation

- [ ] Full gate: `npm test -- print-ui print-components` green; `npm run lint:design` clean; `npx tsc --noEmit -p tsconfig.build.json` clean; `npm run build-storybook` exit 0.
- [ ] Visual verify in Storybook (serve `storybook-static`, Playwright MCP) the four new stories against `2026-05-29-series-plot.html`.
- [ ] `docs/greenfield-component-ledger.md`: flip `MemberList`, `ReadingPanel`, `SeriesPlate` rows → ✅ with paths; add a `SeriesMemberRow` L2 row (✅); add decision rows: heatmap rep-toggle **omitted** (HeatmapChart out-of-scope), phase-strip companion **deferred** (needs L1 renderer aligned to WaterfallChart row geometry), `SeriesMemberRow` **split** from builder `MemberRow` (per line-196 note); bump coverage counts. Reconcile the "MemberList←MemberRow" next-up wording to "MemberList←SeriesMemberRow".
- [ ] Update greenfield memory (`project_greenfield_composite_layer.md` + `MEMORY.md` line): Batch 9 done; gotchas (Swatch-coexistence-gradient-in-ui; heatmap-omitted/phasestrip-deferred; SeriesMemberRow≠builder-MemberRow).
- [ ] Commit ledger + (separately, outside git) memory. `git add docs/greenfield-component-ledger.md`.

---

## Verification (whole batch)

`npm test -- print-ui print-components` green · `npm run lint:design` clean (proves placement-only) · `npx tsc --noEmit -p tsconfig.build.json` clean · `npm run build-storybook` exit 0 · four new `components/*` stories match the series-plot mockup.
