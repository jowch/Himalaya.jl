# Greenfield Batch 13 — Series-builder rail slice (`Field` + `RailBack` + `Dock` + `BuilderRail`)

> Execute with superpowers:subagent-driven-development (fresh implementer per task + spec/quality review).
> Commit trailer (exact last line): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Worktree frontend dir: `packages/HimalayaUI/frontend`. Typecheck: `npx tsc --noEmit -p tsconfig.build.json`.
> Commit ONLY named files. NEVER `git add -A`/`.`. NEVER stage `src/bones/*` or anything under `docs/superpowers/plans/`.

## Context

Branch `worktree-greenfield-ui-rebuild` (NOT merged). Batch 12 (scoping) done. This is the LAST Layer-3 panel slice — the **series-builder editing rail** ("Compose" margin) + its two floating overlays, derived from `docs/redesign-mockups/series-builder.html`. The figure plate (`SeriesPlate`←`WaterfallChart`) already exists (Batch 9); this slice builds the rail that drives it.

**Tasks (dependency-ordered):**
1. `Field` (new Layer-0 primitive, `ui/`) — the `.field` ordering-variable control (value + ▾ chevron).
2. `RailBack` + `Dock` (`components/`) — the two floating overlays (`.railback`, `.dock`).
3. `BuilderRail` (`components/`) + `SeriesBuilderAssembly` story — the `.rail` (needs Field + the overlays for the assembly).

**LOAD-BEARING OMISSIONS ("controls don't lie" — consistent with Batch 9 `SeriesPlate`):** the mockup rail has a **Representation** section (Waterfall/Heatmap `.seg`) and a **Track reflections** toggle (`.check`). Both drive renderers that are out-of-scope / deferred (`HeatmapChart` ⛔ out-of-scope; `TrackingLine` ⏸ deferred — see the ledger). **OMIT both** from `BuilderRail`. The Display section keeps only the offset slider + the log/linear scale toggle.

**NO MOCKUP → NOT IN THIS SLICE:** `ExportSheet` and `ConflictModal` have no file in `docs/redesign-mockups/`. Per the derive-from-mockup rule they are NOT buildable now; leave them ⬜ (they'd need their own mockup/brainstorm).

**Layer-4 deferral:** the topbar (Duplicate/Export buttons), the `.main` grid + rail-collapse *mechanics*, the dropdown MENU behind `Field`, the offset↔dock sync, and durable wiring are page assembly. `BuilderRail`/`RailBack`/`Dock` are PRESENTATIONAL (every interaction datum is a prop with a consumer handler — the Batch 7 contract). The `SeriesBuilderAssembly` story simulates the page (owns offset/scale/collapsed state, the offset↔dock sync, the trace-list data).

### Verified reuse APIs (checked against live source 2026-06-03 — do NOT re-derive)

- `RailSection` (`./RailSection`) — `{ label: ReactNode; count?; note?; children?; className? }`. Renders a `Kicker tone="faint"` label, a children slot, and an optional `text-caption` note. `data-testid="rail-section"`. The rail sections compose this.
- `MemberRow` (`./MemberRow`) — `{ name; sub?: ReactNode; phase: string; coexistWith?: string[]; className? }`. The `.trow`: `GripHandle`+`Swatch`+name/sub+`PhaseChip`, carries `group`, `data-testid="member-row"`. Presentational (no onClick/drag). The builder trace list is `MemberRow`×N slotted into a `flex flex-col gap-0.5` wrapper — NO separate "builder MemberList" component (it'd be a thin flex wrapper; render inline in `BuilderRail`). Drag-reorder is page-deferred (the mockup itself implements none).
- `AutoGroup` (`./AutoGroup`) — `{ variant?: "summary"|"compose"; title?; children; actions?: {label;onClick?;muted?}[]; className? }`. Use `variant="compose"` (bg-plate) + `title="Auto-grouped"` + body children + `actions=[{label:"Confirm series"},{label:"Adjust",muted:true}]`. It renders the accent star + the link-style actions itself.
- `Slider` (`../ui`) — `{ value; min; max; step?; onChange:(v:number)=>void; label?; valueDisplay?; ariaLabel?; className? }`. The offset row: `<Slider label="Trace offset" valueDisplay={`${offset.toFixed(2)}×`} value min={0.4} max={1.4} step={0.05} onChange={onOffsetChange} />`.
- `SegmentedControl<T>` (`../ui`) — `{ options: {value;label;...}[]; value: T; onChange:(t:T)=>void; "aria-label": string; variant?; size?; stretch?; className? }`. The scale toggle: `options=[{value:"log",label:"log q"},{value:"lin",label:"linear q"}]`, `aria-label="q scale"`, `stretch`. Active = `bg-ink text-paper` (built-in).
- `IconButton` (`../ui`) — `{ label: string; tone?: "ghost"|"accent"|"danger"; ...buttonprops }` + an icon child. The rail-head collapse `railtog` (a `›` glyph). `label` is the a11y name (e.g. "Collapse rail").
- `Button` (`../ui`) — `variant?: "solid"|"accent"|"ghost"|"danger"|"outline"|"ghostInverse"`. The rail-foot "+ Add sample" / "Copy as PNG" = `variant="outline"` (the mockup `.btn`: bordered `bg-plate`). NO `size` prop.
- `Kicker` (`../ui`) — `tone="faint"` = the rail-head "Compose" label role + the rail-h section labels (RailSection already uses it).
- `SeriesPlate` (`./SeriesPlate`) — `{ kicker?; title; subtitle?; rows: WaterfallRow[]; scale: "log"|"lin"; onScaleChange; hoveredKey?; onHoverRow?; hoveredQ?; onHoverQ?; legendPhases?; footHint?; footNote?; className? }`. The assembly's FIGURE side — READ `SeriesPlate.stories.tsx` to mirror how it builds `rows: WaterfallRow[]` (from `realMembers`).
- Fixtures — `import { realMembers } from "../fixtures/realSeriesMembers"` (exported name is `realMembers`, NOT `realSeriesMembers`). For the rail trace list, the assembly can define a small inline `{name, sub, phase, coexistWith?}[]` (6 members mirroring the mockup SERIES) — `MemberRow` needs only those.
- House `cx` — copy into each file (no shared export).
- `exactOptionalPropertyTypes: true` — conditional-spread all optional props.

### Design-guard contract
`ui/**` exempt (Field may author appearance). `components/**` NOT — primitives + tokens + layout only. Floating overlays: `fixed`, `right-0`/`bottom-6`, `z-50`, `-translate-y-1/2` are positioning (allowed). `[writing-mode:vertical-rl]` is an arbitrary LAYOUT property (not `text-[`/`rounded-[`/a color) — allowed. `rounded-md`/`rounded-l-lg` allowed; the mockup's 8/11px radii snap to the system step. `shadow-lg` allowed. `npm run lint:design` MUST stay clean. Tests assert via data/role/text, never class strings.

---

## Task 1: `Field` primitive (`ui/`)

Mockup `.field` (CSS ~260–267): a bordered (`border-hair-strong`, `bg-plate`, rounded) clickable row, `justify-between`: a value (`text-ink`, ~12.5px/600) on the left, a `▾` chevron (`text-ink-faint`) on the right. Hover → no strong change (keep simple). It is a `<button type="button">`; clicking opens a menu (page-owned). `data-testid="field"`.

Props: `{ value: string; onClick?: () => void; placeholder?: string; className? }`. (If `value` empty + `placeholder` given, show placeholder in `text-ink-faint`.)

Files: `src/print/ui/Field.tsx` + `.stories.tsx`; modify `src/print/ui/index.ts` (`export { Field } from "./Field"; export type { FieldProps } from "./Field";`); test `test/print-ui/Field.test.tsx`.

TDD: test renders the value, is a `<button>`, fires onClick (3 cases incl. placeholder fallback). Implement with a `<button type="button" data-testid="field" className="w-full flex items-center justify-between border border-hair-strong bg-plate rounded px-3 py-2 text-meta font-semibold text-ink cursor-pointer">` + value span + `<span className="text-ink-faint" aria-hidden>▾</span>`. Explicit `: JSX.Element` return type. Story `ui/Field` (Default + Placeholder). Gate (vitest + lint:design + tsc) + commit the 4 files:
`feat(print/ui): Field — value + chevron dropdown-affordance control`

## Task 2: `RailBack` + `Dock` (`components/`)

**`RailBack`** — mockup `.railback` (CSS ~218–233): a `fixed right-0 top-1/2 -translate-y-1/2 z-50` tab, `bg-plate border border-hair-strong border-r-0 rounded-l-lg px-2 py-3.5 shadow-lg`, holding a `‹` glyph + a vertical-text "Compose" label (`[writing-mode:vertical-rl] text-kicker`/`Kicker`-like uppercase tracked). Props `{ label?: string; onClick?: () => void; className? }` (label default "Compose"). `data-testid="rail-back"`, `type="button"`. It's the affordance to restore a collapsed rail; the consumer decides when to render it (the assembly only mounts it when collapsed).

**`Dock`** — mockup `.dock` (CSS ~342–359): a `fixed right-6 bottom-6 z-50 flex items-center gap-3.5 bg-plate border border-hair-strong rounded-xl px-4 py-3 shadow-lg` floating pill: an uppercase `text-label`/kicker "Offset" label + a `Slider` (`ariaLabel="Trace offset"`, no visible label, `className="w-32"`) + a mono value (`text-data font-bold`). Props `{ offset: number; min?: number; max?: number; step?: number; onOffsetChange:(v:number)=>void; className? }`. `data-testid="dock"`. Presentational.

Files: `src/print/components/RailBack.tsx` + `.stories.tsx`; `src/print/components/Dock.tsx` + `.stories.tsx`; tests `test/print-components/RailBack.test.tsx`, `test/print-components/Dock.test.tsx`.

TDD each (RailBack: renders label + fires onClick; Dock: renders the offset value, renders a slider, fires onOffsetChange). Stories render them over a paper backdrop (note: `fixed` positioning means the story shows them pinned to the viewport — that's expected). Gate + commit (commit the 6 files together OR two commits — prefer ONE commit for the overlay pair):
`feat(print/components): RailBack + Dock — series-builder floating overlays`

## Task 3: `BuilderRail` (`components/`) + `SeriesBuilderAssembly` story

Mockup `.rail` (markup ~416–491). Presentational (NO `useState`; offset/scale/collapsed + handlers are props; the trace rows are a slot). Sections top→bottom:
1. rail-head: a flex-between row — `Kicker tone="faint"` "Compose" + an `IconButton label="Collapse rail" tone="ghost"` holding `›` (`onClick={onCollapse}`).
2. `<AutoGroup variant="compose" title="Auto-grouped" actions=[{label:"Confirm series",onClick:onConfirm},{label:"Adjust",muted:true,onClick:onAdjust}]>` + the body ReactNode (`grouping` prop, with `<strong>` emphasis).
3. `<RailSection label="Ordering variable" note={orderNote}>` → `<Field value={orderedBy} onClick={onChangeOrder} />`.
4. `<RailSection label="Display">` → `<Slider label="Trace offset" valueDisplay={`${offset.toFixed(2)}×`} value={offset} min={0.4} max={1.4} step={0.05} onChange={onOffsetChange} />` + `<SegmentedControl aria-label="q scale" stretch options={[{value:"log",label:"log q"},{value:"lin",label:"linear q"}]} value={scale} onChange={onScaleChange} />`. (NO Representation section, NO Track toggle — see OMISSIONS above.)
5. `<RailSection label="Traces — drag to reorder">` → `<div className="flex flex-col gap-0.5">{traces}</div>` (the `traces` ReactNode slot = `MemberRow`×N from the page).
6. rail-foot: `<div className="flex gap-2">` two `<Button variant="outline" className="flex-1">` — "+ Add sample" (`onAddSample`) + "Copy as PNG" (`onCopyPng`).
7. a closing `text-caption text-ink-faint` hint: "The plate above is the figure as it will export. What you compose is what you publish."

Root: `<aside data-testid="builder-rail" className="flex flex-col gap-5 bg-paper-sunk border-l border-hair px-5 pt-4 pb-7 overflow-y-auto">` (the rail surface). Width is page-owned (the assembly sets it).

Props (presentational):
```ts
import type { ReactNode } from "react";
type Scale = "log" | "lin";
export interface BuilderRailProps {
  grouping: ReactNode;
  onConfirm?: () => void; onAdjust?: () => void;
  orderedBy: string; orderNote?: ReactNode; onChangeOrder?: () => void;
  offset: number; onOffsetChange: (v: number) => void;
  scale: Scale; onScaleChange: (s: Scale) => void;
  traces: ReactNode;                 // MemberRow×N (children-slotting)
  onAddSample?: () => void; onCopyPng?: () => void; onCollapse?: () => void;
  className?: string;                // PLACEMENT-ONLY
}
```

Test (`test/print-components/BuilderRail.test.tsx`): renders the "Compose" head + the AutoGroup title + the Field value + a slider + the scale segmented control + the traces slot (sentinel `<div data-testid="traces-slot" />`) + both rail-foot buttons; fires `onOffsetChange`/`onScaleChange`/`onAddSample`/`onCollapse`; ASSERTS NO "Heatmap"/"Track reflections" text is present (the load-bearing omission). Implement per the section list, conditional-spread optionals, copy `cx`, explicit return type. 

`SeriesBuilderAssembly.stories.tsx` (read `ContactSheetAssembly.stories.tsx` + `SeriesPlate.stories.tsx` first): the page-sim. Owns `useState` for `offset`/`scale`/`collapsed`. Layout = a `grid` with the figure (`SeriesPlate`, built from `realMembers` like SeriesPlate.stories) on the left and `BuilderRail` on the right (`w-[336px]`), collapsing to full-bleed when `collapsed` (rail hidden, `RailBack` shown). Mount `Dock` (offset, synced to the same state) only when `collapsed`. The rail's `scale`/`offset` and the `SeriesPlate` `scale` share the same state (composing the figure↔rail link the page owns). Trace list = a small inline 6-member array → `MemberRow`×N passed as `traces`. Title `"components/SeriesBuilderAssembly"`, one `Page` story.

Gate (vitest BuilderRail + lint:design + tsc) + commit the 3 files:
`feat(print/components): BuilderRail + SeriesBuilderAssembly — series-builder editing rail`

---

## Batch verification (after all 3 tasks)
```
npx vitest run test/print-ui/Field.test.tsx test/print-components/RailBack.test.tsx \
  test/print-components/Dock.test.tsx test/print-components/BuilderRail.test.tsx
npm run lint:design ; npx tsc --noEmit -p tsconfig.build.json ; npm run build-storybook
```
Visually verify `components/SeriesBuilderAssembly` (both rail-open and collapsed/full-bleed) vs `series-builder.html`. Update `docs/greenfield-component-ledger.md` (Batch 13 row + counts; `BuilderRail`/`Field` ⬜→✅; note the Representation/Track omissions + ExportSheet/ConflictModal still ⬜-no-mockup). Ledger committed; plan file NOT.

**After Batch 13: the Layer-3 panel layer is COMPLETE — every page (Focus/Series/Contact/Loupe/Scoping/Builder) is assemblable. The remaining frontier is Layer-4 pages (router/SSE/queue wiring) + the no-mockup modals.**
