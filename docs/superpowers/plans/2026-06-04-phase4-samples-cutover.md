# Phase-4 Samples Contact-Sheet Cutover — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.
> Commit trailer (every commit): `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
> Branch stays UNMERGED/UNPUSHED. Commit ONLY named files — never `git add -A`. Typecheck: `npx tsc --noEmit -p tsconfig.build.json`. Work from `packages/HimalayaUI/frontend/`.

**Goal:** Swap the legacy `/samples` contact sheet (`src/pages/SamplesPage.tsx`) for a greenfield `src/print/pages/SamplesPage.tsx` assembled from `SheetTable` + `SampleTableRow` + `CullBar` + `KbLegend`, wired to the carried corpus/exposure hooks, then delete the legacy body — third route in the `Loupe → Focus → Samples → Series` cutover.

**Architecture:** Incremental-by-route under the carried shell (`AppRoutes` repoints one route; the SSE/queue/QueryClient mechanisms are untouched). The page is the only stateful unit: it owns a **page-global cull selection** (`Set<exposureId>`), the screened-progress aggregate, the keyboard handlers, and all navigation. The presentational composites get minimal **refactor-on-contact** nav seams (they have none today). The cross-sample batch-reject dispatches per-sample via a new `useSetExposureStatusBatch` hook (feasible because `useQueueMutation` merges per-call input over bound scope).

**Tech Stack:** React 18, TanStack Query v5 (`useQueries`), Zustand, react-router v6, Tailwind (design-guard-enforced), boneyard-js skeletons, Vitest + RTL, Playwright (mocked).

**Provenance discipline (hard rules from `docs/superpowers/specs/2026-06-03-phase4-cutover-strategy-design.md`):**
- NEW = `src/print/**`; CARRIED = `queries.ts`/`api.ts`/`state.ts`/`lib/**`/hooks; OLD = `src/pages/**` + legacy `src/components/**` (deleted this cutover).
- The page is assembled from `src/print` composites + the **mockup** (`docs/redesign-mockups/sample-table.html`) content ONLY. Reuse logic, never presentation. Do not re-create a legacy affordance just because legacy had it.
- **Resolved product decision (this plan):** the **Status cell is the Focus door** (`click → /sample/:id`); the **sample name → loupe** and **thumb double-click → loupe** (mockup); thumb click = cull-select. The mockup omits cross-page nav to Focus (it is a screening mockup); the Focus door on Status is the load-bearing nav the user requires (Focus is reached via the samples page), placed on the semantically-correct cell.

---

## Mockup contract (`docs/redesign-mockups/sample-table.html`, sheet view)

- **Head:** `Kicker tone="accent"` "Contact sheet" · serif `<h1>` (the **beamtime/corpus title**, carried filter logic — NOT the mockup's hardcoded string) · `.sub` paragraph · right-aligned `progress` block = `N / M` numeral + "samples screened" + a fill bar (`ProgressBar`).
- **Table:** `SheetTable` (elevated Card; aligned 5-col header Sample/Exposures/Kept/Tags/Status via the shared `SAMPLE_TABLE_COLS`) slotting one `SampleTableRow` per sample.
- **Footer:** `KbLegend` with exactly these 5 shortcuts (`Shortcut = { keyLabel, description }`):
  - `click` — select a frame
  - `⇧ click` — extend the range
  - `X` — drop the selected frames
  - `double-click` — open the loupe
  - `Esc` — clear
- **CullBar** (floating, page root): `count` = selection size, `show` = size>0, Drop (`X`) / Restore / Clear (`Esc`).
- **Interactions (JS in the mockup):** selection is a single page-global `Set`; thumb click toggles, ⇧click extends **within one sample only**, double-click thumb → loupe at that frame, click sample name → loupe at the representative frame, `X` drops selected, `Esc`/Clear clears.
- **Width:** `sheet-shell` `max-width: 1240px` centered → `<PageFrame width="sheet">` (already `max-w-[1240px]`, no change needed).

## Carried logic (REUSE — already verified live)

| Need | Carried source | Shape |
|---|---|---|
| Whole corpus | `useCorpusSamples()` | `{ data: CorpusSample[], isLoading, isError }` |
| Per-sample exposures (O(1) fan-out) | `useCorpusExposures(filtered)` | `{ byId: Map<number, Exposure[]>, isLoading }` |
| Screened aggregate | `useScreenedProgress(filtered)` | `{ screened, total }` |
| Experiment names | `useExperiments()` | `{ data: Experiment[] }` |
| `?beamtime=` filter | `useSearchParams()` | URL param (numeric → `experiment_id` filter) |
| Cull mutation | NEW `useSetExposureStatusBatch()` (Task 4) | `mutate({ sampleId, exposureId, status })` |
| Screened derivation | `isSampleScreened(sample, exposures?)` (`lib/sample/screened.ts`) | `boolean` |
| Display name | `sampleDisplayName(sample)` (`lib/sample/displayName.ts`) | `string` |
| Exposure → filmstrip VM | `toGalleryExposures(exposures)` (`src/print/pages/loupeAdapters.ts`) | `GalleryExposure[]` |
| `SampleTag[] → Tag[]` | `toLoupeTags(tags)` (`loupeAdapters.ts`) | `Tag[]` |

Types: `CorpusSample` (`src/api.ts`) has `id`, `experiment_id`, `display_name: string|null`, `name`, `tags: SampleTag[]`, `phase?: string|null`. `Exposure` has `id`, `status: "accepted"|"rejected"|...`, `selected`, `image_path`, `image_version`.

---

## Task 1: `ThumbnailGallery` + `Thumbnail` — double-click → loupe (`onActivate`)

**Why:** the gallery exposes only `onSelect` (single click). The mockup's "double-click — open the loupe" needs a second handler. Refactor-on-contact in `src/print/` (guard-exempt for `ui/`; `ThumbnailGallery` is in `components/` so stays placement-only — this change adds only a prop + handler, no appearance).

**Files:**
- Read+Modify: `src/print/ui/Thumbnail.tsx` (the leaf — read it first to mirror its `onClick` wiring)
- Modify: `src/print/components/ThumbnailGallery.tsx`
- Test: `test/print-components/ThumbnailGallery.test.tsx` (add cases; create if absent)

- [ ] **Step 1: Read `src/print/ui/Thumbnail.tsx`** to confirm its prop names (`onClick`, `selected`, `size`, `frameNo`) and root element.

- [ ] **Step 2: Failing test** — add to `test/print-components/ThumbnailGallery.test.tsx`:

```tsx
it("fires onActivate(id) on double-click of a thumb", () => {
  const onActivate = vi.fn();
  const { getByTestId } = render(
    <ThumbnailGallery
      exposures={[{ id: 7, src: null, frameNo: 1 }]}
      onActivate={onActivate}
    />,
  );
  fireEvent.doubleClick(getByTestId("thumbnail-7"));
  expect(onActivate).toHaveBeenCalledWith(7);
});
```
(Confirm the thumb's `data-testid` from Step 1; adjust the selector to the real value.)

- [ ] **Step 3: Run → fail.** `npm test -- ThumbnailGallery`

- [ ] **Step 4: Implement.** In `Thumbnail.tsx` add an optional `onDoubleClick?: () => void` and wire it to the root element's `onDoubleClick` (mirror the existing `onClick` conditional-spread pattern). In `ThumbnailGallery.tsx`:
  - Add to `ThumbnailGalleryProps`: `onActivate?: (id: number) => void;`
  - Destructure `onActivate`, and on each `Thumbnail` add: `{...(onActivate ? { onDoubleClick: () => onActivate(exposure.id) } : {})}`

- [ ] **Step 5: Run → pass** + `npm run lint:design` + `npx tsc --noEmit -p tsconfig.build.json`.

- [ ] **Step 6: Commit** `git add src/print/ui/Thumbnail.tsx src/print/components/ThumbnailGallery.tsx test/print-components/ThumbnailGallery.test.tsx`

---

## Task 2: `SpecCell` — sample name → loupe (`onOpenLoupe`)

**Why:** the mockup makes the sample name the loupe trigger. `SpecCell` renders the name as a static `<span>`. Add an optional callback that promotes the name to a keyboard-accessible button; without it the name stays static (Storybook/other consumers unaffected).

**Files:**
- Modify: `src/print/components/SpecCell.tsx`
- Test: `test/print-components/SpecCell.test.tsx` (add cases; create if absent)

- [ ] **Step 1: Failing test:**

```tsx
it("renders the name as a button that fires onOpenLoupe when provided", () => {
  const onOpenLoupe = vi.fn();
  render(<SpecCell name="POPC" sampleId="#42" onOpenLoupe={onOpenLoupe} />);
  fireEvent.click(screen.getByRole("button", { name: /POPC/ }));
  expect(onOpenLoupe).toHaveBeenCalled();
});
it("renders the name as static text when onOpenLoupe is absent", () => {
  render(<SpecCell name="POPC" sampleId="#42" />);
  expect(screen.queryByRole("button", { name: /POPC/ })).toBeNull();
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement.** Add `onOpenLoupe?: () => void;` to `SpecCellProps`. Replace the name `<span data-role="spec-name" …>` with a conditional: when `onOpenLoupe` is set, render
```tsx
<button
  type="button"
  data-role="spec-name"
  onClick={onOpenLoupe}
  className="text-body font-semibold block leading-tight line-clamp-2 text-left hover:text-print-accent"
>
  {name}
</button>
```
else keep the existing `<span>`. (`text-print-accent`/`text-body` are named token/scale utilities — design-guard clean; no `text-[…]` literal.)

- [ ] **Step 4: Run → pass** + `lint:design` + `tsc`.

- [ ] **Step 5: Commit** `git add src/print/components/SpecCell.tsx test/print-components/SpecCell.test.tsx`

---

## Task 3: `SampleTableRow` — navigation seams

**Why:** the row is presentational; the page needs three nav hooks. Forward two to children, own the Status-cell door wrapper here (StatusCell stays pure).

**Files:**
- Modify: `src/print/components/SampleTableRow.tsx`
- Test: `test/print-components/SampleTableRow.test.tsx` (add cases)

- [ ] **Step 1: Failing tests:**

```tsx
it("forwards onOpenLoupe to the SpecCell name button", () => {
  const onOpenLoupe = vi.fn();
  render(<SampleTableRow {...baseProps} onOpenLoupe={onOpenLoupe} />);
  fireEvent.click(screen.getByRole("button", { name: new RegExp(baseProps.name) }));
  expect(onOpenLoupe).toHaveBeenCalled();
});
it("forwards onActivateExposure to the gallery (double-click)", () => {
  const onActivateExposure = vi.fn();
  render(<SampleTableRow {...baseProps} onActivateExposure={onActivateExposure} />);
  fireEvent.doubleClick(screen.getByTestId(`thumbnail-${baseProps.exposures[0].id}`));
  expect(onActivateExposure).toHaveBeenCalledWith(baseProps.exposures[0].id);
});
it("makes the status cell a Focus door when onOpenFocus is set", () => {
  const onOpenFocus = vi.fn();
  render(<SampleTableRow {...baseProps} phase={null} onOpenFocus={onOpenFocus} />);
  fireEvent.click(screen.getByRole("button", { name: /index/i }));
  expect(onOpenFocus).toHaveBeenCalled();
});
```
(`baseProps` = a minimal valid prop set: `name`, `sampleId`, `exposures: [{id:1,src:null,frameNo:1}]`, `kept:1`, `total:1`, `tags:[]`.)

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement.** Add to `SampleTableRowProps`: `onOpenLoupe?: () => void;`, `onOpenFocus?: () => void;`, `onActivateExposure?: (id: number) => void;`.
  - Pass `onOpenLoupe` to `<SpecCell … {...(onOpenLoupe ? { onOpenLoupe } : {})} />`.
  - Pass `{...(onActivateExposure ? { onActivate: onActivateExposure } : {})}` to `<ThumbnailGallery>`.
  - Wrap the **5th cell** (Status). Replace `<div className={CELL}><StatusCell … /></div>` with: when `onOpenFocus` is set, make the cell a button:
```tsx
<div className={CELL}>
  {onOpenFocus ? (
    <button
      type="button"
      onClick={onOpenFocus}
      aria-label={phase ? `Open indexing for ${name} (${phase})` : `Index ${name}`}
      className="flex items-center rounded-md px-2 -mx-2 transition-colors hover:bg-plate/60"
    >
      <StatusCell {...(phase !== undefined ? { phase } : {})} />
    </button>
  ) : (
    <StatusCell {...(phase !== undefined ? { phase } : {})} />
  )}
</div>
```
(`hover:bg-plate/60`, `rounded-md` are named token/radius utilities — guard clean.)

- [ ] **Step 4: Run → pass** + `lint:design` + `tsc`.

- [ ] **Step 5: Commit** `git add src/print/components/SampleTableRow.tsx test/print-components/SampleTableRow.test.tsx`

---

## Task 4: `useSetExposureStatusBatch` — cross-sample cull dispatch

**Why:** page-global selection can span samples, but `useSetExposureStatus(sampleId)` binds one sample. `useQueueMutation.mutate` builds `payload = { kind, clientOpId, ...scope, ...input }` (input spread last — verified in `src/lib/queue/useQueueMutation.ts:124-130`), so `sampleId` passed per-call **overrides** the bound scope and lands in `p.sampleId`. One hook can therefore dispatch to any sample; each `mutate` mints its own `client_op_id` and patches `queryKeys.exposures(sampleId)`.

**Files:**
- Modify: `src/queries.ts` (add hook next to `useSetExposureStatus`, ~line 619)
- Test: `test/queries-batch-exposure-status.test.tsx` (new) — or extend an existing queries test if the suite has one for mutators.

- [ ] **Step 1: Failing test** — assert that a single hook can target two different `sampleId`s. Mock `useQueueMutation` to capture the payloads:

```tsx
// Mock useQueueMutation to record (scope, input) pairs; assert that calling
// mutate({sampleId:2,...}) yields payload.sampleId === 2 even though the hook
// was constructed with scope.sampleId === 0.
```
(Mirror the queue-mutation test style already in `test/`; if no precedent exists, assert at minimum that the returned `mutate` forwards `{sampleId, exposureId, status}` to the inner mutation untouched.)

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement** in `src/queries.ts`:

```ts
/** Batch/cross-sample exposure-status setter for the contact-sheet cull.
 *  Unlike useSetExposureStatus (one bound sample), the sampleId rides in the
 *  per-call input and overrides the placeholder scope (useQueueMutation spreads
 *  input over scope), so one hook dispatches reject/restore to any sample. */
export function useSetExposureStatusBatch() {
  const username = useAppState((s) => s.username);
  const inner = useQueueMutation(setExposureStatusMutator, {
    sampleId: 0, username, clientId: CLIENT_ID,
  });
  return {
    ...inner,
    mutate: (v: { sampleId: number; exposureId: number; status: "accepted" | "rejected" | null }) =>
      // sampleId rides in the input position → payload.sampleId === v.sampleId.
      inner.mutate(v as unknown as { exposureId: number; status: "accepted" | "rejected" | null }),
  };
}
```
(The `as unknown as` cast is required because `SetExposureStatusInput` does not declare `sampleId`; the runtime payload merge is what carries it. Document this inline.)

- [ ] **Step 4: Run → pass** + `tsc`.

- [ ] **Step 5: Commit** `git add src/queries.ts test/queries-batch-exposure-status.test.tsx`

---

## Task 5: `samplesAdapters.ts` — pure row-model derivation

**Why:** keep the page thin; one tested pure function maps `CorpusSample` + its `Exposure[]` → `SampleTableRow` props.

**Files:**
- Create: `src/print/pages/samplesAdapters.ts`
- Test: `test/print-pages/samplesAdapters.test.ts`

- [ ] **Step 1: Failing test** for `toSampleRowModel`:

```ts
import { toSampleRowModel } from "../../src/print/pages/samplesAdapters";
it("derives kept/total/dropped/screened/tags/phase from a sample + exposures", () => {
  const sample = { id: 9, experiment_id: 1, name: "JC009", display_name: "JC009 — LL37",
    tags: [{ id: 1, key: "LL37", value: "" }], phase: null } as any;
  const exposures = [
    { id: 101, status: "accepted", selected: true, image_path: null, image_version: "" },
    { id: 102, status: "rejected", selected: false, image_path: null, image_version: "" },
  ] as any;
  const m = toSampleRowModel(sample, exposures);
  expect(m.sampleId).toBe("#9");
  expect(m.kept).toBe(1);
  expect(m.total).toBe(2);
  expect(m.dropped).toBe(1);
  expect(m.exposures).toHaveLength(2);
  expect(m.tags).toEqual([{ key: "LL37" }]);
  expect(m.phase).toBeNull();
});
it("treats undefined exposures as not-yet-loaded (empty derivation)", () => {
  const sample = { id: 9, tags: [], phase: undefined } as any;
  const m = toSampleRowModel(sample, undefined);
  expect(m.total).toBe(0);
  expect(m.kept).toBe(0);
});
```

- [ ] **Step 2: Run → fail.**

- [ ] **Step 3: Implement:**

```ts
import type { CorpusSample, Exposure } from "../../api";
import type { Tag } from "../ui";
import type { GalleryExposure } from "../components/ThumbnailGallery";
import { sampleDisplayName } from "../../lib/sample/displayName";
import { isSampleScreened } from "../../lib/sample/screened";
import { toGalleryExposures, toLoupeTags } from "./loupeAdapters";

export interface SampleRowModel {
  name: string;
  sampleId: string;
  screened: boolean;
  exposures: GalleryExposure[];
  kept: number;
  total: number;
  dropped: number;
  tags: Tag[];
  phase: string | null | undefined;
}

/** CorpusSample + its (possibly unloaded) exposures → SampleTableRow props.
 *  undefined exposures = "not yet fetched" → empty derivation (the page shows
 *  the boneyard skeleton meanwhile). */
export function toSampleRowModel(
  sample: CorpusSample,
  exposures: Exposure[] | undefined,
): SampleRowModel {
  const exps = exposures ?? [];
  const total = exps.length;
  const kept = exps.filter((e) => e.status !== "rejected").length;
  return {
    name: sampleDisplayName(sample),
    sampleId: `#${sample.id}`,
    screened: isSampleScreened(sample, exposures),
    exposures: toGalleryExposures(exps),
    kept,
    total,
    dropped: total - kept,
    tags: toLoupeTags(sample.tags),
    phase: sample.phase,
  };
}
```
(Verify `sampleDisplayName`/`isSampleScreened` import paths against the live files before finalizing.)

- [ ] **Step 4: Run → pass** + `tsc`.

- [ ] **Step 5: Commit** `git add src/print/pages/samplesAdapters.ts test/print-pages/samplesAdapters.test.ts`

---

## Task 6: Assemble `src/print/pages/SamplesPage.tsx`

**Why:** the cutover body. Composes the head + `SheetTable` + rows + `CullBar` + `KbLegend`, owns page-global selection + keyboard + nav + the `?beamtime` filter, wires the carried hooks. Mirrors `LoupePage.tsx`'s structure (PageFrame, boneyard Skeleton, keyboard `useEffect`, `navigate`).

**Files:**
- Create: `src/print/pages/SamplesPage.tsx`
- Test: `test/print-pages/SamplesPage.test.tsx` (mirror `test/print-pages/FocusPage.test.tsx` mock-plane style)

**Page composition (full structure):**
- `<PageFrame width="sheet">` wrapper, `data-testid="samples-page"`.
- **Filter + title:** read `?beamtime` (numeric → `experiment_id`); `filtered = beamtime ? samples.filter(s.experiment_id===beamtime) : samples`; title = `beamtime === undefined ? "The corpus" : (experiments.find(...)?.name ?? "experiment N")` (carried logic from legacy lines 52-76).
- **Head:** `Kicker tone="accent"` "Contact sheet" + serif `<h1 className="text-display text-ink">{title}</h1>` + `.sub` paragraph (mockup copy) + right block: `{screened} / {total}` numeral (`text-headline-lg`) + `Kicker tone="faint"` "samples screened" + `<ProgressBar value={screened} total={total} label="samples screened" className="w-40 mt-1" />`.
- **Selection state:** `const [selected, setSelected] = useState<Set<number>>(new Set())`. A `useRef` anchor `{ sampleId, exposureId } | null` for ⇧-range (range only extends within the same sample — guard on `sampleId`).
- **exposureId → sampleId lookup:** build a `Map<number, number>` from `corpusExposures.byId` (each entry's exposures → their `sample.id`). Needed by the batch cull.
- **Table:** `<SheetTable empty={<div className="p-10 text-center text-ink-faint">No samples {beamtime===undefined?"in the corpus":"in this beamtime"} yet</div>}>` slotting, for each `filtered` sample:
```tsx
const m = toSampleRowModel(s, corpusExposures.byId.get(s.id));
<SampleTableRow
  key={s.id}
  name={m.name} sampleId={m.sampleId} screened={m.screened}
  exposures={m.exposures} kept={m.kept} total={m.total} dropped={m.dropped}
  tags={m.tags} {...(m.phase !== undefined ? { phase: m.phase } : {})}
  selectedExposureIds={selected}
  onSelectExposure={(id) => toggleSelect(s.id, id)}        // click/⇧click = cull-select (shift read from shiftRef)
  onActivateExposure={(id) => navigate(loupeHref(s.id))}   // double-click thumb → loupe
  onOpenLoupe={() => navigate(loupeHref(s.id))}            // name → loupe
  onOpenFocus={() => navigate(`/sample/${s.id}`)}          // status → focus (DECISION)
/>
```
  where `loupeHref(id)` preserves `?beamtime` like `LoupePage.goBack` does.
- **⇧-range (page-level shift tracking, no composite change):** the gallery's `onSelect` is `(id)=>void` and doesn't carry `shiftKey`, so the page tracks shift state itself: a `shiftRef = useRef(false)` updated by window `keydown`/`keyup` listeners (`e.key === "Shift"`). `toggleSelect(sampleId, id)` then branches: with `shiftRef.current === true` **and** the anchor's `sampleId` matches, select the contiguous span between the anchor's exposure and `id` **within that sample's exposure order** (from `corpusExposures.byId.get(sampleId)`); otherwise toggle the single id and re-anchor `{ sampleId, exposureId: id }`. Range never crosses samples (anchor-sample guard), matching the mockup's `lastClick.si === si` rule.
- **CullBar:** `<CullBar count={selected.size} show={selected.size>0} onReject={()=>batchSet("rejected")} onRestore={()=>batchSet(null)} onClear={()=>setSelected(new Set())} />` where `batchSet(status)` iterates `selected`, looks up each `sampleId`, and calls `batch.mutate({ sampleId, exposureId: id, status })` (skipping ones already in that state), then clears.
- **Keyboard:** a `useEffect` (mirroring `LoupePage`): `X`/`x` → `batchSet("rejected")`; `Escape` → clear; suppressed while typing in `INPUT`/`TEXTAREA`; only meaningful while `selected.size>0`.
- **Footer legend:** `<KbLegend shortcuts={[{keyLabel:"click",description:"select a frame"},{keyLabel:"⇧ click",description:"extend the range"},{keyLabel:"X",description:"drop the selected frames"},{keyLabel:"double-click",description:"open the loupe"},{keyLabel:"Esc",description:"clear"}]} />` (all 5 mockup hints — ⇧-range is wired via shift tracking above).
- **Loading/error:** wrap the table in the boneyard `<Skeleton name="contact-sheet" loading={corpusQuery.isLoading} …>` with a greenfield `CONTACT_SHEET_FIXTURE` (skeleton rows shaped to `SheetTable`); render an error block on `corpusQuery.isError`.

- [ ] **Step 1: Write `test/print-pages/SamplesPage.test.tsx`** (mock-plane like FocusPage.test). Cases:
  - renders head (title from beamtime), `sheet-table`, one `sample-table-row` per sample, `kb-legend`.
  - clicking a thumb toggles selection → `cull-bar` `data-show="true"` and count reflects size.
  - ⇧-click extends the contiguous range within one sample (and never crosses into another sample's frames).
  - CullBar Drop calls the batch mutate with `{sampleId, exposureId, status:"rejected"}` for each selected.
  - pressing `X` with a selection drops; `Escape` clears.
  - clicking a row's Status door navigates to `/sample/:id`; clicking the name navigates to `/samples/loupe/:id`.
  - `?beamtime=` filters rows + sets the title.
  - empty + error states.
- [ ] **Step 2: Run → fail.**
- [ ] **Step 3: Implement `SamplesPage.tsx`** per the composition above. Imports NEW + CARRIED only.
- [ ] **Step 4: Run → pass** + `lint:design` + `tsc`.
- [ ] **Step 5: Storybook check (optional):** the existing `ContactSheetAssembly` story already covers the visual; no new story required for the page (pages aren't storied).
- [ ] **Step 6: Commit** `git add src/print/pages/SamplesPage.tsx test/print-pages/SamplesPage.test.tsx`

---

## Task 7: Repoint route, delete legacy, recapture boneyard

**Files:**
- Modify: `src/components/AppRoutes.tsx` (line 8 import + line 88 element)
- Delete: `src/pages/SamplesPage.tsx`, `src/components/ContactSheetRow.tsx` (+ any now-orphaned legacy: `RejectXMark`, `SampleStatusChip`, legacy `src/components/CullBar.tsx` — **grep before each delete**)
- Modify/replace: `test/contact-sheet.test.tsx` (legacy presentation test — delete; coverage moves to `SamplesPage.test.tsx`), `test/CorpusShell.test.tsx` (drops/updates the legacy `SamplesPage` import)
- Recapture: `src/bones/contact-sheet.bones.json` + `src/bones/registry.ts`

- [ ] **Step 1: Repoint.** In `AppRoutes.tsx` change `import { SamplesPage } from "../pages/SamplesPage";` → `from "../print/pages/SamplesPage";`. App still runs (mixed).
- [ ] **Step 2: Orphan analysis.** For each legacy candidate run `grep -rln "<name>" src` and delete ONLY if every remaining importer is itself being deleted. Known importers to resolve: `ContactSheetRow` is referenced by `src/components/CullBar.tsx`, `src/components/NoUsableExposureNotice.tsx`, `test/contact-sheet.test.tsx` — confirm whether those are legacy-only and deletable, or live (if live, leave them and only delete the page). **Do not** delete anything still imported by surviving code.
- [ ] **Step 3: Delete** the confirmed-orphan legacy files (`git rm`).
- [ ] **Step 4: Fix test fallout.** Delete `test/contact-sheet.test.tsx` (legacy presentation). Update `test/CorpusShell.test.tsx` so it no longer imports the legacy page (it should exercise routing through the greenfield page or a stub). Run `npm test` → green.
- [ ] **Step 5: Recapture the boneyard skeleton** for the greenfield page per `feedback_boneyard_capture_recipe` (backend-proxied dev server + single-URL boneyard-js build to a temp dir + additive `registry.ts` edit). Confirm `contact-sheet.bones.json` reflects the new `SheetTable` shape, not the legacy grid. Stage only `src/bones/contact-sheet.bones.json` + `src/bones/registry.ts`.
- [ ] **Step 6: Commit** the swap: `git add src/components/AppRoutes.tsx src/bones/contact-sheet.bones.json src/bones/registry.ts test/CorpusShell.test.tsx` + the `git rm`'d paths.

---

## Task 8: Gate

- [ ] **Visual fidelity:** serve `sample-table.html` over HTTP + Vite against a **dev-DB copy** backend; screenshot-compare the live `/samples` against the mockup sheet view (head, table alignment, CullBar, legend, 1240 width). Verify the three nav doors by clicking (name→loupe, status→focus, dblclick→loupe). Render with the external `/Volumes/data` mount so thumbnails are real.
- [ ] **Unit:** `npm test` green (the new print-components/print-pages suites + no legacy regressions).
- [ ] **E2E (mocked):** `npm run e2e` — fix any contact-sheet spec that asserted legacy selectors (rewrite to greenfield `data-testid`s, per the assignment-model rule, not legacy groups).
- [ ] **Design guard:** `npm run lint:design` clean.
- [ ] **Types:** `npx tsc --noEmit -p tsconfig.build.json` clean.
- [ ] **Build:** `npm run build` exit 0.
- [ ] **No deferrals:** all 5 mockup interactions (click-select, ⇧-range, X-drop, double-click-loupe, Esc-clear) are wired; the legend matches the live affordances 1:1.

---

## Self-review (writing-plans)

- **Spec coverage:** every per-route recipe step (build gaps → assemble → wire → repoint → delete → gate) maps to Tasks 1-8. The resolved Focus-door decision is encoded in Task 3 + Task 6.
- **Type consistency:** `onActivate`/`onActivateExposure`/`onOpenLoupe`/`onOpenFocus` used identically across Tasks 1/2/3/6; `useSetExposureStatusBatch.mutate({sampleId,exposureId,status})` signature matches its call sites in Task 6; `toSampleRowModel` return shape matches the `SampleTableRow` props it feeds.
- **Verify-at-execution flags:** `Thumbnail.tsx` thumb `data-testid` + `onClick` pattern (Task 1 Step 1), `sampleDisplayName`/`isSampleScreened` import paths (Task 5), and the legacy orphan set (Task 7 Step 2) are the three spots to confirm against live source before editing — each is called out in-step.
