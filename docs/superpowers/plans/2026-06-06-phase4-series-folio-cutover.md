# Phase-4 Plan 4b — Series-folio cutover (`/series` → greenfield)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the legacy Series-folio page body (`/series`) with a greenfield `src/print/pages/SeriesFolioPage.tsx` assembled from the built `FolioHeader`/`Gallery`/`SeriesCard` composites + the d3 waterfall engine, wired to live data via `useSeriesList`/`useSeries`/`useSeriesTraces` (the last shipped in Plan 4a), then delete the OLD.

**Architecture:** The page is the only stateful unit — it owns the search/filter/sort controls (reusing the CARRIED pure `filterSort` from `lib/series/folioFilter.ts`) and maps the filtered `SeriesSummary[]` to a per-card `FolioCard` container that fetches that series' detail + traces and feeds the **presentational** `SeriesCard`. A mini-waterfall row figure is built with the greenfield `toWaterfallRows(members, tracesById)` (real curves), and the phase strip with a new `membersToSegments` adapter that reuses `toWaterfallRows`' own `dominantPhase`/`resolveState` so a card's strip and waterfall never disagree. The body mounts inside the carried `CorpusShell` `<Outlet>`, wrapped in `<PageFrame width="folio">`.

**Tech Stack:** React 18, TypeScript (`exactOptionalPropertyTypes: true`), TanStack Query, react-router-dom, Vitest + RTL, boneyard-js skeletons, Tailwind closed-look design-guard (`npm run lint:design`), d3 (the greenfield waterfall/trace engine).

---

## Standing constraints (do not violate)

- **Provenance (strategy `docs/superpowers/specs/2026-06-03-phase4-cutover-strategy-design.md`):** `src/print/pages/SeriesFolioPage.tsx` + `folioAdapters.ts` import **NEW** (`src/print/**`) + **CARRIED** (`queries.ts`/`api.ts`/`state.ts`/`lib/**`/hooks) **only — never OLD** (`src/components/**`, `src/pages/**`). Enforced by a grep in Task 4. `lib/series/folioFilter.ts` is CARRIED (pure, in `lib/`) — importing it is allowed.
- **Greenfield content, not ported presentation (strategy Rule 3):** assemble from the **mockup** (`docs/redesign-mockups/series-folio.html`) + the `src/print` composites. Do **not** reproduce a legacy affordance just because the legacy folio had it. Scope decisions are fixed in the next section — honor them.
- Commit ONLY specifically-named files. NEVER `git add -A`. The **only** sanctioned `registry.ts`/`*.bones.json` staging is the boneyard recapture in Task 3 (per [[feedback_boneyard_capture_recipe]] — `registry.ts` IS committable for a capture).
- Every commit's exact last line: `Co-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>`
- Work from `packages/HimalayaUI/frontend`. Typecheck `npx tsc --noEmit -p tsconfig.build.json`; design guard `npm run lint:design` (both exit 0). `src/print/pages/**` is NOT design-guard-exempt — placement/layout/named-token classes only.
- Dispatch implementers SEQUENTIALLY (`.git/index.lock`).

---

## Scope decisions (fixed here so they aren't re-litigated)

| Decision | Call | Rationale |
|---|---|---|
| **Pins** (pin/star a series) | ⛔ **OMIT** | No pin affordance anywhere in `series-folio.html`. The `useSeriesPins`/`usePinSeries`/`useUnpinSeries` hooks are unbuilt; the strategy named them as a "folio gap," but Rule 3 binds to the mockup, which has none. Backend routes stay (dormant). |
| **Forks** (fork count / "forked from") | ⛔ **OMIT** | No fork affordance in the mockup. `forked_from_*` fields exist for plumbing but the folio shows none. `useSeriesForksOf` stays unbuilt. |
| **"+N new match" pill** | ⛔ **OMIT** (draft pill only) | The mockup shows `.pill-new` but there is **no backend signal** for "new matches" — the legacy card took it as an unwired prop (always absent). Rendering it would be a control that lies. Only the **draft** pill renders (`content_hash === ""`). If a "new matches" data source is ever added, the pill returns then. |
| **Provenance footer-left** | **Minimal/honest** | The mockup's single-experiment beamtime string ("SSRL · Apr 2026") is **not cheaply derivable** (not in `SeriesSummary`). Render the cross-experiment marker `↔ cross-experiment · q normalized` when `spans_experiments`, else nothing. Richer beamtime provenance is deferred pending a data source. Footer-**right** (edited + author) carries the real metadata. |
| **Figure number** | position-based | `figLabel = draft ? "Recipe" : \`Fig. ${position}\`` where `position` is the 1-based index in the filtered/sorted list (mockup shows "Fig. 1"). Legacy fell back to `series.id`; position reads better and matches the mockup. |
| **Phase-strip derivation** | **fresh adapter, NOT legacy `buildPhaseStrip`** | Legacy `buildPhaseStrip`/`StripSegment` is the OLD single-phase model (`coexistWith: string\|null` from `confirmed_index.phase`); the greenfield `PhaseSegment` is richer (`coexistWith: string[]` + `state`) and `toWaterfallRows` derives phase from `confirmed_phases`. Reusing the legacy helper would make a card's strip disagree with its own waterfall. Write `membersToSegments` reusing `toWaterfallRows`' `dominantPhase`/`resolveState`. |
| **Per-card fetch (N×2 requests)** | accept | Each card does `useSeries(id)` (members) + `useSeriesTraces(id)` (traces) = ~2N requests for N cards — same shape the legacy folio already had (`useSeries` per card), with 4a's batch route collapsing the per-member M. Bounded by card count; TanStack caches/dedupes. Viewport-lazy or list-enrichment is a later optimization (YAGNI). |

---

## Reference: verified APIs (mapped from live source 2026-06-06 — confirm line numbers at execution)

### Greenfield composites (`src/print/components/`, presentational)
```ts
// SeriesCard.tsx (11-36)
interface SeriesCardProps {
  rows: WaterfallRow[]; segments: PhaseSegment[];
  figLabel: string; title: string; sampleCount: number; variable: string;
  provenance: ReactNode; editedLabel: string; author: string;
  notice?: { tone: "new"; count: number } | { tone: "draft" };
  draft?: boolean; onClick?: () => void; className?: string;
}
// Gallery.tsx (8-15): { children?; empty?; className? }  — masonry; page maps & feeds SeriesCard children
// FolioHeader.tsx (4-12): { kicker; title; subtitle?; count; countLabel; className? }
```
### Greenfield engine + types
```ts
// waterfall/waterfallModel.ts
export function toWaterfallRows(members: SeriesMember[], tracesById: Record<number, Trace>): WaterfallRow[]
function dominantPhase(member): string | null   // PRIVATE today → Task 1 exports it
function resolveState(member): AssignmentState  // PRIVATE today → Task 1 exports it
// ui/PhaseStrip.tsx
interface PhaseSegment { phase: string | null; coexistWith?: string[] | null; state?: "form_factor" | "null"; }
```
### CARRIED data + lib
```ts
// queries.ts
useSeriesList()            -> UseQueryResult<SeriesSummary[]>        // 759
useSeries(id?)             -> UseQueryResult<Series>                // 829 (members + samples)
useSeriesTraces(id?)       -> UseQueryResult<Record<number, Trace>> // 843 (Plan 4a)
// lib/series/folioFilter.ts (PURE, CARRIED — import as-is)
filterSort(series: SeriesSummary[], c: FolioControls): SeriesSummary[]
type FolioFilter = "all" | "transition" | "cross";
type FolioSort   = "recent" | "variable" | "size";
interface FolioControls { search: string; filter: FolioFilter; sort: FolioSort; }
// api.ts shapes (verbatim field names)
interface SeriesSummary { id; title; description; content_hash; updated_at; last_event_at; author_username;
  member_count; member_phases: string[]; member_phase_count; has_stale_members; ordering_variable: string|null; spans_experiments: boolean; /* …*/ }
interface Series { id; title; content_hash; ordering_variable; state; members: SeriesMember[]; samples: SeriesSample[]; /*…*/ }
interface SeriesMember { id; exposure_id: number|null; display_order; band_height; y_offset; label_override: string|null; snapshot: MemberSnapshot|null; /*…*/ }
interface MemberSnapshot { effective_peaks; confirmed_index: {…}|null; assignment_state?: AssignmentState; confirmed_phases?: { phase: string; … }[]; analysis_inputs_hash; }
```
### Legacy logic to RE-CREATE (read, do NOT import — it's OLD)
- `formatEdited(iso, now)` — `src/components/SeriesFolioCard.tsx:30-42` (relative-time string; port verbatim into `folioAdapters.ts`).
### Cutover template
- `docs/superpowers/plans/2026-06-04-phase4-samples-cutover.md` — same task shape (adapters → assemble page owning state → repoint+delete OLD → boneyard recapture → gate). `PageFrame width="folio"` already in `PAGE_WIDTHS` (1380).

---

## File map
- **Modify** `src/print/waterfall/waterfallModel.ts` — `export` `dominantPhase` + `resolveState` (no behaviour change).
- **Create** `src/print/pages/folioAdapters.ts` + `test/print-pages/folioAdapters.test.ts` — pure adapters (`membersToSegments`, `toCardChrome`, `formatEdited`).
- **Create** `src/print/pages/SeriesFolioPage.tsx` + `test/print-pages/SeriesFolioPage.test.tsx` — the assembled page + `FolioCard` container.
- **Modify** `src/components/AppRoutes.tsx` — repoint `/series` to the greenfield page.
- **Delete** (grep-guarded) `src/pages/SeriesFolioPage.tsx`, `src/components/SeriesFolioCard.tsx`, `src/components/SeriesMiniWaterfall.tsx`, `src/lib/series/folioFigure.ts` (+ their tests) — all folio-only once `/series` repoints. **Keep** `src/lib/series/folioFilter.ts` (now imported by the greenfield page).
- **Capture/commit** `src/bones/folio.bones.json` (+ the `registry.ts` entry — the sanctioned boneyard edit).

---

## Task 1: Export the phase derivation + the pure folio adapters (TDD)

**Files:** Modify `src/print/waterfall/waterfallModel.ts`; Create `src/print/pages/folioAdapters.ts` + `test/print-pages/folioAdapters.test.ts`.

- [ ] **Step 1: Export the two derivations** — in `waterfallModel.ts`, change `function dominantPhase` → `export function dominantPhase` and `function resolveState` → `export function resolveState` (behaviour unchanged; existing `toWaterfallRows` tests still cover them). Typecheck: `npx tsc --noEmit -p tsconfig.build.json` → 0.

- [ ] **Step 2: Write the failing adapter test** — `test/print-pages/folioAdapters.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import { membersToSegments, toCardChrome, formatEdited } from "../../src/print/pages/folioAdapters";
import type { SeriesMember, SeriesSummary } from "../../src/api";

function member(id: number, snap: SeriesMember["snapshot"]): SeriesMember {
  return { id, series_id: 1, exposure_id: id, display_order: id, band_height: 1, y_offset: 0,
    normalization: "none", color_override: null, label_override: null, q_window_min: null,
    q_window_max: null, peak_display: null, snapshot: snap, is_stale: false, created_by: null, created_at: null };
}
function summary(over: Partial<SeriesSummary>): SeriesSummary {
  return { id: 1, title: "LL37 ratio series", description: null, content_hash: "abc", created_by: null,
    created_at: null, updated_at: "2026-06-04T00:00:00Z", forked_from_id: null, forked_at_hash: null,
    view_grouping_mode: null, view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: "2026-06-04T00:00:00Z", author_username: "alice", member_count: 4,
    member_phases: ["Ia3d", "Im3m"], member_phase_count: 2, has_stale_members: false,
    ordering_variable: "lipid ratio", spans_experiments: false, ...over } as SeriesSummary;
}

describe("membersToSegments", () => {
  it("derives {phase, coexistWith, state} consistent with toWaterfallRows", () => {
    const segs = membersToSegments([
      member(1, { effective_peaks: [], confirmed_index: { phase: "Ia3d" } as any, analysis_inputs_hash: "h",
        confirmed_phases: [{ phase: "Ia3d" } as any, { phase: "Im3m" } as any] }),
      member(2, { effective_peaks: [], confirmed_index: { phase: "Im3m" } as any, analysis_inputs_hash: "h" }),
      member(3, { effective_peaks: [], confirmed_index: null, analysis_inputs_hash: "h", assignment_state: "form_factor" }),
      member(4, null),
    ]);
    expect(segs[0]).toEqual({ phase: "Ia3d", coexistWith: ["Im3m"] });      // coexistence from confirmed_phases tail
    expect(segs[1]).toEqual({ phase: "Im3m" });                            // single phase, no coexist key
    expect(segs[2]).toEqual({ phase: null, state: "form_factor" });        // form-factor → null phase + state
    expect(segs[3]).toEqual({ phase: null });                              // no snapshot → unindexed
  });
});

describe("toCardChrome", () => {
  it("derives saved-card chrome from the summary (position figure, no new-match pill)", () => {
    const c = toCardChrome(summary({}), 2, new Date("2026-06-06T00:00:00Z"));
    expect(c.figLabel).toBe("Fig. 2");
    expect(c.draft).toBe(false);
    expect(c.title).toBe("LL37 ratio series");
    expect(c.sampleCount).toBe(4);
    expect(c.variable).toBe("lipid ratio");
    expect(c.author).toBe("alice");
    expect(c.editedLabel).toBe("2 days ago");
    expect(c.notice).toBeUndefined();   // no "+N new match" — no data source
    expect(c.provenance).toBeNull();    // single-experiment → no marker
  });
  it("marks a draft and a cross-experiment series", () => {
    const c = toCardChrome(summary({ content_hash: "", spans_experiments: true, title: "" }), 1, new Date("2026-06-06T00:00:00Z"));
    expect(c.figLabel).toBe("Recipe");
    expect(c.draft).toBe(true);
    expect(c.notice).toEqual({ tone: "draft" });
    expect(c.title).toBe("Untitled series");
    expect(c.provenance).not.toBeNull();   // cross-experiment marker present
  });
});

describe("formatEdited", () => {
  const now = new Date("2026-06-10T00:00:00Z");
  it("formats relative times", () => {
    expect(formatEdited(null, now)).toBe("recently");
    expect(formatEdited("2026-06-10T00:00:00Z", now)).toBe("just now");
    expect(formatEdited("2026-06-09T00:00:00Z", now)).toBe("yesterday");
    expect(formatEdited("2026-06-06T00:00:00Z", now)).toBe("4 days ago");
    expect(formatEdited("2026-05-27T00:00:00Z", now)).toBe("2 weeks ago");
  });
});
```

- [ ] **Step 3: Run → fail** — `npm test -- print-pages/folioAdapters` → cannot resolve module.

- [ ] **Step 4: Implement** — `src/print/pages/folioAdapters.ts`:

```ts
import type { ReactNode } from "react";
import type { SeriesMember, SeriesSummary } from "../../api";
import type { PhaseSegment } from "../ui/PhaseStrip";
import { dominantPhase, resolveState } from "../waterfall/waterfallModel";

/** Per-member phase-strip segments — SAME phase/state derivation as toWaterfallRows
 *  (imported, not re-derived) so a card's strip and waterfall never disagree.
 *  Coexistence = the confirmed_phases tail beyond the dominant phase. */
export function membersToSegments(members: SeriesMember[]): PhaseSegment[] {
  return members.map((m) => {
    const phase = dominantPhase(m);
    const state = resolveState(m);
    const cp = m.snapshot?.confirmed_phases ?? [];
    const coexistWith = cp.length > 1 ? cp.slice(1).map((p) => p.phase) : null;
    const seg: PhaseSegment = { phase };
    if (coexistWith) seg.coexistWith = coexistWith;             // exactOptionalPropertyTypes: conditional spread
    if (state === "form_factor" || state === "null") seg.state = state;
    return seg;
  });
}

export interface CardChrome {
  figLabel: string; title: string; sampleCount: number; variable: string;
  provenance: ReactNode; editedLabel: string; author: string;
  notice?: { tone: "draft" }; draft: boolean;
}

/** Everything a card shows that is derivable from the LIST summary (no detail fetch).
 *  `position` = 1-based index in the filtered/sorted list (the "Fig. N"). */
export function toCardChrome(s: SeriesSummary, position: number, now: Date): CardChrome {
  const draft = s.content_hash === "";
  const title = s.title.trim() === "" ? "Untitled series" : s.title;
  const chrome: CardChrome = {
    figLabel: draft ? "Recipe" : `Fig. ${position}`,
    title,
    sampleCount: s.member_count,
    variable: s.ordering_variable ?? "",
    // Honest provenance: cross-experiment marker only (single-exp beamtime not derivable). null = omit.
    provenance: s.spans_experiments ? "↔ cross-experiment · q normalized" : null,
    editedLabel: formatEdited(s.updated_at ?? s.last_event_at, now),
    author: s.author_username ?? "",
    draft,
  };
  // No "+N new match" — no backend signal (see plan scope). Only the draft pill.
  if (draft) chrome.notice = { tone: "draft" };
  return chrome;
}

/** Relative-time label (ported verbatim from the legacy SeriesFolioCard). */
export function formatEdited(iso: string | null, now: Date): string {
  if (iso === null) return "recently";
  const then = new Date(iso.replace(" ", "T") + (iso.endsWith("Z") ? "" : "Z"));
  if (Number.isNaN(then.getTime())) return "recently";
  const days = Math.floor((now.getTime() - then.getTime()) / 86_400_000);
  if (days <= 0) return "just now";
  if (days === 1) return "yesterday";
  if (days < 7) return `${days} days ago`;
  if (days < 14) return "1 week ago";
  if (days < 21) return "2 weeks ago";
  return `${Math.floor(days / 7)} weeks ago`;
}
```

> Confirm the legacy `formatEdited` date-parsing (`src/components/SeriesFolioCard.tsx:30-42`) matches — it appends `"Z"` after replacing a space with `T`. The guard above avoids a double-`Z` when the ISO already ends in `Z` (the test feeds `…Z` strings). Mirror legacy behaviour for the space-separated SQLite form.

- [ ] **Step 5: Run → pass** — `npm test -- print-pages/folioAdapters`.
- [ ] **Step 6: Gate** — `npx tsc --noEmit -p tsconfig.build.json`; `npm run lint:design` → 0.
- [ ] **Step 7: Commit**
```bash
git add src/print/waterfall/waterfallModel.ts src/print/pages/folioAdapters.ts test/print-pages/folioAdapters.test.ts
git commit -m "$(printf 'Phase-4 folio: export dominantPhase/resolveState + pure folio adapters\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```

---

## Task 2: Assemble `SeriesFolioPage` + the `FolioCard` container (TDD)

**Files:** Create `src/print/pages/SeriesFolioPage.tsx` + `test/print-pages/SeriesFolioPage.test.tsx`.

**Structure (assemble against `docs/redesign-mockups/series-folio.html`):**
- `<PageFrame width="folio">` wrapping: `FolioHeader` (kicker `"Folio"`, title `"Saved series"`, subtitle = exact `.sub` copy from the mockup, `count={summaries.length}` `countLabel="series in the folio"`) → a controls row (`SearchInput` + three `FilterChip`s `All`/`Has transition`/`Cross-experiment` + a `SegmentedControl` `Recent`/`Variable`/`Largest`) → `Gallery`.
- The page owns `useState` for `search`/`filter`/`sort` (a `FolioControls`), computes `const visible = filterSort(summaries ?? [], controls)`, and maps `visible.map((s, i) => <FolioCard key={s.id} summary={s} position={i + 1} onOpen={() => navigate(\`/series/${s.id}\`)} />)`.
- Empty state: when `visible.length === 0`, pass an `EmptyState` to `Gallery`'s `empty` prop.
- `Skeleton name="folio"` (boneyard, Task 3 captures the fixture) wraps the body while `useSeriesList()` is pending.
- **`FolioCard`** (the per-card data container — NEW, colocated in this file, NOT a print/components composite): calls `useSeries(summary.id)` + `useSeriesTraces(summary.id)`, derives `rows = toWaterfallRows(detail?.members ?? [], traces ?? {})`, `segments = membersToSegments(detail?.members ?? [])`, `chrome = toCardChrome(summary, position, now)`, and renders the presentational `<SeriesCard rows={rows} segments={segments} {...chrome} onClick={onOpen} />`. While detail is pending, `rows`/`segments` are empty (CardFigure renders an empty figure) — a graceful per-card load. `now` is one `new Date()` created at page render and passed down.

- [ ] **Step 1: Write the failing test** — `test/print-pages/SeriesFolioPage.test.tsx`. Mock `../../src/queries` (`useSeriesList`, `useSeries`, `useSeriesTraces`) + `react-router-dom`'s `useNavigate`, mirroring the `LoupePage`/`SamplesPage` page tests. Assert behaviour (testid/role/text, never class strings):
  - renders a `FolioHeader` with the total count and one `SeriesCard` per listed series;
  - typing in the search input filters the visible cards (e.g. 3 series → search narrows to 1);
  - clicking the `Has transition` filter drops single-phase series (`member_phase_count === 1`);
  - clicking a card calls `navigate("/series/:id")`;
  - a series with `content_hash === ""` renders the draft pill ("Draft" / "Recipe");
  - the empty state renders when the filter matches nothing.
  Use real fixtures for `members` so `toWaterfallRows` produces rows (import from `src/print/waterfall/waterfall.fixtures.ts` or build minimal members + a `useSeriesTraces` mock returning a `Record<number, Trace>`).

- [ ] **Step 2: Run → fail** — `npm test -- print-pages/SeriesFolioPage`.
- [ ] **Step 3: Implement** `SeriesFolioPage.tsx` (NEW + CARRIED imports only: `print/components/{FolioHeader,Gallery,SeriesCard}`, `print/ui/{SearchInput,FilterChip,SegmentedControl,EmptyState}`, `print/components/PageFrame`, `print/waterfall/waterfallModel`, `./folioAdapters`, `queries`, `react-router-dom`, `boneyard-js/react`). Placement-only classes.
- [ ] **Step 4: Run → pass** — `npm test -- print-pages/SeriesFolioPage`.
- [ ] **Step 5: Gate** — `npx tsc --noEmit -p tsconfig.build.json`; `npm run lint:design` → 0.
- [ ] **Step 6: Commit** `Phase-4 folio: assemble greenfield SeriesFolioPage + FolioCard container`.

---

## Task 3: Repoint the route, delete the OLD, recapture boneyard

**Files:** Modify `src/components/AppRoutes.tsx`; delete legacy folio files; capture `src/bones/folio.bones.json`.

- [ ] **Step 1: Repoint** — in `AppRoutes.tsx`, change `import { SeriesFolioPage } from "../pages/SeriesFolioPage";` → `import { SeriesFolioPage } from "../print/pages/SeriesFolioPage";` (the `/series` route element is unchanged). `npx tsc --noEmit -p tsconfig.build.json` → 0; app still runs (mixed state).

- [ ] **Step 2: Grep-guard the deletions** —
```bash
grep -rn "SeriesFolioCard\|SeriesMiniWaterfall\|folioFigure\|buildMiniWaterfall\|buildPhaseStrip\|memberPhase" src test e2e
```
Confirm each legacy folio symbol's only remaining importers are the legacy files themselves / their tests (the greenfield page does **not** import any of them — it uses `membersToSegments` + `toWaterfallRows`). `lib/series/folioFilter.ts` must show the **greenfield** page as an importer now — that file is KEPT.

- [ ] **Step 3: Delete the orphans** —
```bash
git rm src/pages/SeriesFolioPage.tsx test/SeriesFolioPage.test.tsx \
       src/components/SeriesFolioCard.tsx test/SeriesFolioCard.test.tsx \
       src/components/SeriesMiniWaterfall.tsx test/SeriesMiniWaterfall.test.tsx \
       src/lib/series/folioFigure.ts test/folioFigure.test.ts
```
(Adjust test paths to the real ones from Step 2's grep. `folioFigure.ts` is deleted **whole** — `buildPhaseStrip`/`memberPhase` were used only by the now-deleted card/mini; the greenfield uses `membersToSegments` instead. If the grep shows any **other** surface still importing `buildPhaseStrip`/`memberPhase`, DEFER `folioFigure.ts` and note it — do not break another page.)

- [ ] **Step 4: Verify whole** — `npx tsc --noEmit -p tsconfig.build.json` (no dangling imports); `npm test -- print-pages folio` (greenfield folio tests pass; the deleted legacy suites are gone).
- [ ] **Step 5: Commit** `Phase-4 folio: route /series -> greenfield + delete legacy folio`.

- [ ] **Step 6: Recapture boneyard** — per [[feedback_boneyard_capture_recipe]]: run the backend-proxied dev server, open `/series` against a cold cache so the boneyard plugin writes `src/bones/folio.bones.json`, verify the skeleton geometry matches the folio masonry, then:
```bash
git add src/bones/folio.bones.json src/bones/registry.ts   # registry.ts IS staged here (sanctioned capture, Rule 6)
git commit -m "$(printf 'Phase-4 folio: capture greenfield folio skeleton\n\nCo-Authored-By: Claude Opus 4.8 (1M context) <noreply@anthropic.com>')"
```
If wiring a `Skeleton name="folio"` requires a `registry.ts` change you did not capture, STOP and surface to the human (do not hand-edit `registry.ts`).

---

## Task 4: Final gate + visual fidelity

- [ ] **Step 1: Provenance grep** — `grep -rn 'from "\.\./\.\./components\|from "\.\./\.\./pages\|/src/components/\|/src/pages/' src/print/pages/SeriesFolioPage.tsx src/print/pages/folioAdapters.ts` → no matches (NEW + CARRIED only).
- [ ] **Step 2: Full local gate** — `npm run lint:design`; `npx tsc --noEmit -p tsconfig.build.json`; `npm test -- print-pages print-components`; `npm run e2e` (mocked — run any series-folio spec; update it if it asserted legacy DOM); `npm run build` → all green.
- [ ] **Step 3: Manual visual fidelity** — `npm run dev`, open `/series`, compare against `docs/redesign-mockups/series-folio.html`: the masonry of distinct-height cards, each card's mini-waterfall (real curves) + kick-row (Fig. N / Recipe + draft pill) + title + meta + phase strip + footer; the header count; search + filter chips + sort segmented control all functional; clicking a card lands on `/series/:id` (the still-legacy builder during migration — expected). Fix placement-only issues in the page, appearance in `ui/` primitives.

---

## Self-review checklist (run before declaring done)
- **Spec coverage:** adapters (T1: `membersToSegments` consistent-with-waterfall, `toCardChrome`, `formatEdited`) · page assembly + per-card container + search/filter/sort wiring (T2) · repoint + grep-guarded OLD deletion + boneyard (T3) · provenance grep + full gate + fidelity (T4). ✅
- **Scope honesty:** pins ⛔, forks ⛔, "+N new match" ⛔ (draft pill only), provenance minimal — all per the Scope-decisions table; the page renders no control without a data source. ✅
- **Single source of truth:** the strip's phase derivation imports `dominantPhase`/`resolveState` from `waterfallModel` — the same functions `toWaterfallRows` uses — so a card's strip and waterfall agree. ✅
- **Provenance:** `SeriesFolioPage`/`folioAdapters` import NEW + CARRIED only (T4.1 grep enforces); `lib/series/folioFilter.ts` is CARRIED (kept); `lib/series/folioFigure.ts` is OLD-only (deleted whole, grep-guarded). ✅
- **Type consistency:** `SeriesCardProps`/`PhaseSegment`/`WaterfallRow`/`SeriesSummary`/`SeriesMember`/`FolioControls` field names match the verified `src/print/**` + `api.ts` + `lib/series/folioFilter.ts` sources; `exactOptionalPropertyTypes` honored via conditional spread of `coexistWith`/`state`/`notice`. ✅
- **No placeholders:** adapters fully coded; the page task gives the structure + load-bearing wiring + a behavioral test (the JSX is assembled against the cited mockup + verified component props, the established page-task pattern). ✅
