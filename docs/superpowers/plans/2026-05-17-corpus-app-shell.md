# Corpus App Shell — Hoisted Routing + Nav-Bridge Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Ship the corpus-scoped HimalayaUI app shell (topbar + stage-tabs + Beamtime chip) and a single hoisted `<Routes>` table so the new shell and the legacy `AppShell` coexist under one router, with a nav-bridge that keeps the URL/`activePage` models consistent and coerces a stale persisted `activePage`.

**Architecture:** One `<Routes>` table (in a new `AppRoutes.tsx`) with two react-router v6 **layout routes** — a pathless `<Route element={<CorpusShell/>}>` for new corpus surfaces and a pathless `<Route element={<AppShell/>}>` for legacy Index/Inspect/Compare. `AppShell` is refactored from a router into a layout component that renders `<Outlet/>`; the legacy URL-sync hooks move into it so they only run on legacy routes (structural isolation — the new shell can never be flagged stale). Stale-`activePage` coercion lands in `state.ts` via a `persist` `merge` callback (no version bump).

**Tech Stack:** React 18, react-router-dom v6.30, Zustand v4.5 (`persist` middleware), TanStack Query v5, Tailwind CSS v4 (`@theme`), Vitest + React Testing Library (JSDOM).

**Spec:** [`docs/superpowers/specs/2026-05-17-corpus-app-shell-design.md`](../specs/2026-05-17-corpus-app-shell-design.md). Issue #155 (I1.1).

**Conventions:**
- All `npm` commands run from `packages/HimalayaUI/frontend/`. All file paths below are relative to the repo root.
- Single-file test run: `node_modules/.bin/vitest run test/<file>` (from the frontend dir).
- Never assert on Tailwind class strings — assert on `data-testid` / `data-*` attributes (see `packages/HimalayaUI/frontend/test/AGENTS.md`).
- Each task ends with one commit. Commit messages reference `#155`.

---

## File structure

| File | Status | Responsibility |
|---|---|---|
| `packages/HimalayaUI/frontend/src/styles.css` | modify | Add "The Print" palette `@theme` color tokens. |
| `packages/HimalayaUI/frontend/src/state.ts` | modify | `coerceActivePage` + `mergePersistedState`; wire `merge` into `persist`. |
| `packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx` | create | New-shell topbar: wordmark, three stage-tabs, Beamtime chip. |
| `packages/HimalayaUI/frontend/src/components/CorpusShell.tsx` | create | New-shell layout route: `CorpusTopbar` + `<Outlet/>`. |
| `packages/HimalayaUI/frontend/src/pages/SamplesPage.tsx` | create | Placeholder body at `/samples`; #160 builds the contact sheet here. |
| `packages/HimalayaUI/frontend/src/components/AppShell.tsx` | modify | Refactor router → layout component; host the relocated URL-sync hooks. |
| `packages/HimalayaUI/frontend/src/components/AppRoutes.tsx` | create | The single hoisted `<Routes>` + shared root effects (theme, shortcuts). |
| `packages/HimalayaUI/frontend/src/App.tsx` | modify | Render `<AppRoutes/>`; drop the relocated URL-sync hook calls. |
| `packages/HimalayaUI/frontend/test/coerceActivePage.test.ts` | create | Unit tests for the coercion + `merge`. |
| `packages/HimalayaUI/frontend/test/CorpusTopbar.test.tsx` | create | Topbar render tests. |
| `packages/HimalayaUI/frontend/test/CorpusShell.test.tsx` | create | Shell + placeholder render tests. |
| `packages/HimalayaUI/frontend/test/AppRoutes.test.tsx` | create | Nav-bridge integration + relocated #77 tests. |
| `packages/HimalayaUI/frontend/test/AppShell.test.tsx` | delete | Cases relocate to `AppRoutes.test.tsx` (AppShell is no longer a router). |
| `packages/HimalayaUI/frontend/test/Compare.routing.test.tsx` | modify | Render `<AppRoutes/>` instead of bare `<AppShell/>`. |
| `packages/HimalayaUI/frontend/test/compareRouteFlatten.test.tsx` | modify | Render `<AppRoutes/>` instead of bare `<AppShell/>`. |

**Not touched:** `parseLocation.ts` (the spec forbids extending its 3-page union), `TabRocker.tsx` (the legacy rocker is untouched until Phase 5), `useStateFromUrl.ts` / `useUrlFromState.ts` (only their call site moves — zero internal change).

---

## Task 1: "The Print" palette tokens

Add the corpus-redesign light-aesthetic color tokens to the Tailwind `@theme` block so the new shell can use `bg-paper`, `text-ink`, `border-hair`, etc. The redesign is light-first; these are single-valued (no dark override). Values are from `docs/redesign-mockups/sample-table.html`.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/styles.css` (the `@theme` block, ~line 23)

- [ ] **Step 1: Add the palette tokens**

In `packages/HimalayaUI/frontend/src/styles.css`, find this exact text inside the `@theme` block:

```css
  --color-peak-manual: oklch(0.78 0.22 340);
  --radius: 6px;
```

Replace it with:

```css
  --color-peak-manual: oklch(0.78 0.22 340);

  /* "The Print" palette — the corpus-redesign light aesthetic (redesign
     master plan §4.1). Consumed by the new CorpusShell topbar. The redesign
     is light-first, so these are single-valued (no theme-light override).
     Values from docs/redesign-mockups/sample-table.html. */
  --color-paper:        oklch(0.978 0.006 85);
  --color-paper-sunk:   oklch(0.951 0.008 84);
  --color-plate:        oklch(0.992 0.004 90);
  --color-ink:          oklch(0.265 0.013 68);
  --color-ink-soft:     oklch(0.467 0.012 70);
  --color-ink-faint:    oklch(0.640 0.010 74);
  --color-hair:         oklch(0.882 0.008 80);
  --color-hair-strong:  oklch(0.806 0.010 78);
  --color-print-accent: oklch(0.555 0.150 38);

  --radius: 6px;
```

- [ ] **Step 2: Verify Tailwind still compiles**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vite build 2>&1 | tail -5
```

Expected: build completes with no CSS error (a `✓ built in …` line). The new tokens are valid CSS custom properties; Tailwind v4 generates `bg-paper` / `text-ink` / `border-hair` utilities from them automatically.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/styles.css
git commit -m "feat(ui): add 'The Print' palette @theme tokens (#155)"
```

---

## Task 2: Stale-`activePage` coercion in `state.ts`

Add a pure `coerceActivePage` and a `mergePersistedState` callback, and wire `merge` into the `persist` config. A stale persisted `activePage` (one naming a retired surface) is coerced to a valid `PageId` at rehydration, so it can never strand the user on an empty `PageBody`. Adding a `merge` callback is **not** a `persist` version bump — `version` stays `3`.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/state.ts`
- Test: `packages/HimalayaUI/frontend/test/coerceActivePage.test.ts`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/frontend/test/coerceActivePage.test.ts`:

```ts
import { describe, it, expect } from "vitest";
import {
  coerceActivePage,
  mergePersistedState,
  useAppState,
} from "../src/state";

describe("coerceActivePage", () => {
  it("passes valid PageIds through unchanged", () => {
    expect(coerceActivePage("index")).toBe("index");
    expect(coerceActivePage("inspect")).toBe("inspect");
    expect(coerceActivePage("compare")).toBe("compare");
  });

  it("coerces an unknown surface name to 'index'", () => {
    // Forward-looking: when #1.7 retires Inspect, a persisted "inspect"
    // would land here. Today any non-PageId string exercises the path.
    expect(coerceActivePage("loupe")).toBe("index");
    expect(coerceActivePage("samples")).toBe("index");
    expect(coerceActivePage("")).toBe("index");
  });

  it("coerces non-string / missing values to 'index'", () => {
    expect(coerceActivePage(undefined)).toBe("index");
    expect(coerceActivePage(null)).toBe("index");
    expect(coerceActivePage(42)).toBe("index");
  });
});

describe("mergePersistedState (persist merge)", () => {
  it("coerces a stale persisted activePage during the merge", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "ghost", theme: "light" },
      current,
    );
    expect(merged.activePage).toBe("index");
  });

  it("keeps a valid persisted activePage and other persisted keys", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(
      { activePage: "compare", theme: "light" },
      current,
    );
    expect(merged.activePage).toBe("compare");
    expect(merged.theme).toBe("light");
  });

  it("preserves store actions through the merge", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState({ activePage: "inspect" }, current);
    expect(typeof merged.setActivePage).toBe("function");
  });

  it("handles undefined persisted state (first run)", () => {
    const current = useAppState.getState();
    const merged = mergePersistedState(undefined, current);
    expect(merged.activePage).toBe("index");
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vitest run test/coerceActivePage.test.ts
```

Expected: FAIL — `coerceActivePage` and `mergePersistedState` are not exported from `../src/state`.

- [ ] **Step 3: Add the coercion helpers to `state.ts`**

In `packages/HimalayaUI/frontend/src/state.ts`, find this exact line:

```ts
export type PageId = "index" | "compare" | "inspect";
```

Immediately after it, add:

```ts

/** The set of `activePage` values that name a live legacy surface. As each
 *  surface is retired (#1.7 Inspect, #3.6 Compare, #4.4 Index), shrink this
 *  set — `coerceActivePage` then redirects a stale persisted value. */
export const VALID_PAGE_IDS: ReadonlySet<PageId> = new Set<PageId>([
  "index",
  "inspect",
  "compare",
]);

/** Coerce an arbitrary persisted value to a valid `PageId`. A value naming a
 *  surface that no longer exists (or a non-string) falls back to "index", so
 *  it can never strand the user on an empty `PageBody` (issue-#77 class). */
export function coerceActivePage(raw: unknown): PageId {
  return typeof raw === "string" && VALID_PAGE_IDS.has(raw as PageId)
    ? (raw as PageId)
    : "index";
}
```

- [ ] **Step 4: Add `mergePersistedState` and wire it into `persist`**

In `packages/HimalayaUI/frontend/src/state.ts`, find the end of the `persist` options object — this exact text:

```ts
      partialize: (s) => ({
        username: s.username,
        firstName: s.firstName,
        lastName: s.lastName,
        activeExperimentId: s.activeExperimentId,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
        activePage: s.activePage,
        tutorialSeen: s.tutorialSeen,
        theme: s.theme,
      }),
    },
  ),
);
```

Replace it with:

```ts
      partialize: (s) => ({
        username: s.username,
        firstName: s.firstName,
        lastName: s.lastName,
        activeExperimentId: s.activeExperimentId,
        activeSampleId: s.activeSampleId,
        activeExposureId: s.activeExposureId,
        activePage: s.activePage,
        tutorialSeen: s.tutorialSeen,
        theme: s.theme,
      }),
      merge: mergePersistedState,
    },
  ),
);

/** persist `merge` — replicates zustand's default shallow merge
 *  ({ ...current, ...persisted }), then coerces a stale persisted
 *  `activePage` so it never enters the store. Adding a `merge` callback is
 *  NOT a `persist` version bump — `version` stays 3. */
export function mergePersistedState(
  persisted: unknown,
  current: AppState,
): AppState {
  const merged = { ...current, ...(persisted as Partial<AppState> | undefined) };
  return { ...merged, activePage: coerceActivePage(merged.activePage) };
}
```

Note: `mergePersistedState` is referenced in the `persist` config before its declaration — that is fine, JavaScript hoists `function` declarations.

- [ ] **Step 5: Run the test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vitest run test/coerceActivePage.test.ts
```

Expected: PASS — all 7 tests green.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/state.ts packages/HimalayaUI/frontend/test/coerceActivePage.test.ts
git commit -m "feat(state): coerce stale persisted activePage via persist merge (#155)"
```

---

## Task 3: `CorpusTopbar` component

The new-shell topbar: corpus wordmark, three stage-tabs, Beamtime facet chip. Per the brainstorming decision, **only the Samples stage-tab is active** (a `<Link to="/samples">`); Index and Series are rendered as disabled `<button>`s until Phases 4 and 3 build those surfaces. The Beamtime chip is a presentational placeholder — `?beamtime=` URL state belongs to #160.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx`
- Test: `packages/HimalayaUI/frontend/test/CorpusTopbar.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/frontend/test/CorpusTopbar.test.tsx`:

```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { CorpusTopbar } from "../src/components/CorpusTopbar";

function renderTopbar() {
  // MemoryRouter — the Samples stage-tab is a react-router <Link>.
  return render(
    <MemoryRouter>
      <CorpusTopbar />
    </MemoryRouter>,
  );
}

describe("CorpusTopbar", () => {
  it("shows the corpus wordmark", () => {
    renderTopbar();
    const wordmark = screen.getByTestId("corpus-wordmark");
    expect(wordmark).toHaveTextContent("Himalaya");
    expect(wordmark).toHaveTextContent("SAXS");
  });

  it("renders three stage-tabs", () => {
    renderTopbar();
    expect(screen.getByTestId("stage-tab-samples")).toBeInTheDocument();
    expect(screen.getByTestId("stage-tab-index")).toBeInTheDocument();
    expect(screen.getByTestId("stage-tab-series")).toBeInTheDocument();
  });

  it("makes Samples the active tab and links it to /samples", () => {
    renderTopbar();
    const samples = screen.getByTestId("stage-tab-samples");
    expect(samples).toHaveAttribute("href", "/samples");
    expect(samples).toHaveAttribute("data-active", "true");
  });

  it("renders Index and Series as inert (disabled) tabs", () => {
    renderTopbar();
    expect(screen.getByTestId("stage-tab-index")).toBeDisabled();
    expect(screen.getByTestId("stage-tab-series")).toBeDisabled();
  });

  it("shows the Beamtime facet chip", () => {
    renderTopbar();
    expect(screen.getByTestId("beamtime-chip")).toHaveTextContent("Beamtime");
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vitest run test/CorpusTopbar.test.tsx
```

Expected: FAIL — cannot resolve `../src/components/CorpusTopbar`.

- [ ] **Step 3: Create the component**

Create `packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx`:

```tsx
import { Link } from "react-router-dom";

interface Stage {
  id: "samples" | "index" | "series";
  label: string;
  /** Absolute path the tab links to. Omitted = inert (disabled) tab. */
  to?: string;
}

// During the Phase-1 migration only the Samples surface exists in the new
// shell. Index and Series are inert tabs until Phases 4 and 3 build those
// surfaces (redesign master plan §2.4).
const STAGES: readonly Stage[] = [
  { id: "samples", label: "Samples", to: "/samples" },
  { id: "index", label: "Index" },
  { id: "series", label: "Series" },
];

/**
 * CorpusTopbar — the topbar for the redesigned corpus-scoped shell: a
 * corpus-level wordmark (experiment is a facet now, not the crumb), the
 * three workflow stage-tabs, and the Beamtime facet chip.
 *
 * The Beamtime chip is a presentational placeholder in #155 — `?beamtime=`
 * URL query state is owned by the `/samples` route (#160).
 */
export function CorpusTopbar(): JSX.Element {
  return (
    <header
      data-testid="corpus-topbar"
      className="h-14 shrink-0 flex items-center gap-4 px-6 border-b border-hair bg-paper"
    >
      <span
        data-testid="corpus-wordmark"
        className="text-sm font-bold uppercase tracking-[0.16em] text-ink"
      >
        Himalaya <span className="font-semibold text-ink-faint">· SAXS</span>
      </span>

      <nav data-testid="stage-tabs" className="flex gap-0.5">
        {STAGES.map((s) => {
          const dot = (
            <span
              aria-hidden="true"
              className={
                "inline-block w-1 h-1 rounded-full mr-1.5 align-middle " +
                (s.to !== undefined ? "bg-print-accent" : "bg-hair-strong")
              }
            />
          );
          return s.to !== undefined ? (
            <Link
              key={s.id}
              to={s.to}
              data-testid={`stage-tab-${s.id}`}
              data-active="true"
              className="px-2.5 py-1.5 rounded text-xs font-semibold uppercase
                         tracking-wide text-ink bg-paper-sunk no-underline"
            >
              {dot}
              {s.label}
            </Link>
          ) : (
            <button
              key={s.id}
              type="button"
              disabled
              data-testid={`stage-tab-${s.id}`}
              className="px-2.5 py-1.5 rounded text-xs font-semibold uppercase
                         tracking-wide text-ink-faint cursor-not-allowed"
            >
              {dot}
              {s.label}
            </button>
          );
        })}
      </nav>

      <button
        type="button"
        data-testid="beamtime-chip"
        title="Filter to a beamtime (coming soon)"
        className="flex items-center gap-1.5 rounded-full border border-hair-strong
                   bg-plate px-2.5 py-1 text-xs font-semibold text-ink"
      >
        <span>Beamtime</span>
        <span aria-hidden="true" className="text-ink-faint">
          ▾
        </span>
      </button>

      <span className="flex-1" />
    </header>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vitest run test/CorpusTopbar.test.tsx
```

Expected: PASS — all 5 tests green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx packages/HimalayaUI/frontend/test/CorpusTopbar.test.tsx
git commit -m "feat(ui): add CorpusTopbar — wordmark, stage-tabs, Beamtime chip (#155)"
```

---

## Task 4: `CorpusShell` layout component + `SamplesPage` placeholder

`CorpusShell` is the layout-route element for new corpus surfaces — it renders the topbar plus an `<Outlet/>` for the matched child route. `SamplesPage` is the placeholder body at `/samples`; #160 (I1.4) builds the real contact sheet into it.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/CorpusShell.tsx`
- Create: `packages/HimalayaUI/frontend/src/pages/SamplesPage.tsx`
- Test: `packages/HimalayaUI/frontend/test/CorpusShell.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `packages/HimalayaUI/frontend/test/CorpusShell.test.tsx`:

```tsx
import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { CorpusShell } from "../src/components/CorpusShell";
import { SamplesPage } from "../src/pages/SamplesPage";

describe("CorpusShell", () => {
  it("renders the topbar and the matched child route in its Outlet", () => {
    render(
      <MemoryRouter initialEntries={["/samples"]}>
        <Routes>
          <Route element={<CorpusShell />}>
            <Route path="/samples" element={<SamplesPage />} />
          </Route>
        </Routes>
      </MemoryRouter>,
    );
    expect(screen.getByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.getByTestId("corpus-topbar")).toBeInTheDocument();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
  });
});

describe("SamplesPage", () => {
  it("renders the placeholder body", () => {
    render(<SamplesPage />);
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vitest run test/CorpusShell.test.tsx
```

Expected: FAIL — cannot resolve `../src/components/CorpusShell` / `../src/pages/SamplesPage`.

- [ ] **Step 3: Create `SamplesPage`**

Create `packages/HimalayaUI/frontend/src/pages/SamplesPage.tsx`:

```tsx
/**
 * SamplesPage — placeholder body for the corpus contact sheet.
 *
 * Issue #155 (I1.1) only establishes the route slot so the corpus shell has
 * a body to mount. Issue #160 (I1.4) builds the real contact-sheet table
 * into this file.
 */
export function SamplesPage(): JSX.Element {
  return (
    <div data-testid="samples-page" className="p-10 text-sm text-ink-faint">
      Contact sheet — coming soon.
    </div>
  );
}
```

- [ ] **Step 4: Create `CorpusShell`**

Create `packages/HimalayaUI/frontend/src/components/CorpusShell.tsx`:

```tsx
import { Outlet } from "react-router-dom";
import { CorpusTopbar } from "./CorpusTopbar";

/**
 * CorpusShell — the layout-route element for the redesigned corpus-scoped
 * surfaces. Renders the topbar and an <Outlet/> for the matched child route.
 *
 * Registered as a pathless layout route in `AppRoutes`; later issues add
 * child routes under it (#161 the loupe, #179 the focus workspace).
 */
export function CorpusShell(): JSX.Element {
  return (
    <div
      data-testid="corpus-shell"
      className="h-full w-full flex flex-col min-h-0 bg-paper text-ink"
    >
      <CorpusTopbar />
      <main className="flex-1 min-h-0 overflow-auto">
        <Outlet />
      </main>
    </div>
  );
}
```

- [ ] **Step 5: Run the test to verify it passes**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vitest run test/CorpusShell.test.tsx
```

Expected: PASS — both tests green.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/CorpusShell.tsx packages/HimalayaUI/frontend/src/pages/SamplesPage.tsx packages/HimalayaUI/frontend/test/CorpusShell.test.tsx
git commit -m "feat(ui): add CorpusShell layout + SamplesPage placeholder (#155)"
```

---

## Task 5: Hoist the router — `AppRoutes`, refactor `AppShell`, rewire `App`

This is the atomic structural change: extract one hoisted `<Routes>` into a new `AppRoutes.tsx`, refactor `AppShell` from a router into a layout component, and point `App.tsx` at `AppRoutes`. The legacy URL-sync hooks (`useStateFromUrl`/`useUrlFromState`) move from `App.tsx` into `AppShell` so they only mount under legacy routes — the new shell can never trigger the legacy stale-path logic. Three test files that render bare `<AppShell/>` as a router are updated in the same task because the change is atomic.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/AppRoutes.tsx`
- Modify: `packages/HimalayaUI/frontend/src/components/AppShell.tsx`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx`
- Create: `packages/HimalayaUI/frontend/test/AppRoutes.test.tsx`
- Delete: `packages/HimalayaUI/frontend/test/AppShell.test.tsx`
- Modify: `packages/HimalayaUI/frontend/test/Compare.routing.test.tsx`
- Modify: `packages/HimalayaUI/frontend/test/compareRouteFlatten.test.tsx`

- [ ] **Step 1: Write the new nav-bridge integration test**

Create `packages/HimalayaUI/frontend/test/AppRoutes.test.tsx`:

```tsx
/**
 * AppRoutes — the single hoisted route table. Tests the nav-bridge: new
 * routes mount the corpus shell, legacy routes mount AppShell, and the two
 * nav models do not desync. Includes the relocated #77 compare-sync tests
 * (formerly AppShell.test.tsx — AppShell is no longer a router).
 */
import { describe, it, expect, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { useAppState } from "../src/state";
import { AppRoutes } from "../src/components/AppRoutes";

function makeQc() {
  return new QueryClient({
    defaultOptions: {
      queries: { retry: false, gcTime: Infinity, staleTime: 0 },
      mutations: { retry: false },
    },
  });
}

function renderRoutes(initialPath: string, initialIndex?: number) {
  const qc = makeQc();
  const entries = initialIndex !== undefined
    ? { initialEntries: [initialPath, "/"], initialIndex }
    : { initialEntries: [initialPath] };
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter {...entries}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("AppRoutes — nav-bridge shell selection", () => {
  beforeEach(() => {
    useAppState.setState({ activePage: "index", activeExperimentId: undefined });
  });

  it("mounts the corpus shell (not AppShell) at /samples", async () => {
    renderRoutes("/samples");
    expect(await screen.findByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("app-shell")).toBeNull();
  });

  it("mounts AppShell (not the corpus shell) at /index", async () => {
    renderRoutes("/index");
    expect(await screen.findByTestId("app-shell")).toBeInTheDocument();
    expect(screen.queryByTestId("corpus-shell")).toBeNull();
  });

  it("does not flag /samples as a stale path", async () => {
    // The legacy URL-sync hooks live in AppShell, which is not mounted on a
    // corpus route — so /samples cannot be parsed as a stale path.
    renderRoutes("/samples");
    await screen.findByTestId("corpus-shell");
    expect(useAppState.getState().staleUrlContext).toBeNull();
  });
});

describe("AppRoutes — Zustand → URL compare-sync (#77)", () => {
  beforeEach(() => {
    useAppState.setState({ activePage: "index", activeExperimentId: undefined });
  });

  it("activePage='compare' + URL '/' navigates to /compare/all", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: undefined });
    renderRoutes("/");
    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });

  it("activePage='compare' + URL '/' navigates to /experiments/:eid/compare when an experiment is set", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: 7 });
    renderRoutes("/");
    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });

  it("activePage='index' + URL '/' does NOT navigate to a compare route", async () => {
    useAppState.setState({ activePage: "index", activeExperimentId: undefined });
    renderRoutes("/");
    await new Promise((resolve) => setTimeout(resolve, 20));
    expect(screen.queryByTestId("compare-page")).toBeNull();
  });

  it("back-nav from /experiments/:eid/compare/:id to '/' bounces to a compare URL (intentional)", async () => {
    useAppState.setState({ activePage: "compare", activeExperimentId: 7 });
    renderRoutes("/experiments/7/compare/123", 1);
    await waitFor(() => {
      expect(screen.getByTestId("compare-page")).toBeInTheDocument();
    });
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vitest run test/AppRoutes.test.tsx
```

Expected: FAIL — cannot resolve `../src/components/AppRoutes`.

- [ ] **Step 3: Refactor `AppShell` into a layout component**

Replace the **entire contents** of `packages/HimalayaUI/frontend/src/components/AppShell.tsx` with:

```tsx
import { useEffect } from "react";
import { Outlet, useLocation, useNavigate } from "react-router-dom";
import { useAppState } from "../state";
import { AppHeader } from "./AppHeader";
import { TabRocker } from "./TabRocker";
import { NavModal } from "./NavModal";
import { useStateFromUrl } from "../hooks/useStateFromUrl";
import { useUrlFromState } from "../hooks/useUrlFromState";

/**
 * AppShell — the legacy layout: the grain background, the app header, the
 * page-nav rocker, and an <Outlet/> for the matched legacy route (Index /
 * Inspect / Compare). Registered as a pathless layout route in `AppRoutes`.
 *
 * Phase 5 (#182) deletes this; until then it is the migration rollback path.
 *
 * Routing model: AppShell is the element of the *legacy* layout route, so it
 * mounts only when a legacy URL matches. The URL↔Zustand sync hooks live
 * here (not at the app root) precisely so they never run on a corpus route —
 * the new shell owns its own URL and the legacy stale-path logic cannot
 * strand the user there.
 */
export function AppShell(): JSX.Element {
  // URL ↔ Zustand sync, relocated from App.tsx. Order matters:
  // useStateFromUrl populates Zustand before useUrlFromState reflects it
  // back — the `resolving`-flag handshake between the two depends on
  // useStateFromUrl's effect running first.
  useStateFromUrl();
  useUrlFromState();

  const experimentId = useAppState((s) => s.activeExperimentId);
  const setActivePage = useAppState((s) => s.setActivePage);
  const activePage = useAppState((s) => s.activePage);
  const location = useLocation();
  const navigate = useNavigate();

  // Sync URL → Zustand activePage. When the URL is a compare path, mark the
  // page tab as "compare". On other paths we leave activePage alone.
  const onComparePath =
    location.pathname.startsWith("/compare") ||
    /^\/experiments\/\d+\/compare(\/|$)/.test(location.pathname);
  useEffect(() => {
    if (onComparePath) setActivePage("compare");
  }, [onComparePath, setActivePage]);

  // Symmetric: when activePage is "compare" but the URL isn't on a compare
  // path, navigate so the URL-routed Compare page mounts. Without this, a
  // reload at "/" with persisted activePage='compare' renders the rocker but
  // no page body (issue #77).
  useEffect(() => {
    if (activePage !== "compare") return;
    if (onComparePath) return;
    const url =
      experimentId !== undefined
        ? `/experiments/${experimentId}/compare`
        : "/compare/all";
    navigate(url, { replace: true });
  }, [activePage, onComparePath, experimentId, navigate]);

  return (
    <div
      data-testid="app-shell"
      className="h-full w-full max-w-[1600px] mx-auto flex flex-col min-h-0 relative"
      // --chrome-h: AppHeader (h-11 = 44px) + TabRocker row (~40px).
      style={{ "--chrome-h": "84px" } as React.CSSProperties}
    >
      <AppHeader />
      <div className="shrink-0 flex justify-center pt-1 pb-2">
        <TabRocker
          experimentId={experimentId}
          onNavigateAway={(target) => {
            // Leaving Compare → return to "/" so the legacy body renders the
            // chosen page.
            if (onComparePath && target !== "compare") navigate("/");
          }}
        />
      </div>
      <Outlet />
      <NavModal />
    </div>
  );
}
```

This drops from `AppShell`: the internal `<Routes>`, the `theme` effect, `useGlobalShortcuts` + `useSamples` (all move to `AppRoutes`), and `PageBody` / `EditToBareRedirect` (move to `AppRoutes`).

- [ ] **Step 4: Create `AppRoutes`**

Create `packages/HimalayaUI/frontend/src/components/AppRoutes.tsx`:

```tsx
import { useEffect } from "react";
import { Routes, Route, Navigate, useLocation } from "react-router-dom";
import { useAppState } from "../state";
import { useSamples } from "../queries";
import { useGlobalShortcuts } from "../hooks/useGlobalShortcuts";
import { CorpusShell } from "./CorpusShell";
import { SamplesPage } from "../pages/SamplesPage";
import { AppShell } from "./AppShell";
import { IndexPage } from "../pages/IndexPage";
import { InspectPage } from "../pages/InspectPage";
import { Compare } from "../pages/Compare";
import { ResolvingFallback } from "./ResolvingFallback";
import { StaleUrlPage } from "./StaleUrlPage";

/**
 * PageBody — the legacy Index/Inspect body switcher, driven by Zustand
 * `activePage`. Compare URLs are matched by their explicit <Route> entries
 * below, so activePage === "compare" never reaches here.
 */
function PageBody(): JSX.Element {
  const activePage = useAppState((s) => s.activePage);
  const resolving = useAppState((s) => s.resolving);
  const staleUrlContext = useAppState((s) => s.staleUrlContext);

  if (resolving) return <ResolvingFallback />;
  if (staleUrlContext !== null) {
    return <StaleUrlPage staleUrlContext={staleUrlContext} />;
  }
  if (activePage === "index") return <IndexPage />;
  if (activePage === "inspect") return <InspectPage />;
  return <></>;
}

/**
 * EditToBareRedirect — the `/edit` URL segment is gone (Compare UX Phase B).
 * Old `/edit` deep-links resolve by redirecting to the bare path. Reads
 * `useLocation()` (the router store), not `window.location` — under
 * MemoryRouter the latter stays at "/".
 */
function EditToBareRedirect(): JSX.Element {
  const loc = useLocation();
  const pathname = loc.pathname.replace(/\/edit\/?$/, "");
  return (
    <Navigate
      to={{ pathname, search: loc.search, hash: loc.hash }}
      replace
    />
  );
}

/**
 * AppRoutes — the single hoisted top-level <Routes> table, plus the shared
 * root effects (theme, global shortcuts) that sit above both shell bodies.
 *
 * Two pathless layout routes: <CorpusShell> for new corpus surfaces and
 * <AppShell> for the legacy Index/Inspect/Compare surfaces. Later redesign
 * issues register their route slot under the corpus layout route (#161 the
 * loupe, #179 the focus workspace).
 */
export function AppRoutes(): JSX.Element {
  const theme = useAppState((s) => s.theme);
  const experimentId = useAppState((s) => s.activeExperimentId);

  // Apply theme to <html> so CSS can key off `html.theme-light`. Hoisted
  // here (above both shell bodies) so the theme works under either shell.
  useEffect(() => {
    document.documentElement.className =
      theme === "light" ? "theme-light" : "";
    return () => {
      document.documentElement.className = "";
    };
  }, [theme]);

  // Global keyboard shortcuts — hoisted above both shell bodies so they work
  // under either. `useSamples(experimentId ?? 0)` matches the prior AppShell
  // call site; the `,`/`.` sample-step shortcut needs the active
  // experiment's samples.
  const samplesQ = useSamples(experimentId ?? 0);
  useGlobalShortcuts(experimentId === undefined ? undefined : samplesQ.data);

  return (
    <Routes>
      <Route element={<CorpusShell />}>
        <Route path="/samples" element={<SamplesPage />} />
      </Route>
      <Route element={<AppShell />}>
        <Route path="/" element={<PageBody />} />
        <Route path="/index" element={<PageBody />} />
        <Route path="/index/:experiment" element={<PageBody />} />
        <Route path="/index/:experiment/:sample" element={<PageBody />} />
        <Route path="/inspect" element={<PageBody />} />
        <Route path="/inspect/:experiment" element={<PageBody />} />
        <Route path="/inspect/:experiment/:sample" element={<PageBody />} />
        <Route path="/experiments/:eid/compare" element={<Compare />} />
        <Route path="/experiments/:eid/compare/new" element={<Compare />} />
        <Route path="/experiments/:eid/compare/:id" element={<Compare />} />
        <Route
          path="/experiments/:eid/compare/:id/edit"
          element={<EditToBareRedirect />}
        />
        <Route path="/compare/all" element={<Compare />} />
        <Route path="/compare/all/new" element={<Compare />} />
        <Route path="/compare/all/:id" element={<Compare />} />
        <Route path="/compare/all/:id/edit" element={<EditToBareRedirect />} />
        <Route path="*" element={<PageBody />} />
      </Route>
    </Routes>
  );
}
```

- [ ] **Step 5: Rewire `App.tsx`**

In `packages/HimalayaUI/frontend/src/App.tsx`:

(a) Replace the import line:

```tsx
import { AppShell } from "./components/AppShell";
```

with:

```tsx
import { AppRoutes } from "./components/AppRoutes";
```

(b) Delete these two import lines entirely:

```tsx
import { useStateFromUrl } from "./hooks/useStateFromUrl";
import { useUrlFromState } from "./hooks/useUrlFromState";
```

(c) Delete this block from the body of `App()` (the comment plus the two hook calls):

```tsx
  // Permalink URL ↔ Zustand sync. Order matters on cold mount —
  // useStateFromUrl populates Zustand from the address bar first; then
  // useUrlFromState's equality guard makes ordering irrelevant after the
  // first render.
  useStateFromUrl();    // URL → Zustand
  useUrlFromState();    // Zustand → URL

```

(d) In the returned JSX, replace `<AppShell />` with `<AppRoutes />`. The final return becomes:

```tsx
  return (
    <>
      <AppRoutes />
      <OnboardingFlow />
      <ConflictModal />
      <ToastContainer />
      <InfrastructureBanner />
    </>
  );
```

- [ ] **Step 6: Delete the obsolete `AppShell.test.tsx`**

`AppShell.test.tsx` renders bare `<AppShell/>` expecting it to be a router; its four cases were relocated into `AppRoutes.test.tsx` in Step 1.

```bash
git rm packages/HimalayaUI/frontend/test/AppShell.test.tsx
```

- [ ] **Step 7: Update `Compare.routing.test.tsx` to render `<AppRoutes/>`**

In `packages/HimalayaUI/frontend/test/Compare.routing.test.tsx`:

(a) Replace the import line:

```tsx
import { AppShell } from "../src/components/AppShell";
```

with:

```tsx
import { AppRoutes } from "../src/components/AppRoutes";
```

(b) In the `renderShell` helper, replace `<AppShell />` with `<AppRoutes />`:

```tsx
function renderShell(initialPath: string) {
  const qc = makeQc();
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[initialPath]}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}
```

- [ ] **Step 8: Update `compareRouteFlatten.test.tsx` to render `<AppRoutes/>`**

In `packages/HimalayaUI/frontend/test/compareRouteFlatten.test.tsx`:

(a) Replace the import line:

```tsx
import { AppShell } from "../src/components/AppShell";
```

with:

```tsx
import { AppRoutes } from "../src/components/AppRoutes";
```

(b) Replace the entire `renderShell` helper. The old helper wrapped `<AppShell/>` in its own `<Routes><Route path="*">` — `AppRoutes` now owns the `<Routes>`, so `LocationProbe` becomes a plain sibling under `MemoryRouter`:

```tsx
function renderShell(initialPath: string) {
  const qc = makeQc();
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[initialPath]}>
        <LocationProbe />
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}
```

The unused `Routes` and `Route` imports from `react-router-dom` can stay or be removed; if `tsc` (Step 11) flags them as unused, remove them from the import on line 10, leaving `import { MemoryRouter, useLocation } from "react-router-dom";`.

- [ ] **Step 9: Run the new and updated routing tests**

Run (from `packages/HimalayaUI/frontend/`):

```bash
node_modules/.bin/vitest run test/AppRoutes.test.tsx test/Compare.routing.test.tsx test/compareRouteFlatten.test.tsx
```

Expected: PASS — all three files green. The nav-bridge tests confirm `/samples` mounts `CorpusShell` and `/index` mounts `AppShell`; the relocated #77 tests confirm the compare-sync still works through the hoisted router.

- [ ] **Step 10: Run the full Vitest suite**

The router hoist touches the app root; `smoke.test.tsx` and `sampleSwitchKeypress.test.tsx` render `<App/>` and must stay green (legacy routing behavior is unchanged).

Run (from `packages/HimalayaUI/frontend/`):

```bash
npm test > /tmp/vitest-t5.out 2>&1; grep -E "Test Files|Tests |FAIL" /tmp/vitest-t5.out
```

Expected: no `FAIL` lines; the `Test Files` / `Tests` summary lines show all passing. If `smoke.test.tsx` or `sampleSwitchKeypress.test.tsx` fails, the refactor changed legacy behavior — stop and reconcile before committing (the legacy routes and the relocated hooks must behave exactly as before).

- [ ] **Step 11: Run the type-check + build**

Run (from `packages/HimalayaUI/frontend/`):

```bash
npm run build 2>&1 | tail -8
```

Expected: `tsc --noEmit` reports no errors and `vite build` completes (`✓ built in …`). If `tsc` flags unused `Routes`/`Route` imports in `compareRouteFlatten.test.tsx`, apply the fix noted in Step 8(b).

- [ ] **Step 12: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/AppRoutes.tsx \
        packages/HimalayaUI/frontend/src/components/AppShell.tsx \
        packages/HimalayaUI/frontend/src/App.tsx \
        packages/HimalayaUI/frontend/test/AppRoutes.test.tsx \
        packages/HimalayaUI/frontend/test/Compare.routing.test.tsx \
        packages/HimalayaUI/frontend/test/compareRouteFlatten.test.tsx
git commit -m "feat(ui): hoist a single Routes table over corpus + legacy shells (#155)"
```

---

## Verification (whole feature)

After Task 5, confirm against the issue #155 acceptance criteria:

- [ ] A single top-level `<Routes>` (`AppRoutes.tsx`) renders `AppShell` for legacy paths and `CorpusShell` for `/samples` — proven by `AppRoutes.test.tsx`.
- [ ] The new shell shows the wordmark, three stage-tabs, and the Beamtime chip — proven by `CorpusTopbar.test.tsx`.
- [ ] `useGlobalShortcuts` and the theme effect run in `AppRoutes`, above both shell bodies.
- [ ] URL↔`activePage` stays consistent; the relocated #77 tests are green.
- [ ] A stale persisted `activePage` is coerced (`coerceActivePage` / `mergePersistedState`) — proven by `coerceActivePage.test.ts`.
- [ ] `persist` `version` is still `3` — no bump (verify in `state.ts`).
- [ ] `npm test` is fully green and `npm run build` passes.
