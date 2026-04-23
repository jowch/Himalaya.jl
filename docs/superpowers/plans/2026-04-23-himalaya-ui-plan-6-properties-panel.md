# HimalayaUI Plan 6 — Tabbed Properties Panel Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Replace the single-purpose exposure list in the bottom-center pane with a tabbed Properties panel containing Exposures, Peaks, Tags, and Notes.

**Architecture:** One new `PropertiesPanel` component owns tab state and renders a header of tab buttons plus the active tab's content. Four small tab components — `ExposuresTab` (wraps the existing `ExposureList`), `PeaksTab`, `TagsTab`, `NotesTab` — each use TanStack Query for server state and Zustand selectors for context. Sample-level mutations (`useUpdateSample`, `useAddSampleTag`, `useRemoveSampleTag`) invalidate the samples query so the navbar breadcrumb reacts to renames.

**Tech Stack:** React 18, TypeScript strict (`exactOptionalPropertyTypes: true`), Vite, Zustand, TanStack Query v5, Tailwind v4, Vitest + RTL, Playwright. No backend changes — all four REST endpoints already exist.

**Assumptions about the starting state** (verified at plan-write time):
- `main` is at `160e5c9` — Plan 5 merged. Backend: 198 tests pass. Frontend: 67 unit tests + 6 E2E pass.
- Backend endpoints (all existing): `GET /api/experiments/:id/samples` → `[{id, label, name, notes, tags}]`; `PATCH /api/samples/:id` with body `{name?, notes?}`; `POST /api/samples/:id/tags` with body `{key, value}`; `DELETE /api/samples/:id/tags/:tag_id`. All mutating endpoints read `X-Username`.
- `Sample` type in `api.ts` already includes `notes: string | null` and `tags: SampleTag[]`.
- `api.ts` already exports `updateSample`, `addSampleTag`, `removeSampleTag` fetchers. No sample-level query hooks yet (samples-for-experiment is a single list hook; there's no per-sample refetch hook).
- `queries.ts` has `queryKeys.samples(experimentId)`. Sample mutations must invalidate this key.
- `state.ts` has `activeSampleId` and `activeExposureId` in Zustand.
- `components/ExposureList.tsx` already renders an exposures list. Plan 6 wraps it in an `ExposuresTab` thin shim — no logic change.

**Out of scope for Plan 6 (deferred):**
- The Phase panel's **Recent** section (audit trail per exposure). Spec calls it "not interactive in v1". Separate plan.
- Editing exposure tags (only sample tags here — exposure tags get their own UI later).
- Deriving averaged/background-subtracted exposures (Plan 4's backend supports file-kind only; derived construction is a future feature).

**Commands run from:** `packages/HimalayaUI/frontend/`.

---

## File Map

- Modify: `packages/HimalayaUI/frontend/src/queries.ts` (add `useUpdateSample`, `useAddSampleTag`, `useRemoveSampleTag`)
- Create: `packages/HimalayaUI/frontend/src/components/PropertiesPanel.tsx`
- Create: `packages/HimalayaUI/frontend/src/components/ExposuresTab.tsx` (thin wrapper around existing `ExposureList`)
- Create: `packages/HimalayaUI/frontend/src/components/PeaksTab.tsx`
- Create: `packages/HimalayaUI/frontend/src/components/TagsTab.tsx`
- Create: `packages/HimalayaUI/frontend/src/components/NotesTab.tsx`
- Modify: `packages/HimalayaUI/frontend/src/App.tsx` (swap `centerBottom={<ExposureList />}` → `<PropertiesPanel />`)
- Create tests: `test/{queries-samples,PropertiesPanel,PeaksTab,TagsTab,NotesTab}.test.tsx`
- Modify: `packages/HimalayaUI/frontend/e2e/smoke.spec.ts` (tab switching + add tag)

No backend changes. No CLAUDE.md changes.

---

## Task 1: Sample mutation hooks

**Motivation:** Enable the Tags and Notes tabs to trigger PATCH/POST/DELETE against the samples endpoints with proper cache invalidation. Mutations invalidate `queryKeys.samples(EXPERIMENT_ID)` — samples-for-experiment is the single query that holds all sample data including tags and notes.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/queries.ts`
- Create: `packages/HimalayaUI/frontend/test/queries-samples.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/queries-samples.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { renderHook, act } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import {
  useUpdateSample, useAddSampleTag, useRemoveSampleTag, queryKeys,
} from "../src/queries";

const EXPERIMENT_ID = 1;

function withClient() {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return { client, wrapper };
}

function mockOnce(status: number, body: unknown): void {
  vi.spyOn(global, "fetch").mockResolvedValueOnce(
    new Response(status === 204 ? null : JSON.stringify(body), {
      status, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("queries — sample mutations", () => {
  beforeEach(() => { vi.restoreAllMocks(); });

  it("useUpdateSample invalidates samples on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(200, { id: 10, experiment_id: 1, label: "A1", name: "s1",
                    notes: "updated", tags: [] });
    const { result } = renderHook(
      () => useUpdateSample(EXPERIMENT_ID, 10), { wrapper },
    );
    await act(async () => { await result.current.mutateAsync({ notes: "updated" }); });
    expect(invalidate).toHaveBeenCalledWith({
      queryKey: queryKeys.samples(EXPERIMENT_ID),
    });
  });

  it("useAddSampleTag invalidates samples on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(201, { id: 5, sample_id: 10, key: "lipid", value: "DOPC",
                    source: "manual" });
    const { result } = renderHook(
      () => useAddSampleTag(EXPERIMENT_ID, 10), { wrapper },
    );
    await act(async () => {
      await result.current.mutateAsync({ key: "lipid", value: "DOPC" });
    });
    expect(invalidate).toHaveBeenCalledWith({
      queryKey: queryKeys.samples(EXPERIMENT_ID),
    });
  });

  it("useRemoveSampleTag invalidates samples on success", async () => {
    const { client, wrapper } = withClient();
    const invalidate = vi.spyOn(client, "invalidateQueries");
    mockOnce(204, null);
    const { result } = renderHook(
      () => useRemoveSampleTag(EXPERIMENT_ID, 10), { wrapper },
    );
    await act(async () => { await result.current.mutateAsync(5); });
    expect(invalidate).toHaveBeenCalledWith({
      queryKey: queryKeys.samples(EXPERIMENT_ID),
    });
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- queries-samples`
Expected: 3 tests fail (hooks not defined).

- [ ] **Step 3: Extend `queries.ts`**

Append to `packages/HimalayaUI/frontend/src/queries.ts` (below the existing mutation hooks). Reuse the existing `authOpts(username)` helper.

```ts
export function useUpdateSample(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (patch: { name?: string; notes?: string }) =>
      api.updateSample(sampleId, patch, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) }),
  });
}

export function useAddSampleTag(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: ({ key, value }: { key: string; value: string }) =>
      api.addSampleTag(sampleId, key, value, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) }),
  });
}

export function useRemoveSampleTag(experimentId: number, sampleId: number) {
  const qc = useQueryClient();
  const username = useAppState((s) => s.username);
  return useMutation({
    mutationFn: (tagId: number) =>
      api.removeSampleTag(sampleId, tagId, authOpts(username)),
    onSuccess: () =>
      qc.invalidateQueries({ queryKey: queryKeys.samples(experimentId) }),
  });
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test && npm run build`
Expected: all green.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/queries.ts \
        packages/HimalayaUI/frontend/test/queries-samples.test.tsx
git commit -m "feat(ui): sample mutation hooks (update, add/remove tag)"
```

---

## Task 2: PropertiesPanel shell + tab state

**Motivation:** Stand up the skeleton — tab header, active-tab state, placeholder body per tab. Individual tabs get filled in by later tasks but this lets the app compile with the new component in place.

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/PropertiesPanel.tsx`
- Create: `packages/HimalayaUI/frontend/test/PropertiesPanel.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/PropertiesPanel.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PropertiesPanel } from "../src/components/PropertiesPanel";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({
    activeSampleId: undefined,
    activeExposureId: undefined,
  });
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response("[]", { status: 200, headers: { "Content-Type": "application/json" } }),
  );
});

describe("<PropertiesPanel>", () => {
  it("renders four tab buttons", () => {
    renderWithProviders(<PropertiesPanel />);
    expect(screen.getByRole("tab", { name: /exposures/i })).toBeInTheDocument();
    expect(screen.getByRole("tab", { name: /peaks/i })).toBeInTheDocument();
    expect(screen.getByRole("tab", { name: /tags/i })).toBeInTheDocument();
    expect(screen.getByRole("tab", { name: /notes/i })).toBeInTheDocument();
  });

  it("starts on Exposures tab", () => {
    renderWithProviders(<PropertiesPanel />);
    expect(screen.getByRole("tab", { name: /exposures/i }))
      .toHaveAttribute("aria-selected", "true");
  });

  it("switching tabs updates aria-selected", () => {
    renderWithProviders(<PropertiesPanel />);
    fireEvent.click(screen.getByRole("tab", { name: /peaks/i }));
    expect(screen.getByRole("tab", { name: /peaks/i }))
      .toHaveAttribute("aria-selected", "true");
    expect(screen.getByRole("tab", { name: /exposures/i }))
      .toHaveAttribute("aria-selected", "false");
  });

  it("tab panel has matching aria-labelledby for the active tab", () => {
    renderWithProviders(<PropertiesPanel />);
    fireEvent.click(screen.getByRole("tab", { name: /tags/i }));
    const panel = screen.getByRole("tabpanel");
    const activeTab = screen.getByRole("tab", { name: /tags/i });
    expect(panel.getAttribute("aria-labelledby")).toBe(activeTab.id);
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- PropertiesPanel`
Expected: 4 tests fail (component not defined).

- [ ] **Step 3: Implement `PropertiesPanel.tsx`**

Create `packages/HimalayaUI/frontend/src/components/PropertiesPanel.tsx`:

```tsx
import { useState } from "react";
import { ExposuresTab } from "./ExposuresTab";
import { PeaksTab } from "./PeaksTab";
import { TagsTab } from "./TagsTab";
import { NotesTab } from "./NotesTab";

type TabKey = "exposures" | "peaks" | "tags" | "notes";

interface TabDef {
  key: TabKey;
  label: string;
  Component: () => JSX.Element;
}

const TABS: readonly TabDef[] = [
  { key: "exposures", label: "Exposures", Component: ExposuresTab },
  { key: "peaks",     label: "Peaks",     Component: PeaksTab     },
  { key: "tags",      label: "Tags",      Component: TagsTab      },
  { key: "notes",     label: "Notes",     Component: NotesTab     },
] as const;

export function PropertiesPanel(): JSX.Element {
  const [active, setActive] = useState<TabKey>("exposures");
  const ActiveBody = TABS.find((t) => t.key === active)!.Component;

  return (
    <div className="flex flex-col h-full">
      <div role="tablist" className="flex gap-1 border-b border-border pb-1 mb-2">
        {TABS.map((t) => {
          const selected = t.key === active;
          return (
            <button
              key={t.key}
              id={`tab-${t.key}`}
              role="tab"
              aria-selected={selected}
              aria-controls={`tabpanel-${t.key}`}
              className={
                "px-3 py-1 rounded-md text-[13px] focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent " +
                (selected
                  ? "bg-bg-elevated text-fg border-b-2 border-accent"
                  : "text-fg-muted hover:text-fg hover:bg-bg-hover")
              }
              onClick={() => setActive(t.key)}
            >
              {t.label}
            </button>
          );
        })}
      </div>
      <div
        role="tabpanel"
        id={`tabpanel-${active}`}
        aria-labelledby={`tab-${active}`}
        className="flex-1 overflow-auto"
      >
        <ActiveBody />
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Create stub tab files**

To get the tests passing before implementing the bodies, create four stubs. Each is one line.

Create `packages/HimalayaUI/frontend/src/components/ExposuresTab.tsx`:

```tsx
export function ExposuresTab(): JSX.Element { return <p className="text-fg-muted italic">Exposures.</p>; }
```

Create `packages/HimalayaUI/frontend/src/components/PeaksTab.tsx`:

```tsx
export function PeaksTab(): JSX.Element { return <p className="text-fg-muted italic">Peaks.</p>; }
```

Create `packages/HimalayaUI/frontend/src/components/TagsTab.tsx`:

```tsx
export function TagsTab(): JSX.Element { return <p className="text-fg-muted italic">Tags.</p>; }
```

Create `packages/HimalayaUI/frontend/src/components/NotesTab.tsx`:

```tsx
export function NotesTab(): JSX.Element { return <p className="text-fg-muted italic">Notes.</p>; }
```

- [ ] **Step 5: Run to verify pass**

Run: `npm test -- PropertiesPanel && npm run build`
Expected: 4 PropertiesPanel tests pass; build clean.

- [ ] **Step 6: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PropertiesPanel.tsx \
        packages/HimalayaUI/frontend/src/components/ExposuresTab.tsx \
        packages/HimalayaUI/frontend/src/components/PeaksTab.tsx \
        packages/HimalayaUI/frontend/src/components/TagsTab.tsx \
        packages/HimalayaUI/frontend/src/components/NotesTab.tsx \
        packages/HimalayaUI/frontend/test/PropertiesPanel.test.tsx
git commit -m "feat(ui): PropertiesPanel shell with tabbed navigation"
```

---

## Task 3: ExposuresTab — reuse existing ExposureList

**Motivation:** `ExposureList` already renders the body we want. The tab is a thin shim so `PropertiesPanel` can treat all tabs uniformly as `() => JSX.Element`.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/ExposuresTab.tsx`

- [ ] **Step 1: Replace the stub**

Replace `packages/HimalayaUI/frontend/src/components/ExposuresTab.tsx`:

```tsx
import { ExposureList } from "./ExposureList";

export function ExposuresTab(): JSX.Element {
  return <ExposureList />;
}
```

No new tests — `ExposureList.test.tsx` already covers the underlying behavior.

- [ ] **Step 2: Run full suite**

Run: `npm test && npm run build`
Expected: all green.

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/ExposuresTab.tsx
git commit -m "feat(ui): ExposuresTab wraps ExposureList"
```

---

## Task 4: PeaksTab — table of q, prominence, sharpness, source

**Motivation:** Read-only view into the current exposure's peaks. Gives the user the numeric data behind the markers on the trace.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/PeaksTab.tsx`
- Create: `packages/HimalayaUI/frontend/test/PeaksTab.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/PeaksTab.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PeaksTab } from "../src/components/PeaksTab";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeExposureId: undefined });
});

function mockPeaks(items: unknown[]): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(items), {
      status: 200, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("<PeaksTab>", () => {
  it("shows a hint when no exposure is active", () => {
    renderWithProviders(<PeaksTab />);
    expect(screen.getByText(/no exposure selected/i)).toBeInTheDocument();
  });

  it("shows an empty-state message when the exposure has no peaks", async () => {
    useAppState.setState({ activeExposureId: 42 });
    mockPeaks([]);
    renderWithProviders(<PeaksTab />);
    await waitFor(() => expect(screen.getByText(/no peaks/i)).toBeInTheDocument());
  });

  it("renders one table row per peak with q, prominence, sharpness, source", async () => {
    useAppState.setState({ activeExposureId: 42 });
    mockPeaks([
      { id: 1, exposure_id: 42, q: 0.123, intensity: 10,
        prominence: 0.8, sharpness: 2.5, source: "auto" },
      { id: 2, exposure_id: 42, q: 0.456, intensity: 20,
        prominence: null, sharpness: null, source: "manual" },
    ]);
    renderWithProviders(<PeaksTab />);
    await waitFor(() => expect(screen.getByText("0.1230")).toBeInTheDocument());
    expect(screen.getByText("0.4560")).toBeInTheDocument();
    expect(screen.getByText("0.80")).toBeInTheDocument();
    expect(screen.getByText("2.50")).toBeInTheDocument();
    expect(screen.getByText("auto")).toBeInTheDocument();
    expect(screen.getByText("manual")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- PeaksTab`

- [ ] **Step 3: Implement**

Replace `packages/HimalayaUI/frontend/src/components/PeaksTab.tsx`:

```tsx
import { useAppState } from "../state";
import { usePeaks } from "../queries";

function num(value: number | null, digits: number): string {
  return value != null ? value.toFixed(digits) : "—";
}

export function PeaksTab(): JSX.Element {
  const activeExposureId = useAppState((s) => s.activeExposureId);
  const q = usePeaks(activeExposureId);

  if (activeExposureId === undefined) {
    return <p className="text-fg-muted italic">No exposure selected.</p>;
  }
  if (q.isPending) return <p className="text-fg-muted">Loading peaks…</p>;
  if (q.error) {
    return <p className="text-error">Error: {(q.error as Error).message}</p>;
  }

  const peaks = q.data ?? [];
  if (peaks.length === 0) {
    return <p className="text-fg-muted italic">No peaks for this exposure.</p>;
  }

  return (
    <table className="w-full text-[13px] font-mono">
      <thead className="text-left text-fg-muted uppercase text-[11px] tracking-wide">
        <tr>
          <th className="py-1 pr-3">q</th>
          <th className="py-1 pr-3">prominence</th>
          <th className="py-1 pr-3">sharpness</th>
          <th className="py-1">source</th>
        </tr>
      </thead>
      <tbody>
        {peaks.map((p) => (
          <tr key={p.id} data-peak-id={p.id} className="border-t border-border">
            <td className="py-1 pr-3">{num(p.q, 4)}</td>
            <td className="py-1 pr-3">{num(p.prominence, 2)}</td>
            <td className="py-1 pr-3">{num(p.sharpness, 2)}</td>
            <td className="py-1">{p.source}</td>
          </tr>
        ))}
      </tbody>
    </table>
  );
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- PeaksTab && npm run build`
Expected: 3/3 pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/PeaksTab.tsx \
        packages/HimalayaUI/frontend/test/PeaksTab.test.tsx
git commit -m "feat(ui): PeaksTab renders peaks table"
```

---

## Task 5: TagsTab — sample tag editor

**Motivation:** Let the user add and remove tags on the active sample. Existing tags live on the `Sample` object returned by `listSamples`; we derive the active sample and pass its tags to chip-style buttons with × to remove.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/TagsTab.tsx`
- Create: `packages/HimalayaUI/frontend/test/TagsTab.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/TagsTab.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { TagsTab } from "../src/components/TagsTab";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeSampleId: undefined });
});

function mockSamples(samples: unknown[]): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(samples), {
      status: 200, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("<TagsTab>", () => {
  it("shows a hint when no sample is active", () => {
    renderWithProviders(<TagsTab />);
    expect(screen.getByText(/no sample selected/i)).toBeInTheDocument();
  });

  it("renders existing tags as chips with a remove button each", async () => {
    useAppState.setState({ activeSampleId: 10 });
    mockSamples([{
      id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null,
      tags: [
        { id: 1, key: "lipid",   value: "DOPC",     source: "manifest" },
        { id: 2, key: "peptide", value: "melittin", source: "manual"   },
      ],
    }]);
    renderWithProviders(<TagsTab />);
    await waitFor(() =>
      expect(screen.getByText(/lipid/i)).toBeInTheDocument(),
    );
    expect(screen.getByText(/DOPC/i)).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /remove lipid/i })).toBeInTheDocument();
  });

  it("submitting the add form posts {key, value}", async () => {
    const user = userEvent.setup();
    useAppState.setState({ activeSampleId: 10 });
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url.endsWith("/api/experiments/1/samples")) {
        return new Response(JSON.stringify([{
          id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [],
        }]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url.endsWith("/api/samples/10/tags")) {
        return new Response(JSON.stringify({
          id: 5, sample_id: 10, key: "lipid", value: "DOPC", source: "manual",
        }), { status: 201, headers: { "Content-Type": "application/json" } });
      }
      return new Response("not found", { status: 404 });
    });

    renderWithProviders(<TagsTab />);
    await user.type(screen.getByPlaceholderText(/key/i), "lipid");
    await user.type(screen.getByPlaceholderText(/value/i), "DOPC");
    fireEvent.click(screen.getByRole("button", { name: /add tag/i }));

    await waitFor(() => {
      const urls = fetchSpy.mock.calls.map((c) =>
        typeof c[0] === "string" ? c[0] : (c[0] as Request).url);
      const postUrls = urls.filter((u) => u === "/api/samples/10/tags");
      expect(postUrls.length).toBeGreaterThan(0);
    });
  });

  it("clicking remove on an existing tag calls DELETE", async () => {
    useAppState.setState({ activeSampleId: 10 });
    const fetchSpy = vi.spyOn(global, "fetch").mockImplementation(async (input) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      const method = typeof input === "string" ? "GET" : (input as Request).method;
      if (url.endsWith("/api/experiments/1/samples") && method === "GET") {
        return new Response(JSON.stringify([{
          id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null,
          tags: [{ id: 1, key: "lipid", value: "DOPC", source: "manual" }],
        }]), { status: 200, headers: { "Content-Type": "application/json" } });
      }
      if (url.endsWith("/api/samples/10/tags/1")) {
        return new Response(null, { status: 204 });
      }
      return new Response("not found", { status: 404 });
    });

    renderWithProviders(<TagsTab />);
    const rm = await screen.findByRole("button", { name: /remove lipid/i });
    fireEvent.click(rm);
    await waitFor(() => {
      const calls = fetchSpy.mock.calls.filter((c) => {
        const url = typeof c[0] === "string" ? c[0] : (c[0] as Request).url;
        return url === "/api/samples/10/tags/1";
      });
      expect(calls.length).toBeGreaterThan(0);
    });
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- TagsTab`

- [ ] **Step 3: Implement `TagsTab.tsx`**

Replace `packages/HimalayaUI/frontend/src/components/TagsTab.tsx`:

```tsx
import { useState } from "react";
import { useAppState } from "../state";
import { useSamples, useAddSampleTag, useRemoveSampleTag } from "../queries";

const EXPERIMENT_ID = 1;

export function TagsTab(): JSX.Element {
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const samplesQ = useSamples(EXPERIMENT_ID);
  const addTag    = useAddSampleTag(EXPERIMENT_ID, activeSampleId ?? 0);
  const removeTag = useRemoveSampleTag(EXPERIMENT_ID, activeSampleId ?? 0);

  const [key, setKey]     = useState("");
  const [value, setValue] = useState("");

  if (activeSampleId === undefined) {
    return <p className="text-fg-muted italic">No sample selected.</p>;
  }
  if (samplesQ.isPending) return <p className="text-fg-muted">Loading…</p>;

  const sample = samplesQ.data?.find((s) => s.id === activeSampleId);
  const tags = sample?.tags ?? [];

  const inputClass =
    "bg-bg-elevated border border-border rounded-md px-2 py-1 focus:outline focus:outline-1 focus:outline-accent focus:border-accent";

  function onSubmit(ev: React.FormEvent): void {
    ev.preventDefault();
    const k = key.trim();
    const v = value.trim();
    if (!k || !v) return;
    addTag.mutate({ key: k, value: v }, {
      onSuccess: () => { setKey(""); setValue(""); },
    });
  }

  return (
    <div className="flex flex-col gap-3">
      <ul className="flex flex-wrap gap-1">
        {tags.length === 0 ? (
          <li className="text-fg-muted italic text-[13px]">No tags.</li>
        ) : tags.map((t) => (
          <li
            key={t.id}
            data-tag-id={t.id}
            className="flex items-center gap-1 px-2 py-1 bg-bg-elevated border border-border rounded-md text-[13px]"
          >
            <span className="font-mono">{t.key}</span>
            <span className="text-fg-muted">=</span>
            <span>{t.value}</span>
            {t.source === "manifest" && (
              <span className="text-fg-muted text-[11px] ml-1">(manifest)</span>
            )}
            <button
              className="ml-1 text-fg-muted hover:text-error rounded-md px-1 focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent"
              aria-label={`Remove ${t.key}`}
              onClick={() => removeTag.mutate(t.id)}
            >
              ×
            </button>
          </li>
        ))}
      </ul>

      <form onSubmit={onSubmit} className="flex gap-2 items-center">
        <input
          className={inputClass + " w-28"}
          type="text"
          placeholder="key"
          value={key}
          onChange={(e) => setKey(e.target.value)}
        />
        <input
          className={inputClass + " w-32"}
          type="text"
          placeholder="value"
          value={value}
          onChange={(e) => setValue(e.target.value)}
        />
        <button
          type="submit"
          className="bg-accent border border-accent text-white rounded-md px-2.5 py-1 hover:brightness-110 focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
          aria-label="Add tag"
          disabled={addTag.isPending}
        >
          {addTag.isPending ? "Adding…" : "Add"}
        </button>
      </form>
    </div>
  );
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- TagsTab && npm run build`

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/TagsTab.tsx \
        packages/HimalayaUI/frontend/test/TagsTab.test.tsx
git commit -m "feat(ui): TagsTab edits sample tags"
```

---

## Task 6: NotesTab — free-text editor with save-on-blur

**Motivation:** Let the user edit `sample.notes`. Save on blur if the value changed — simpler than debounced autosave and avoids spurious writes while typing.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/NotesTab.tsx`
- Create: `packages/HimalayaUI/frontend/test/NotesTab.test.tsx`

- [ ] **Step 1: Write failing tests**

Create `packages/HimalayaUI/frontend/test/NotesTab.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, waitFor } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { NotesTab } from "../src/components/NotesTab";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeSampleId: undefined });
});

function mockSample(notes: string | null): ReturnType<typeof vi.spyOn> {
  return vi.spyOn(global, "fetch").mockImplementation(async (input) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    if (url.endsWith("/api/experiments/1/samples")) {
      return new Response(JSON.stringify([{
        id: 10, experiment_id: 1, label: "A1", name: "s1", notes, tags: [],
      }]), { status: 200, headers: { "Content-Type": "application/json" } });
    }
    if (url.endsWith("/api/samples/10")) {
      return new Response(JSON.stringify({
        id: 10, experiment_id: 1, label: "A1", name: "s1",
        notes: "new notes", tags: [],
      }), { status: 200, headers: { "Content-Type": "application/json" } });
    }
    return new Response("not found", { status: 404 });
  });
}

describe("<NotesTab>", () => {
  it("shows a hint when no sample is active", () => {
    renderWithProviders(<NotesTab />);
    expect(screen.getByText(/no sample selected/i)).toBeInTheDocument();
  });

  it("pre-fills the textarea with existing notes", async () => {
    useAppState.setState({ activeSampleId: 10 });
    mockSample("existing value");
    renderWithProviders(<NotesTab />);
    const ta = await screen.findByRole("textbox");
    await waitFor(() => expect((ta as HTMLTextAreaElement).value).toBe("existing value"));
  });

  it("PATCHes /api/samples/:id on blur when value changed", async () => {
    useAppState.setState({ activeSampleId: 10 });
    const spy = mockSample("");
    renderWithProviders(<NotesTab />);
    const ta = await screen.findByRole("textbox") as HTMLTextAreaElement;
    fireEvent.change(ta, { target: { value: "new notes" } });
    fireEvent.blur(ta);
    await waitFor(() => {
      const patchCalls = spy.mock.calls.filter((c) => {
        const url = typeof c[0] === "string" ? c[0] : (c[0] as Request).url;
        const init = c[1] as RequestInit | undefined;
        return url.endsWith("/api/samples/10") && init?.method === "PATCH";
      });
      expect(patchCalls.length).toBeGreaterThan(0);
    });
  });

  it("does NOT PATCH on blur if the value is unchanged", async () => {
    useAppState.setState({ activeSampleId: 10 });
    const spy = mockSample("existing");
    renderWithProviders(<NotesTab />);
    const ta = await screen.findByRole("textbox") as HTMLTextAreaElement;
    await waitFor(() => expect(ta.value).toBe("existing"));
    fireEvent.blur(ta);
    // Give any pending network a tick; then assert no PATCH happened.
    await new Promise((r) => setTimeout(r, 10));
    const patches = spy.mock.calls.filter((c) =>
      (c[1] as RequestInit | undefined)?.method === "PATCH");
    expect(patches.length).toBe(0);
  });
});
```

- [ ] **Step 2: Run to verify failure**

Run: `npm test -- NotesTab`

- [ ] **Step 3: Implement `NotesTab.tsx`**

Replace `packages/HimalayaUI/frontend/src/components/NotesTab.tsx`:

```tsx
import { useEffect, useRef, useState } from "react";
import { useAppState } from "../state";
import { useSamples, useUpdateSample } from "../queries";

const EXPERIMENT_ID = 1;

export function NotesTab(): JSX.Element {
  const activeSampleId = useAppState((s) => s.activeSampleId);
  const samplesQ = useSamples(EXPERIMENT_ID);
  const updateSample = useUpdateSample(EXPERIMENT_ID, activeSampleId ?? 0);

  const sample = samplesQ.data?.find((s) => s.id === activeSampleId);
  const serverNotes = sample?.notes ?? "";

  const [value, setValue] = useState(serverNotes);
  const lastServerRef = useRef(serverNotes);

  // Sync local state when the server value changes (e.g., sample switch or
  // remote update). Avoid clobbering local edits the user hasn't blurred yet.
  useEffect(() => {
    if (serverNotes !== lastServerRef.current) {
      setValue(serverNotes);
      lastServerRef.current = serverNotes;
    }
  }, [serverNotes]);

  if (activeSampleId === undefined) {
    return <p className="text-fg-muted italic">No sample selected.</p>;
  }

  function onBlur(): void {
    if (value === lastServerRef.current) return;
    updateSample.mutate({ notes: value });
    lastServerRef.current = value;
  }

  return (
    <textarea
      aria-label="Sample notes"
      className="w-full h-full min-h-24 bg-bg-elevated border border-border rounded-md p-2 text-[13px] focus:outline focus:outline-1 focus:outline-accent focus:border-accent resize-none"
      value={value}
      onChange={(e) => setValue(e.target.value)}
      onBlur={onBlur}
      placeholder="Notes about this sample…"
    />
  );
}
```

- [ ] **Step 4: Run to verify pass**

Run: `npm test -- NotesTab && npm run build`
Expected: 4/4 pass.

- [ ] **Step 5: Commit**

```bash
git add packages/HimalayaUI/frontend/src/components/NotesTab.tsx \
        packages/HimalayaUI/frontend/test/NotesTab.test.tsx
git commit -m "feat(ui): NotesTab edits sample.notes on blur"
```

---

## Task 7: Wire PropertiesPanel into App.tsx

**Motivation:** Swap the bottom-center pane from `<ExposureList />` to `<PropertiesPanel />`. The tab content already includes an `ExposuresTab` wrapping `ExposureList`, so no UX loss.

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/App.tsx`

- [ ] **Step 1: Edit `App.tsx`**

In `packages/HimalayaUI/frontend/src/App.tsx`:

1. Remove the import line `import { ExposureList } from "./components/ExposureList";` (no longer referenced directly from App).
2. Add: `import { PropertiesPanel } from "./components/PropertiesPanel";`
3. Replace `centerBottom={<ExposureList />}` with `centerBottom={<PropertiesPanel />}`.

- [ ] **Step 2: Verify the full suite**

Run: `npm test && npm run build`
Expected: all green. `ExposureList.test.tsx` still runs (the component still exists, just called via ExposuresTab instead of App directly).

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/src/App.tsx
git commit -m "feat(ui): wire PropertiesPanel into bottom-center pane"
```

---

## Task 8: E2E — tab switching + add tag

**Motivation:** Lock in the new tabbed surface with one end-to-end flow: switch to Tags tab, add a tag, confirm the POST fires. Also verifies default tab = Exposures.

**Files:**
- Modify: `packages/HimalayaUI/frontend/e2e/smoke.spec.ts`

- [ ] **Step 1: Append test**

Append to `packages/HimalayaUI/frontend/e2e/smoke.spec.ts`:

```ts
test("switching to Tags tab and adding a tag posts to /api/samples/:id/tags", async ({ page }) => {
  let posted: { key: string; value: string } | null = null;

  await page.route("**/api/users", (r) => r.fulfill({ status: 200, body: "[]" }));
  await page.route("**/api/experiments/1", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify({
      id: 1, name: "demo", path: "/x", data_dir: "/x/data",
      analysis_dir: "/x/analysis", manifest_path: null,
      created_at: "2026-04-22T00:00:00Z",
    }),
  }));
  await page.route("**/api/experiments/1/samples", (r) => r.fulfill({
    status: 200,
    body: JSON.stringify([
      { id: 10, experiment_id: 1, label: "A1", name: "s1", notes: null, tags: [] },
    ]),
  }));
  await page.route("**/api/samples/10/exposures", (r) => r.fulfill({
    status: 200, body: "[]",
  }));
  await page.route("**/api/samples/10/tags", async (r) => {
    if (r.request().method() === "POST") {
      posted = r.request().postDataJSON() as { key: string; value: string };
      return r.fulfill({
        status: 201,
        body: JSON.stringify({
          id: 99, sample_id: 10,
          key: posted.key, value: posted.value, source: "manual",
        }),
      });
    }
    return r.fulfill({ status: 404, body: "nope" });
  });

  await page.addInitScript(() => {
    localStorage.setItem("himalaya-ui:state", JSON.stringify({
      state: { username: "alice", activeSampleId: 10, activeExposureId: undefined },
      version: 1,
    }));
  });

  await page.goto("/");
  // Default tab is Exposures
  await expect(page.getByRole("tab", { name: /exposures/i }))
    .toHaveAttribute("aria-selected", "true");

  await page.getByRole("tab", { name: /tags/i }).click();
  await expect(page.getByRole("tab", { name: /tags/i }))
    .toHaveAttribute("aria-selected", "true");

  await page.getByPlaceholder(/key/i).fill("lipid");
  await page.getByPlaceholder(/value/i).fill("DOPC");
  await page.getByRole("button", { name: /add tag/i }).click();

  await expect.poll(() => posted?.key, { timeout: 2000 }).toBe("lipid");
  await expect.poll(() => posted?.value, { timeout: 2000 }).toBe("DOPC");
});
```

- [ ] **Step 2: Build + E2E**

Run: `npm run build && npm run e2e`
Expected: all tests pass (existing 6 + new 1).

- [ ] **Step 3: Commit**

```bash
git add packages/HimalayaUI/frontend/e2e/smoke.spec.ts
git commit -m "test(e2e): tab switching + add tag"
```

---

## Final verification

From `packages/HimalayaUI/frontend/`:

```bash
npm test && npm run build && npm run e2e
```

Expected: all unit tests pass (67 existing + ~15 new ≈ 82), clean build, 7 E2E pass.

Plan 7+ ideas (deferred): Recent section in PhasePanel, export UI, per-user action audit view, beamline-config editor, derived-exposure construction UI.
