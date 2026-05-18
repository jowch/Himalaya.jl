# Loupe View Implementation Plan (#161 / I1.5)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the sample loupe view at `/samples/loupe/:sampleId` — a focused single-sample inspection surface with a full detector image, an exposure-thumbnail strip, and a metadata sidebar (per-exposure facts, drop/restore, set-representative, read-only sample tags).

**Architecture:** A URL-routed page (`LoupePage`) owns data fetching, the active-exposure state, keyboard handling, and the exposure mutations; two prop-driven presentational components (`LoupeFrame`, `LoupeSidebar`) render the surface. The page mounts under the existing `CorpusShell` layout route. It rebuilds the mockup layout reusing the leaf primitives `DetectorImage` and `ThumbnailGallery` — it does not re-compose the Inspect cards (which stay untouched until #163 deletes Inspect).

**Tech Stack:** React 18 + TypeScript (strict, `exactOptionalPropertyTypes`), react-router-dom, TanStack Query, Tailwind, `boneyard-js` skeletons, Vitest + React Testing Library (JSDOM).

**Design spec:** [`docs/superpowers/specs/2026-05-18-loupe-view-design.md`](../specs/2026-05-18-loupe-view-design.md). **Mockup:** `docs/redesign-mockups/sample-table.html` (`loupe-shell` markup).

**Working directory for all commands:** `packages/HimalayaUI/frontend/`.

---

## Conventions for every task

- **Test runner:** `node_modules/.bin/vitest run test/<file>` (one-shot). Tests live in `packages/HimalayaUI/frontend/test/`; they import source as `../src/...`.
- **Assert on `data-testid` / text, never Tailwind class strings** (`test/AGENTS.md`).
- **Skeleton gating:** always `query.isLoading`, never `isPending` (`docs/boneyard.md` Rule 1).
- **TypeScript:** optional props are modelled `T | null` or with explicit `undefined` in the union (`exactOptionalPropertyTypes` is on).
- **Palette:** the loupe mounts under `CorpusShell` — style it with the corpus "Print" palette (`paper`, `paper-sunk`, `ink`, `ink-soft`, `ink-faint`) plus the shared `accent` / `success` / `border` tokens, mirroring `CorpusShell.tsx` / `SamplesPage.tsx`. The Tailwind classes in each task are a working starting point; tests assert on `data-testid` / text, never class strings, so palette refinement never breaks a test.
- **Type-check:** run `npx tsc --noEmit` before each commit. The project has *pre-existing* errors in legacy test files that rely on vitest globals (`Cannot find name 'test'/'expect'`) — these are unrelated to this feature. Where a task step says `npx tsc --noEmit` → "Expected: clean", read it as **"no error line references a file this task creates or modifies"** (`LoupePage.tsx`, `LoupeFrame.tsx`, `LoupeSidebar.tsx`, `AppRoutes.tsx`, `LoupePage.test.tsx`, `LoupeFrame.test.tsx`, `LoupeSidebar.test.tsx`). The loupe test files use explicit `import { … } from "vitest"`, so they stay error-free. Task 8's `npm run build` (`tsc -p tsconfig.build.json`, src-scoped) is the real build gate.

## File structure

| File | New / Modify | Responsibility |
|---|---|---|
| `src/pages/LoupePage.tsx` | new | Route entry. Reads `:sampleId`; fetches the sample + exposures + experiment; owns active-exposure state, keyboard handling, back-nav, the drop/restore + set-representative mutations; renders loading / not-found; composes the plate. |
| `src/components/LoupeFrame.tsx` | new | Prop-driven. The big `DetectorImage size="full"` (with a "Dropped" overlay when rejected) + the exposure-thumbnail strip (`ThumbnailGallery`). |
| `src/components/LoupeSidebar.tsx` | new | Prop-driven. The "This exposure" meta-list, the Verdict box, the Representative box, the read-only Sample-tags section, the keyboard legend. |
| `src/components/AppRoutes.tsx` | modify | Add the `/samples/loupe/:sampleId` route slot under the existing `<CorpusShell>` group. |
| `src/bones/loupe.bones.json` | new (captured) | Boneyard skeleton fixture for the loupe body, captured via the dev-server HMR plugin. |
| `test/LoupePage.test.tsx` | new | Page-level tests (route param, not-found, flip, mutations, keyboard, skeleton, file-per-exposure guard). |
| `test/LoupeFrame.test.tsx` | new | `LoupeFrame` prop-driven tests. |
| `test/LoupeSidebar.test.tsx` | new | `LoupeSidebar` prop-driven tests. |

---

## Task 1: Route registration + LoupePage read scaffold

Render the loupe at `/samples/loupe/:sampleId`, showing the sample's identity in the head. Mutations and the exposure body come in later tasks.

**Files:**
- Create: `src/pages/LoupePage.tsx`
- Modify: `src/components/AppRoutes.tsx`
- Test: `test/LoupePage.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/LoupePage.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import type { CorpusSample, Exposure } from "../src/api";
import { LoupePage } from "../src/pages/LoupePage";

// LoupePage is the only module under test that imports from ../src/queries.
// Mutable holders let each test set query state before render.
const h = vi.hoisted(() => ({
  corpusQ: {} as { data?: CorpusSample[]; isLoading: boolean },
  exposuresQ: {} as { data?: Exposure[]; isLoading: boolean },
  experimentQ: {} as { data?: { id: number; name: string | null; path: string } },
  setStatusMutate: vi.fn(),
  setRepMutate: vi.fn(),
}));

vi.mock("../src/queries", () => ({
  useCorpusSamples: () => h.corpusQ,
  useExposures: () => h.exposuresQ,
  useExperiment: () => h.experimentQ,
  useSetExposureStatus: () => ({
    mutate: h.setStatusMutate, isPending: false, error: null, reset: () => {},
  }),
  useSelectExposure: () => ({
    mutate: h.setRepMutate, isPending: false, error: null, reset: () => {},
  }),
}));

function sample(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 7, experiment_id: 3, name: "JC042", display_name: "DOPE 80%",
    notes: null, tags: [], q_units: "A-1", ...over,
  };
}

function exposure(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 7, filename: "JC042-001.dat", kind: "file",
    selected: false, status: null, image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
    ...over,
  };
}

function renderAt(path: string) {
  const client = makeClient();
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route path="/samples" element={<div data-testid="samples-marker" />} />
          <Route path="/samples/loupe/:sampleId" element={<LoupePage />} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("LoupePage — identity", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.exposuresQ = { data: [exposure()], isLoading: false };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
  });

  it("renders the sample identified by the :sampleId route param", () => {
    renderAt("/samples/loupe/7");
    expect(screen.getByTestId("loupe-page")).toBeInTheDocument();
    expect(screen.getByText("DOPE 80%")).toBeInTheDocument();
    expect(screen.getByText(/Beamtime March/)).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: FAIL — `Failed to resolve import "../src/pages/LoupePage"`.

- [ ] **Step 3: Create the LoupePage read scaffold**

Create `src/pages/LoupePage.tsx`:

```tsx
import { useParams } from "react-router-dom";
import { useCorpusSamples, useExposures, useExperiment } from "../queries";

/**
 * LoupePage — the sample loupe at /samples/loupe/:sampleId. A focused
 * single-sample inspection surface: full detector image, exposure strip,
 * metadata sidebar. Mounted under the CorpusShell layout route (#161 / I1.5).
 *
 * URL-owned: the sample id comes from the route param, never from the
 * Zustand `activeSampleId` (master plan §2.3 — new surfaces own their URL).
 */
export function LoupePage(): JSX.Element {
  const { sampleId: sampleIdParam } = useParams<{ sampleId: string }>();
  const sampleId = Number(sampleIdParam);
  const hasValidId = Number.isFinite(sampleId);

  const corpusQ = useCorpusSamples();
  useExposures(hasValidId ? sampleId : undefined);

  const sample = corpusQ.data?.find((s) => s.id === sampleId);
  const experimentQ = useExperiment(sample?.experiment_id ?? 0);
  const experimentName =
    experimentQ.data?.name ?? experimentQ.data?.path ?? undefined;

  return (
    <div
      data-testid="loupe-page"
      className="mx-auto max-w-[1100px] px-8 py-7"
    >
      <header className="mb-5">
        <h2 className="text-2xl text-ink">
          {sample?.display_name ?? sample?.name ?? "—"}
        </h2>
        <p className="mt-1 font-mono text-xs text-ink-faint">
          {experimentName ?? "—"}
          {sample?.name ? ` · ${sample.name}` : ""}
        </p>
      </header>
    </div>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: PASS.

- [ ] **Step 5: Register the route**

In `src/components/AppRoutes.tsx`, add the import after the `SamplesPage` import (line 9):

```tsx
import { LoupePage } from "../pages/LoupePage";
```

Then add the loupe route inside the existing `<CorpusShell>` group — change:

```tsx
      <Route element={<CorpusShell />}>
        <Route path="/samples" element={<SamplesPage />} />
      </Route>
```

to:

```tsx
      <Route element={<CorpusShell />}>
        <Route path="/samples" element={<SamplesPage />} />
        <Route path="/samples/loupe/:sampleId" element={<LoupePage />} />
      </Route>
```

- [ ] **Step 6: Verify types and tests**

Run: `npx tsc --noEmit`
Expected: clean.
Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx test/AppRoutes.test.tsx`
Expected: PASS (the existing `AppRoutes` suite stays green — the add is append-only).

- [ ] **Step 7: Commit**

```bash
git add src/pages/LoupePage.tsx src/components/AppRoutes.tsx test/LoupePage.test.tsx
git commit -m "Add loupe route + LoupePage read scaffold (#161)"
```

---

## Task 2: LoupeFrame component

The big detector frame plus the exposure-thumbnail strip. Prop-driven — no hooks, no data fetching.

**Files:**
- Create: `src/components/LoupeFrame.tsx`
- Test: `test/LoupeFrame.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/LoupeFrame.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import type { Exposure } from "../src/api";
import { LoupeFrame } from "../src/components/LoupeFrame";

function exposure(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 7, filename: "JC042-001.dat", kind: "file",
    selected: false, status: null, image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
    ...over,
  };
}

describe("LoupeFrame", () => {
  it("renders one strip thumbnail per exposure", () => {
    const exposures = [
      exposure({ id: 1 }), exposure({ id: 2 }), exposure({ id: 3 }),
    ];
    render(
      <LoupeFrame
        exposure={exposures[0]}
        exposures={exposures}
        onSelectExposure={() => {}}
      />,
    );
    expect(screen.getByTestId("thumb-cell-1")).toBeInTheDocument();
    expect(screen.getByTestId("thumb-cell-2")).toBeInTheDocument();
    expect(screen.getByTestId("thumb-cell-3")).toBeInTheDocument();
  });

  it("calls onSelectExposure with the clicked exposure id", () => {
    const onSelect = vi.fn();
    const exposures = [exposure({ id: 1 }), exposure({ id: 2 })];
    render(
      <LoupeFrame
        exposure={exposures[0]}
        exposures={exposures}
        onSelectExposure={onSelect}
      />,
    );
    fireEvent.click(screen.getByTestId("thumb-cell-2"));
    expect(onSelect).toHaveBeenCalledWith(2);
  });

  it("shows the Dropped badge when the active exposure is rejected", () => {
    const e = exposure({ id: 1, status: "rejected" });
    const { rerender } = render(
      <LoupeFrame exposure={e} exposures={[e]} onSelectExposure={() => {}} />,
    );
    expect(screen.getByTestId("loupe-dropped-badge")).toBeInTheDocument();

    const kept = exposure({ id: 1, status: "accepted" });
    rerender(
      <LoupeFrame exposure={kept} exposures={[kept]} onSelectExposure={() => {}} />,
    );
    expect(screen.queryByTestId("loupe-dropped-badge")).not.toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `node_modules/.bin/vitest run test/LoupeFrame.test.tsx`
Expected: FAIL — `Failed to resolve import "../src/components/LoupeFrame"`.

- [ ] **Step 3: Create the component**

Create `src/components/LoupeFrame.tsx`:

```tsx
import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";
import { ThumbnailGallery } from "./ThumbnailGallery";

interface Props {
  /** The exposure shown in the big frame. */
  exposure: Exposure;
  /** All exposures for the sample — drives the strip. */
  exposures: Exposure[];
  onSelectExposure: (id: number) => void;
}

/**
 * LoupeFrame — the left column of the loupe: the full-size detector image
 * for the active exposure plus a thumbnail strip of every exposure.
 *
 * Iterates `exposures` by id and renders every one regardless of `kind`
 * (`file` / `averaged` / `background_subtracted`) — the loupe must not
 * assume one file per exposure (master plan §11; the deferred
 * derived-exposure feature must remain possible).
 */
export function LoupeFrame({
  exposure,
  exposures,
  onSelectExposure,
}: Props): JSX.Element {
  const isRejected = exposure.status === "rejected";
  const caption = exposure.filename ?? `Exposure #${exposure.id}`;

  return (
    <div data-testid="loupe-frame" className="flex flex-col gap-3">
      <div
        className={[
          "relative mx-auto aspect-square w-full max-w-[500px]",
          "overflow-hidden rounded border border-border bg-bg",
          isRejected ? "opacity-40" : "",
        ].join(" ")}
      >
        <DetectorImage
          exposureId={exposure.id}
          imagePath={exposure.image_path}
          imageVersion={exposure.image_version}
          size="full"
          className="h-full w-full"
        />
        {isRejected && (
          <span
            data-testid="loupe-dropped-badge"
            className="absolute left-3 top-3 rounded bg-accent px-2 py-0.5
                       text-[10px] font-bold uppercase tracking-wide text-bg"
          >
            Dropped
          </span>
        )}
        <span className="absolute bottom-2 left-3 font-mono text-[11px] text-ink-faint">
          {caption}
        </span>
      </div>
      <div className="h-[88px]">
        <ThumbnailGallery
          exposures={exposures}
          selectedId={exposure.id}
          onSelect={onSelectExposure}
          className="justify-center"
        />
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `node_modules/.bin/vitest run test/LoupeFrame.test.tsx`
Expected: PASS.

- [ ] **Step 5: Verify types**

Run: `npx tsc --noEmit`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add src/components/LoupeFrame.tsx test/LoupeFrame.test.tsx
git commit -m "Add LoupeFrame component (#161)"
```

---

## Task 3: LoupeSidebar component

The metadata sidebar: "This exposure" meta-list, Verdict box, Representative box, read-only Sample-tags section, keyboard legend. Prop-driven.

**Files:**
- Create: `src/components/LoupeSidebar.tsx`
- Test: `test/LoupeSidebar.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/LoupeSidebar.test.tsx`:

```tsx
import { describe, it, expect, vi } from "vitest";
import { render, screen, fireEvent } from "@testing-library/react";
import type { CorpusSample, Exposure } from "../src/api";
import { LoupeSidebar } from "../src/components/LoupeSidebar";

function sample(over: Partial<CorpusSample> = {}): CorpusSample {
  return {
    id: 7, experiment_id: 3, name: "JC042", display_name: "DOPE 80%",
    notes: null, tags: [], q_units: "A-1", ...over,
  };
}

function exposure(over: Partial<Exposure> = {}): Exposure {
  return {
    id: 100, sample_id: 7, filename: "JC042-001.dat", kind: "file",
    selected: false, status: null, image_path: null, image_version: "",
    tags: [], sources: [], trace_hash: null, analysis_inputs_hash: null,
    ...over,
  };
}

function defaultProps() {
  const e = exposure();
  return {
    exposure: e,
    exposures: [e],
    sample: sample(),
    onDropToggle: vi.fn(),
    onSetRepresentative: vi.fn(),
  };
}

describe("LoupeSidebar — meta-list", () => {
  it("shows filename, kind, frame position and status", () => {
    const exposures = [
      exposure({ id: 1 }), exposure({ id: 2, kind: "averaged", filename: null }),
    ];
    render(
      <LoupeSidebar
        {...defaultProps()}
        exposure={exposures[1]}
        exposures={exposures}
      />,
    );
    expect(screen.getByTestId("loupe-meta-filename")).toHaveTextContent("—");
    expect(screen.getByTestId("loupe-meta-kind")).toHaveTextContent("averaged");
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("2 of 2");
    expect(screen.getByTestId("loupe-meta-status")).toHaveTextContent("pending");
  });
});

describe("LoupeSidebar — verdict", () => {
  it("offers Drop for a kept exposure and calls onDropToggle", () => {
    const props = defaultProps();
    render(<LoupeSidebar {...props} />);
    const toggle = screen.getByTestId("loupe-drop-toggle");
    expect(toggle).toHaveTextContent("Drop");
    fireEvent.click(toggle);
    expect(props.onDropToggle).toHaveBeenCalledTimes(1);
  });

  it("offers Restore for a dropped exposure", () => {
    const props = defaultProps();
    render(<LoupeSidebar {...props} exposure={exposure({ status: "rejected" })} />);
    expect(screen.getByTestId("loupe-drop-toggle")).toHaveTextContent("Restore");
  });

  it("displays an existing rejection reason for a dropped exposure", () => {
    const props = defaultProps();
    render(
      <LoupeSidebar
        {...props}
        exposure={exposure({
          status: "rejected",
          tags: [
            { id: 5, key: "rejection_reason", value: "Beam flare", source: "manual" },
          ],
        })}
      />,
    );
    expect(screen.getByTestId("loupe-verdict")).toHaveTextContent("Beam flare");
  });
});

describe("LoupeSidebar — representative", () => {
  it("offers Set as representative and calls onSetRepresentative", () => {
    const props = defaultProps();
    render(<LoupeSidebar {...props} />);
    const btn = screen.getByTestId("loupe-set-representative");
    fireEvent.click(btn);
    expect(props.onSetRepresentative).toHaveBeenCalledTimes(1);
  });

  it("marks an already-representative exposure and hides the button", () => {
    const props = defaultProps();
    render(<LoupeSidebar {...props} exposure={exposure({ selected: true })} />);
    expect(screen.getByTestId("loupe-rep")).toHaveTextContent(/Representative/);
    expect(screen.queryByTestId("loupe-set-representative")).not.toBeInTheDocument();
  });
});

describe("LoupeSidebar — sample tags", () => {
  it("renders tags as read-only chips with no remove control", () => {
    const props = defaultProps();
    render(
      <LoupeSidebar
        {...props}
        sample={sample({
          tags: [{ id: 9, key: "lipid", value: "DOPE", source: "manual" }],
        })}
      />,
    );
    const tags = screen.getByTestId("loupe-tags");
    expect(tags).toHaveTextContent("lipid: DOPE");
    // Read-only: SampleMetadataCard's remove button is aria-labelled
    // "Remove <key> tag" — the loupe must not render it (editing is #159).
    expect(screen.queryByLabelText("Remove lipid tag")).not.toBeInTheDocument();
  });

  it("shows an empty-state hint when the sample has no tags", () => {
    render(<LoupeSidebar {...defaultProps()} />);
    expect(screen.getByTestId("loupe-tags")).toHaveTextContent("No tags yet");
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `node_modules/.bin/vitest run test/LoupeSidebar.test.tsx`
Expected: FAIL — `Failed to resolve import "../src/components/LoupeSidebar"`.

- [ ] **Step 3: Create the component**

Create `src/components/LoupeSidebar.tsx`:

```tsx
import type { CorpusSample, Exposure } from "../api";

interface Props {
  /** The active exposure. */
  exposure: Exposure;
  /** All exposures for the sample — drives the "frame N of M" position. */
  exposures: Exposure[];
  /** The sample, for the read-only tags section. */
  sample: CorpusSample;
  /** Drop the exposure if kept, restore it if dropped. */
  onDropToggle: () => void;
  /** Mark the active exposure as the indexing representative. */
  onSetRepresentative: () => void;
}

function MetaRow(
  { label, value, testid }: { label: string; value: string; testid: string },
): JSX.Element {
  return (
    <div className="flex justify-between font-mono text-[11.5px]">
      <span className="text-ink-faint">{label}</span>
      <span data-testid={testid} className="text-ink">{value}</span>
    </div>
  );
}

function SectionHeading({ children }: { children: string }): JSX.Element {
  return (
    <div className="mb-2 text-[10px] font-bold uppercase tracking-wide text-ink-faint">
      {children}
    </div>
  );
}

/**
 * LoupeSidebar — the right column of the loupe. Per-exposure facts, a
 * keep/drop verdict, the indexing-representative control, a read-only
 * sample-tags section, and a keyboard legend.
 *
 * The tags section is intentionally read-only: the corpus sample-tag
 * add/remove round-trip is #159 (I1.3), out of scope for #161.
 */
export function LoupeSidebar({
  exposure,
  exposures,
  sample,
  onDropToggle,
  onSetRepresentative,
}: Props): JSX.Element {
  const isRejected = exposure.status === "rejected";
  const isRepresentative = exposure.selected;
  const frameIndex = exposures.findIndex((e) => e.id === exposure.id);
  const frameLabel =
    frameIndex >= 0 ? `${frameIndex + 1} of ${exposures.length}` : "—";
  // Display an existing rejection reason if one was authored elsewhere (the
  // Inspect card / a future culling surface). The loupe never *authors* a
  // reason — its drop is a plain status toggle (spec §6).
  const rejectionReason = exposure.tags.find(
    (t) => t.key === "rejection_reason",
  )?.value;

  return (
    <aside data-testid="loupe-sidebar" className="flex flex-col gap-5">
      {/* This exposure */}
      <section>
        <SectionHeading>This exposure</SectionHeading>
        <div className="flex flex-col gap-1.5">
          <MetaRow
            label="Filename"
            value={exposure.filename ?? "—"}
            testid="loupe-meta-filename"
          />
          <MetaRow label="Kind" value={exposure.kind} testid="loupe-meta-kind" />
          <MetaRow label="Frame" value={frameLabel} testid="loupe-meta-frame" />
          <MetaRow
            label="Status"
            value={exposure.status ?? "pending"}
            testid="loupe-meta-status"
          />
        </div>
      </section>

      {/* Verdict */}
      <section
        data-testid="loupe-verdict"
        className="flex items-center gap-3 rounded border border-border bg-bg-subtle p-3"
      >
        <span
          className={[
            "h-2.5 w-2.5 shrink-0 rounded-full",
            isRejected ? "bg-accent" : "bg-success",
          ].join(" ")}
        />
        <div className="flex-1">
          <div className="text-[13px] font-bold text-ink">
            {isRejected ? "Dropped" : "Kept"}
          </div>
          <div className="text-[10.5px] text-ink-faint">
            {isRejected
              ? (rejectionReason ?? "Dropped from this sample.")
              : "Everything is kept until you drop it."}
          </div>
        </div>
        <button
          data-testid="loupe-drop-toggle"
          onClick={onDropToggle}
          className="rounded border border-border bg-paper px-2.5 py-1.5
                     text-[11.5px] font-semibold text-ink hover:bg-bg-subtle"
        >
          {isRejected ? "Restore" : "Drop"}
        </button>
      </section>

      {/* Representative */}
      <section
        data-testid="loupe-rep"
        className="rounded border border-border p-3"
      >
        {isRepresentative ? (
          <div className="flex items-center gap-2 text-xs font-bold text-accent">
            <span className="h-2 w-2 rounded-full bg-accent" />
            Representative for indexing
          </div>
        ) : (
          <>
            <div className="text-[11.5px] text-ink-soft">
              One exposure per sample carries forward to the Index stage.
              Pick the cleanest, strongest frame.
            </div>
            <button
              data-testid="loupe-set-representative"
              onClick={onSetRepresentative}
              className="mt-2 rounded border border-border bg-paper px-2.5 py-1.5
                         text-[11.5px] font-semibold text-ink hover:bg-bg-subtle"
            >
              Set as representative
            </button>
          </>
        )}
      </section>

      {/* Sample tags — read-only (editing is #159) */}
      <section>
        <SectionHeading>Sample tags</SectionHeading>
        <div data-testid="loupe-tags" className="flex flex-wrap gap-1">
          {sample.tags.length === 0 ? (
            <span className="text-[11.5px] text-ink-faint">No tags yet</span>
          ) : (
            sample.tags.map((tag) => (
              <span
                key={tag.id}
                className="inline-flex items-center rounded-full border border-border
                           bg-bg-subtle px-2 py-0.5 text-xs text-ink-soft"
              >
                {tag.key}: {tag.value}
              </span>
            ))
          )}
        </div>
      </section>

      {/* Keyboard legend */}
      <section>
        <SectionHeading>Keys</SectionHeading>
        <div className="flex flex-col gap-1 text-[11px] text-ink-faint">
          <div><kbd className="font-mono">← →</kbd> flip frames</div>
          <div><kbd className="font-mono">X</kbd> drop / restore</div>
          <div><kbd className="font-mono">R</kbd> set representative</div>
          <div><kbd className="font-mono">Esc</kbd> back to the sheet</div>
        </div>
      </section>
    </aside>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `node_modules/.bin/vitest run test/LoupeSidebar.test.tsx`
Expected: PASS.

- [ ] **Step 5: Verify types**

Run: `npx tsc --noEmit`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add src/components/LoupeSidebar.tsx test/LoupeSidebar.test.tsx
git commit -m "Add LoupeSidebar component (#161)"
```

---

## Task 4: LoupePage composition — exposure body, default selection, not-found

Wire the active-exposure state, the default-exposure selection rule, the not-found panel, the strip-flip, the back button — composing `LoupeFrame` and `LoupeSidebar`.

**Files:**
- Modify: `src/pages/LoupePage.tsx`
- Test: `test/LoupePage.test.tsx`

- [ ] **Step 1: Write the failing tests**

First, update the React Testing Library import at the top of `test/LoupePage.test.tsx` — Tasks 4–7 need `fireEvent`. Change:

```tsx
import { render, screen } from "@testing-library/react";
```

to:

```tsx
import { render, screen, fireEvent } from "@testing-library/react";
```

Then append a new `describe` block to `test/LoupePage.test.tsx`, after the existing one:

```tsx
describe("LoupePage — composition", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.exposuresQ = {
      data: [
        exposure({ id: 100, status: "accepted" }),
        exposure({ id: 101, status: "rejected" }),
      ],
      isLoading: false,
    };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
  });

  it("renders the frame and sidebar for a known sample", () => {
    renderAt("/samples/loupe/7");
    expect(screen.getByTestId("loupe-frame")).toBeInTheDocument();
    expect(screen.getByTestId("loupe-sidebar")).toBeInTheDocument();
  });

  it("default-selects the first accepted exposure", () => {
    renderAt("/samples/loupe/7");
    // Exposure 100 is accepted, 101 rejected → 100 is the default.
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("1 of 2");
    expect(screen.getByTestId("loupe-meta-status")).toHaveTextContent("accepted");
  });

  it("flips the active exposure when a strip thumbnail is clicked", () => {
    renderAt("/samples/loupe/7");
    fireEvent.click(screen.getByTestId("thumb-cell-101"));
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("2 of 2");
    expect(screen.getByTestId("loupe-meta-status")).toHaveTextContent("rejected");
  });

  it("shows a not-found panel for a sample id absent from the corpus", () => {
    renderAt("/samples/loupe/999");
    expect(screen.getByTestId("loupe-not-found")).toBeInTheDocument();
    expect(screen.queryByTestId("loupe-frame")).not.toBeInTheDocument();
  });

  it("navigates back to /samples when the back button is clicked", () => {
    renderAt("/samples/loupe/7");
    fireEvent.click(screen.getByTestId("loupe-back"));
    expect(screen.getByTestId("samples-marker")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: FAIL — `loupe-frame` / `loupe-not-found` / `loupe-back` not found.

- [ ] **Step 3: Rewrite LoupePage with the exposure body**

Replace the entire contents of `src/pages/LoupePage.tsx` with:

```tsx
import { useEffect, useMemo, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { useCorpusSamples, useExposures, useExperiment } from "../queries";
import { LoupeFrame } from "../components/LoupeFrame";
import { LoupeSidebar } from "../components/LoupeSidebar";
import type { Exposure } from "../api";

/**
 * Default exposure when the loupe opens: the indexing representative, else
 * the first accepted exposure, else the first exposure. Lifted from
 * InspectPage so the loupe opens on the same frame Inspect would have.
 */
function defaultExposureId(exposures: Exposure[]): number | undefined {
  const representative = exposures.find((e) => e.selected);
  if (representative) return representative.id;
  const firstAccepted = exposures.find((e) => e.status === "accepted");
  if (firstAccepted) return firstAccepted.id;
  return exposures[0]?.id;
}

/**
 * LoupePage — the sample loupe at /samples/loupe/:sampleId. A focused
 * single-sample inspection surface: full detector image, exposure strip,
 * metadata sidebar. Mounted under the CorpusShell layout route (#161 / I1.5).
 *
 * URL-owned: the sample id comes from the route param, never from the
 * Zustand `activeSampleId` (master plan §2.3 — new surfaces own their URL).
 */
export function LoupePage(): JSX.Element {
  const { sampleId: sampleIdParam } = useParams<{ sampleId: string }>();
  const navigate = useNavigate();
  const sampleId = Number(sampleIdParam);
  const hasValidId = Number.isFinite(sampleId);

  const corpusQ = useCorpusSamples();
  const exposuresQ = useExposures(hasValidId ? sampleId : undefined);

  const sample = corpusQ.data?.find((s) => s.id === sampleId);
  const experimentQ = useExperiment(sample?.experiment_id ?? 0);
  const experimentName =
    experimentQ.data?.name ?? experimentQ.data?.path ?? undefined;

  const exposures = useMemo(() => exposuresQ.data ?? [], [exposuresQ.data]);

  // Active exposure — local state, defaulted by `defaultExposureId`. Reset
  // when the sample changes so the next sample picks its own default.
  const [activeId, setActiveId] = useState<number | undefined>(undefined);
  useEffect(() => {
    setActiveId(undefined);
  }, [sampleId]);
  const computedDefault = defaultExposureId(exposures);
  useEffect(() => {
    if (activeId === undefined && computedDefault !== undefined) {
      setActiveId(computedDefault);
    }
  }, [activeId, computedDefault]);
  const activeExposure = exposures.find((e) => e.id === activeId);

  function goBack(): void {
    navigate("/samples");
  }

  if (!corpusQ.isLoading && !sample) {
    return (
      <div data-testid="loupe-page" className="mx-auto max-w-[1100px] px-8 py-7">
        <div
          data-testid="loupe-not-found"
          className="rounded border border-border p-8 text-sm text-ink-faint"
        >
          Sample not found.{" "}
          <button
            onClick={goBack}
            className="font-semibold text-accent hover:underline"
          >
            Back to the sheet
          </button>
        </div>
      </div>
    );
  }

  return (
    <div data-testid="loupe-page" className="mx-auto max-w-[1100px] px-8 py-7">
      <button
        data-testid="loupe-back"
        onClick={goBack}
        className="mb-3.5 text-[11.5px] font-semibold text-accent hover:underline"
      >
        ← Back to the sheet
      </button>
      <header className="mb-5">
        <h2 className="text-2xl text-ink">
          {sample?.display_name ?? sample?.name ?? "—"}
        </h2>
        <p className="mt-1 font-mono text-xs text-ink-faint">
          {experimentName ?? "—"}
          {sample?.name ? ` · ${sample.name}` : ""}
        </p>
      </header>

      <div className="grid grid-cols-[1fr_286px] gap-7">
        {sample && activeExposure ? (
          <>
            <LoupeFrame
              exposure={activeExposure}
              exposures={exposures}
              onSelectExposure={setActiveId}
            />
            <LoupeSidebar
              exposure={activeExposure}
              exposures={exposures}
              sample={sample}
              onDropToggle={() => {}}
              onSetRepresentative={() => {}}
            />
          </>
        ) : (
          <div className="col-span-2 p-8 text-sm text-ink-faint">
            This sample has no exposures.
          </div>
        )}
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: PASS (all tests, including Task 1's identity test).

- [ ] **Step 5: Verify types**

Run: `npx tsc --noEmit`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add src/pages/LoupePage.tsx test/LoupePage.test.tsx
git commit -m "Compose LoupeFrame + LoupeSidebar into LoupePage (#161)"
```

---

## Task 5: Interactions — drop/restore, set-representative, keyboard

Wire the exposure mutations to the sidebar controls and add the keyboard shortcuts (arrows flip, X drops/restores, R sets representative, Esc goes back).

**Files:**
- Modify: `src/pages/LoupePage.tsx`
- Test: `test/LoupePage.test.tsx`

Background — the mutation hook call shapes (from `src/queries.ts`):
- `useSetExposureStatus(sampleId)` → `.mutate({ exposureId, status })` where `status` is `"accepted" | "rejected" | null`.
- `useSelectExposure(sampleId)` → `.mutate(exposureId)` (a bare number).

- [ ] **Step 1: Write the failing tests**

Append to `test/LoupePage.test.tsx` (a new `describe` block):

```tsx
describe("LoupePage — interactions", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.exposuresQ = {
      data: [
        exposure({ id: 100, status: "accepted" }),
        exposure({ id: 101, status: null }),
      ],
      isLoading: false,
    };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
  });

  it("drops a kept exposure via the verdict toggle", () => {
    renderAt("/samples/loupe/7"); // default active = 100 (accepted)
    fireEvent.click(screen.getByTestId("loupe-drop-toggle"));
    expect(h.setStatusMutate).toHaveBeenCalledWith({
      exposureId: 100, status: "rejected",
    });
  });

  it("restores a dropped exposure via the verdict toggle", () => {
    h.exposuresQ = {
      data: [exposure({ id: 100, status: "rejected" })],
      isLoading: false,
    };
    renderAt("/samples/loupe/7");
    fireEvent.click(screen.getByTestId("loupe-drop-toggle"));
    expect(h.setStatusMutate).toHaveBeenCalledWith({
      exposureId: 100, status: null,
    });
  });

  it("sets the representative via the rep button", () => {
    renderAt("/samples/loupe/7");
    fireEvent.click(screen.getByTestId("loupe-set-representative"));
    expect(h.setRepMutate).toHaveBeenCalledWith(100);
  });

  it("flips frames with the arrow keys", () => {
    renderAt("/samples/loupe/7"); // default active = 100
    fireEvent.keyDown(document.body, { key: "ArrowRight" });
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("2 of 2");
    fireEvent.keyDown(document.body, { key: "ArrowLeft" });
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("1 of 2");
  });

  it("drops the active exposure with the X key", () => {
    renderAt("/samples/loupe/7");
    fireEvent.keyDown(document.body, { key: "x" });
    expect(h.setStatusMutate).toHaveBeenCalledWith({
      exposureId: 100, status: "rejected",
    });
  });

  it("sets the representative with the R key", () => {
    renderAt("/samples/loupe/7");
    fireEvent.keyDown(document.body, { key: "r" });
    expect(h.setRepMutate).toHaveBeenCalledWith(100);
  });

  it("goes back to /samples with the Escape key", () => {
    renderAt("/samples/loupe/7");
    fireEvent.keyDown(document.body, { key: "Escape" });
    expect(screen.getByTestId("samples-marker")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run the tests to verify they fail**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: FAIL — `h.setStatusMutate` / `h.setRepMutate` not called; arrow keys do nothing.

- [ ] **Step 3: Wire the mutations and keyboard handler**

In `src/pages/LoupePage.tsx`, update the imports — change:

```tsx
import { useEffect, useMemo, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import { useCorpusSamples, useExposures, useExperiment } from "../queries";
```

to:

```tsx
import { useCallback, useEffect, useMemo, useState } from "react";
import { useNavigate, useParams } from "react-router-dom";
import {
  useCorpusSamples,
  useExposures,
  useExperiment,
  useSetExposureStatus,
  useSelectExposure,
} from "../queries";
```

**Hook-ordering rule — load-bearing.** Every hook added in this task — `useSetExposureStatus`, `useSelectExposure`, the three `useCallback`s, the `useCallback` form of `goBack`, and the keyboard `useEffect` — must appear **above** the `if (!corpusQ.isLoading && !sample) return (…)` early return. The not-found branch returns before the JSX; if any hook is placed after it, React throws "rendered more hooks than during the previous render" on the not-found path. The insertion points below all sit above that early return — keep them there.

Then, inside the component, after the `activeExposure` line:

```tsx
  const activeExposure = exposures.find((e) => e.id === activeId);
```

add the mutation hooks and handlers:

```tsx
  const setStatus = useSetExposureStatus(hasValidId ? sampleId : 0);
  const setRepresentative = useSelectExposure(hasValidId ? sampleId : 0);

  const handleDropToggle = useCallback(() => {
    if (!activeExposure) return;
    setStatus.mutate({
      exposureId: activeExposure.id,
      status: activeExposure.status === "rejected" ? null : "rejected",
    });
  }, [activeExposure, setStatus]);

  const handleSetRepresentative = useCallback(() => {
    if (!activeExposure) return;
    setRepresentative.mutate(activeExposure.id);
  }, [activeExposure, setRepresentative]);

  const flip = useCallback(
    (delta: number) => {
      if (activeId === undefined || exposures.length === 0) return;
      const idx = exposures.findIndex((e) => e.id === activeId);
      if (idx < 0) return;
      const next = Math.min(Math.max(idx + delta, 0), exposures.length - 1);
      setActiveId(exposures[next].id);
    },
    [activeId, exposures],
  );
```

Then change `goBack` from a function declaration to a `useCallback` (so the keyboard effect can depend on it) — replace:

```tsx
  function goBack(): void {
    navigate("/samples");
  }
```

with:

```tsx
  const goBack = useCallback(() => {
    navigate("/samples");
  }, [navigate]);
```

Then add the keyboard effect immediately after `goBack`:

```tsx
  // Loupe keyboard shortcuts. The input/textarea guard is in from the start
  // so it survives #159 adding tag-edit inputs to the sidebar.
  useEffect(() => {
    function onKeyDown(e: KeyboardEvent): void {
      const tag = (e.target as HTMLElement | null)?.tagName;
      if (tag === "INPUT" || tag === "TEXTAREA") return;
      if (e.key === "ArrowLeft") flip(-1);
      else if (e.key === "ArrowRight") flip(1);
      else if (e.key === "x" || e.key === "X") handleDropToggle();
      else if (e.key === "r" || e.key === "R") handleSetRepresentative();
      else if (e.key === "Escape") goBack();
    }
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [flip, handleDropToggle, handleSetRepresentative, goBack]);
```

**Known benign overlap with `useGlobalShortcuts`.** `useGlobalShortcuts` (mounted in `AppRoutes`, app-wide) also binds `ArrowLeft`/`ArrowRight` — to legacy page-tab switching via `setActivePage`. On a corpus route both `window` listeners fire, so an arrow press flips the loupe exposure *and* mutates the Zustand `activePage`. This is harmless: `setActivePage` is a pure state write, and the URL-sync hook (`useUrlFromState`) is mounted only inside the legacy `AppShell` — it does **not** run on `/samples/loupe/*`, so no navigation occurs. The `activePage` drift is coerced by the Phase 1 nav-bridge. No guard is needed in #161; scoping `useGlobalShortcuts` to legacy routes is a dual-nav-retirement concern (Phase 5), not this issue.

Finally, wire the two handlers into `LoupeSidebar` — change:

```tsx
            <LoupeSidebar
              exposure={activeExposure}
              exposures={exposures}
              sample={sample}
              onDropToggle={() => {}}
              onSetRepresentative={() => {}}
            />
```

to:

```tsx
            <LoupeSidebar
              exposure={activeExposure}
              exposures={exposures}
              sample={sample}
              onDropToggle={handleDropToggle}
              onSetRepresentative={handleSetRepresentative}
            />
```

- [ ] **Step 4: Run the tests to verify they pass**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: PASS (all describe blocks).

- [ ] **Step 5: Verify types**

Run: `npx tsc --noEmit`
Expected: clean.

- [ ] **Step 6: Commit**

```bash
git add src/pages/LoupePage.tsx test/LoupePage.test.tsx
git commit -m "Wire loupe drop/restore, representative, keyboard (#161)"
```

---

## Task 6: Boneyard skeleton

Show a skeleton while the sample or exposures query is in its cold-load state.

**Files:**
- Modify: `src/pages/LoupePage.tsx`
- Create (captured): `src/bones/loupe.bones.json`
- Test: `test/LoupePage.test.tsx`

Background — `<Skeleton>` falls through to `fallback` when no bones are committed for `name` (`docs/boneyard.md` Rule 3). The unit test asserts the fallback appears; the captured `loupe.bones.json` is produced by the dev-server HMR plugin in Step 6.

- [ ] **Step 1: Write the failing test**

Append to `test/LoupePage.test.tsx` (a new `describe` block):

```tsx
describe("LoupePage — loading", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: undefined, isLoading: true };
    h.exposuresQ = { data: undefined, isLoading: true };
    h.experimentQ = { data: undefined };
  });

  it("shows the loupe skeleton while data is loading", () => {
    renderAt("/samples/loupe/7");
    expect(screen.getByTestId("loupe-skeleton")).toBeInTheDocument();
    // Body content must not render while loading.
    expect(screen.queryByTestId("loupe-frame")).not.toBeInTheDocument();
    expect(screen.queryByTestId("loupe-not-found")).not.toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: FAIL — `loupe-skeleton` not found.

- [ ] **Step 3: Add the skeleton wrapper**

In `src/pages/LoupePage.tsx`, add the `Skeleton` import after the react-router import:

```tsx
import { Skeleton } from "boneyard-js/react";
```

Add a module-level loading flag and a skeleton fixture. Place this `LOUPE_FIXTURE` constant above the `LoupePage` function, after the imports:

The fixture's exposure has `image_path: null`, so `DetectorImage` takes its placeholder branch — the captured detector-frame bone is a plain rectangle sized by `LoupeFrame`'s `aspect-square` wrapper rather than a measured image. That is the correct skeleton geometry (the frame is a fixed-aspect box regardless of image content), so leaving `image_path: null` is deliberate.

```tsx
// Fixture for boneyard's headless capture — a real render with mock props so
// the capture CLI can measure the loupe body layout (docs/boneyard.md Rule 2).
const LOUPE_FIXTURE_EXPOSURE = {
  id: 0, sample_id: 0, filename: "JC000-001.dat", kind: "file" as const,
  selected: false, status: "accepted" as const, image_path: null,
  image_version: "", tags: [], sources: [],
  trace_hash: null, analysis_inputs_hash: null,
};
const LOUPE_FIXTURE = (
  <div className="grid grid-cols-[1fr_286px] gap-7">
    <LoupeFrame
      exposure={LOUPE_FIXTURE_EXPOSURE}
      exposures={[LOUPE_FIXTURE_EXPOSURE]}
      onSelectExposure={() => {}}
    />
    <LoupeSidebar
      exposure={LOUPE_FIXTURE_EXPOSURE}
      exposures={[LOUPE_FIXTURE_EXPOSURE]}
      sample={{
        id: 0, experiment_id: 0, name: "JC000", display_name: "Sample",
        notes: null, tags: [], q_units: "A-1",
      }}
      onDropToggle={() => {}}
      onSetRepresentative={() => {}}
    />
  </div>
);
```

Inside the component, compute the loading flag — add after the `exposures` `useMemo`:

```tsx
  const isLoading = corpusQ.isLoading || exposuresQ.isLoading;
```

Then wrap the body grid. Replace the `<div className="grid grid-cols-[1fr_286px] gap-7">…</div>` block with a `<Skeleton>`-wrapped version:

```tsx
      <Skeleton
        name="loupe"
        className="block"
        loading={isLoading}
        stagger={50}
        transition={200}
        fixture={LOUPE_FIXTURE}
        fallback={
          <div
            data-testid="loupe-skeleton"
            className="p-8 text-sm italic text-ink-faint"
          >
            Loading sample…
          </div>
        }
      >
        <div className="grid grid-cols-[1fr_286px] gap-7">
          {sample && activeExposure ? (
            <>
              <LoupeFrame
                exposure={activeExposure}
                exposures={exposures}
                onSelectExposure={setActiveId}
              />
              <LoupeSidebar
                exposure={activeExposure}
                exposures={exposures}
                sample={sample}
                onDropToggle={handleDropToggle}
                onSetRepresentative={handleSetRepresentative}
              />
            </>
          ) : (
            <div className="col-span-2 p-8 text-sm text-ink-faint">
              This sample has no exposures.
            </div>
          )}
        </div>
      </Skeleton>
```

Note: the `data-testid="loupe-skeleton"` lives on the `fallback`. When committed bones exist, `<Skeleton>` renders the captured geometry instead — the test asserts the fallback because bones are not committed at unit-test time, which is the documented behaviour (`docs/boneyard.md` Rule 3).

- [ ] **Step 4: Run the test to verify it passes**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: PASS (all describe blocks).

- [ ] **Step 5: Verify types**

Run: `npx tsc --noEmit`
Expected: clean.

- [ ] **Step 6: Capture the bones**

Start the dev server and trigger the loupe's loading state so the boneyard HMR plugin captures geometry (`docs/boneyard.md` "Adding a new skeleton"):

```bash
npm run dev
```

In the browser, open a loupe URL (e.g. `/samples/loupe/<an existing sample id>`) and reload so the cold-load path renders. The Vite plugin writes `src/bones/loupe.bones.json` and regenerates `src/bones/registry.ts`. Stop the dev server. Verify both files changed:

```bash
git status --short src/bones/
```

Expected: `src/bones/loupe.bones.json` (new) and `src/bones/registry.ts` (modified) appear.

If the dev server is not available in this environment, skip the capture: `<Skeleton>` falls through to the `fallback` (the documented no-bones behaviour) and the acceptance criterion "a boneyard skeleton shows while loading" is still met by the wrapper + fallback. Note the skip in the commit message.

- [ ] **Step 7: Commit**

```bash
git add src/pages/LoupePage.tsx test/LoupePage.test.tsx src/bones/
git commit -m "Add boneyard skeleton to the loupe (#161)"
```

---

## Task 7: File-per-exposure regression guard

The loupe must not assume one file per exposure (master plan §11; the deferred derived-exposure feature). Add an explicit regression test exercising a mixed-`kind` exposure list, including an exposure with `filename: null`.

**Files:**
- Test: `test/LoupePage.test.tsx`

- [ ] **Step 1: Write the regression test**

Append to `test/LoupePage.test.tsx` (a new `describe` block):

```tsx
describe("LoupePage — not file-per-exposure", () => {
  beforeEach(() => {
    vi.clearAllMocks();
    h.corpusQ = { data: [sample()], isLoading: false };
    h.experimentQ = { data: { id: 3, name: "Beamtime March", path: "/x" } };
  });

  it("renders and flips every exposure regardless of kind or missing filename", () => {
    // A file exposure, an averaged exposure with no filename, and a
    // background-subtracted exposure — the loupe keys off exposure id, never
    // off a non-null filename.
    h.exposuresQ = {
      data: [
        exposure({ id: 200, kind: "file", filename: "JC042-001.dat" }),
        exposure({ id: 201, kind: "averaged", filename: null }),
        exposure({ id: 202, kind: "background_subtracted", filename: null }),
      ],
      isLoading: false,
    };
    renderAt("/samples/loupe/7");

    // All three render in the strip.
    expect(screen.getByTestId("thumb-cell-200")).toBeInTheDocument();
    expect(screen.getByTestId("thumb-cell-201")).toBeInTheDocument();
    expect(screen.getByTestId("thumb-cell-202")).toBeInTheDocument();

    // The filename-less averaged exposure is flippable and renders cleanly.
    fireEvent.click(screen.getByTestId("thumb-cell-201"));
    expect(screen.getByTestId("loupe-meta-kind")).toHaveTextContent("averaged");
    expect(screen.getByTestId("loupe-meta-filename")).toHaveTextContent("—");
    expect(screen.getByTestId("loupe-meta-frame")).toHaveTextContent("2 of 3");

    // And the background-subtracted one too.
    fireEvent.click(screen.getByTestId("thumb-cell-202"));
    expect(screen.getByTestId("loupe-meta-kind"))
      .toHaveTextContent("background_subtracted");
  });
});
```

- [ ] **Step 2: Run the test**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx`
Expected: PASS — the loupe was built to iterate exposures by id (Tasks 2 & 4), so this guard passes without code changes. If it fails, the failure points at a `filename`-non-null assumption that must be removed.

- [ ] **Step 3: Run the full loupe suite + types**

Run: `node_modules/.bin/vitest run test/LoupePage.test.tsx test/LoupeFrame.test.tsx test/LoupeSidebar.test.tsx`
Expected: PASS.
Run: `npx tsc --noEmit`
Expected: clean.

- [ ] **Step 4: Commit**

```bash
git add test/LoupePage.test.tsx
git commit -m "Add file-per-exposure regression guard for the loupe (#161)"
```

---

## Task 8: Full build verification

**Files:** none (verification only).

- [ ] **Step 1: Run the full frontend build**

Run (from `packages/HimalayaUI/frontend/`): `npm run build`
Expected: `tsc --noEmit` and `vite build` both succeed — no type errors, no build errors.

- [ ] **Step 2: Run the full Vitest suite**

Run: `npm test > /tmp/vitest.out 2>&1; grep -E "FAIL|passed|failed" /tmp/vitest.out`
Expected: no `FAIL` lines; the loupe suites and all pre-existing suites pass. The `AppRoutes` suite in particular must stay green (the route add is append-only).

- [ ] **Step 3: Commit (only if a fix was needed)**

If Steps 1–2 surfaced a regression, fix it, then:

```bash
git add -A
git commit -m "Fix loupe build/test regressions (#161)"
```

If Steps 1–2 were clean, no commit is needed — the feature is complete.

---

## Acceptance criteria check (issue #161)

| Criterion | Task |
|---|---|
| `/samples/loupe/:sampleId` renders a full detector image, an exposure strip, and a metadata sidebar | 1, 2, 3, 4 |
| The metadata sidebar includes a Sample-tags section | 3 |
| The loupe flips between exposure frames | 4 (strip click), 5 (arrow keys) |
| The loupe does not hard-code a file-per-exposure assumption | 2, 4, 7 (regression guard) |
| A boneyard skeleton shows while loading | 6 |
| Vitest covers the loupe | 2, 3, 4, 5, 6, 7 |

## Out of scope (deferred to sibling issues)

- Corpus sample-tag add/remove round-trip — #159 (I1.3). The loupe's tags section is read-only (Task 3).
- The contact-sheet table and inbound links to the loupe — #160 (I1.4).
- Contact-sheet bulk culling — #162 (I1.6).
- Inspect deletion — #163 (I1.7).
- Sample `display_name` / `notes` editing and "index this exposure" navigation — #179 (I4.2, Phase 4).
