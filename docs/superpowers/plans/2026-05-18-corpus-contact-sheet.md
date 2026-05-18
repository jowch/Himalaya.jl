# Corpus Contact-Sheet Sample Table Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build the corpus contact-sheet sample table at `/samples` (#160 / I1.4) — a corpus-wide table where each row shows a sample's identity, exposure-thumbnail strip, Kept count, Tags, and Status, with a `?beamtime=` experiment filter.

**Architecture:** A `SamplesPage` container calls `useCorpusSamples()` once, filters client-side by a `?beamtime=` URL query, and renders one `ContactSheetRow` per sample inside a boneyard skeleton. Each `ContactSheetRow` owns its own `useExposures(sampleId)` query (per-sample fan-out) — the exposure strip and Kept count derive from it. The `CorpusTopbar` Beamtime chip becomes a functional experiment picker that writes `?beamtime=` to the URL.

**Tech Stack:** React 18 + TypeScript, TanStack Query v5, react-router-dom v6, Tailwind (`@theme` "The Print" palette), boneyard-js skeletons, Vitest + Testing Library.

---

## Background — read before starting

- **Spec:** `docs/superpowers/specs/2026-05-18-corpus-contact-sheet-design.md` — the design this plan implements. Read it.
- **Frontend conventions:** `packages/HimalayaUI/frontend/src/AGENTS.md` and `.../src/components/AGENTS.md`.
- All frontend commands run from `packages/HimalayaUI/frontend/`.
- **`npm test` is one-shot** (a settings hook injects vitest's `--run` flag). `npm test -- <substring>` runs only matching files.
- **Verified facts the code below depends on:**
  - `useCorpusSamples()` (`src/queries.ts`) returns `CorpusSample[]`; `CorpusSample` (`src/api.ts:40`) = `{ id, experiment_id, name, display_name, notes, tags, q_units }`. No exposures.
  - `useExposures(sampleId)` (`src/queries.ts:145`) returns `Exposure[]`; `Exposure` (`src/api.ts:129`) has `id, sample_id, filename, kind, selected, status, image_path, image_version, tags, sources, trace_hash, analysis_inputs_hash`.
  - `useExperiments()` (`src/queries.ts:106`) returns `Experiment[]`; `Experiment` (`src/api.ts:10`) has `id, name, path, …, q_units`.
  - `SampleTag` (`src/api.ts:21`) = `{ id, key, value, source }`.
  - `DetectorImage` (`src/components/DetectorImage.tsx`) props: `exposureId`, `imagePath`, `imageVersion`, `size: "thumb" | "full"`, `className?`.
  - The boneyard `<Skeleton>` is imported from `boneyard-js/react`; gate on `query.isLoading` (boneyard.md Rule 1).
  - `data-testid="samples-page"` on the `SamplesPage` root is asserted by `test/AppRoutes.test.tsx` and `test/CorpusShell.test.tsx` — **it must survive the rewrite.**
  - "The Print" palette tokens (`src/styles.css`): `paper`, `paper-sunk`, `plate`, `ink`, `ink-soft`, `ink-faint`, `hair`, `hair-strong`, `print-accent` → Tailwind classes `bg-paper`, `text-ink`, `border-hair`, `text-print-accent`, etc.

---

## File Structure

| File | Responsibility |
|---|---|
| `src/components/ContactSheetRow.tsx` | **Create.** One sample row: 5 cells, owns its `useExposures` query. Exports the `SampleStatus` type and the `CONTACT_SHEET_COLS` grid-template class string. |
| `src/pages/SamplesPage.tsx` | **Rewrite.** Contact-sheet container: corpus query, `?beamtime=` filter, column header, skeleton, row list. |
| `src/components/CorpusTopbar.tsx` | **Modify.** Turn the inert Beamtime chip into an experiment-picker `<select>` writing `?beamtime=`. |
| `test/contact-sheet.test.tsx` | **Create.** Vitest for `ContactSheetRow`, `SamplesPage`, and the `?beamtime=` round-trip. |
| `test/CorpusTopbar.test.tsx` | **Modify.** Add QueryClient/fetch-mock wrapper; cover the experiment picker. |
| `test/CorpusShell.test.tsx` | **Modify.** Add QueryClient + fetch mock — `SamplesPage` now issues real queries. |

---

## Task 1: `ContactSheetRow` component

**Files:**
- Create: `packages/HimalayaUI/frontend/src/components/ContactSheetRow.tsx`
- Create: `packages/HimalayaUI/frontend/test/contact-sheet.test.tsx`

- [ ] **Step 1: Write the failing test**

Create `test/contact-sheet.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { makeClient } from "./test-utils";
import { ContactSheetRow } from "../src/components/ContactSheetRow";
import type { CorpusSample, Exposure } from "../src/api";

/** Route fetch by path so per-sample exposure fan-out is order-independent. */
function mockFetch(routes: Record<string, unknown>): void {
  vi.spyOn(global, "fetch").mockImplementation((input: RequestInfo | URL) => {
    const url = typeof input === "string" ? input : (input as Request).url;
    const path = url.split("?")[0];
    if (!(path in routes)) {
      return Promise.resolve(new Response(`unmocked: ${url}`, { status: 404 }));
    }
    return Promise.resolve(
      new Response(JSON.stringify(routes[path]), {
        status: 200,
        headers: { "Content-Type": "application/json" },
      }),
    );
  });
}

function makeExposure(
  over: Partial<Exposure> & { id: number; sample_id: number },
): Exposure {
  return {
    filename: `f${over.id}.dat`,
    kind: "file",
    selected: false,
    status: null,
    image_path: null,
    image_version: "",
    tags: [],
    sources: [],
    trace_hash: null,
    analysis_inputs_hash: null,
    ...over,
  };
}

function makeSample(over: Partial<CorpusSample> & { id: number }): CorpusSample {
  return {
    experiment_id: 1,
    name: `sample ${over.id}`,
    display_name: null,
    notes: null,
    tags: [],
    q_units: "A-1",
    ...over,
  };
}

function renderRow(sample: CorpusSample) {
  const client = makeClient();
  const wrapper = ({ children }: { children: ReactNode }) => (
    <QueryClientProvider client={client}>{children}</QueryClientProvider>
  );
  return render(<ContactSheetRow sample={sample} />, { wrapper });
}

describe("ContactSheetRow", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("renders the sample identity (name and id)", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(makeSample({ id: 7, name: "Lipid 1-2" }));
    const cell = await screen.findByTestId("sample-cell");
    expect(cell).toHaveTextContent("Lipid 1-2");
    expect(cell).toHaveTextContent("#7");
  });

  it("prefers display_name over name", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(makeSample({ id: 7, name: "raw", display_name: "Pretty Name" }));
    expect(await screen.findByTestId("sample-cell")).toHaveTextContent(
      "Pretty Name",
    );
  });

  it("renders one thumbnail per exposure", async () => {
    mockFetch({
      "/api/samples/7/exposures": [
        makeExposure({ id: 1, sample_id: 7 }),
        makeExposure({ id: 2, sample_id: 7 }),
        makeExposure({ id: 3, sample_id: 7 }),
      ],
    });
    renderRow(makeSample({ id: 7 }));
    await waitFor(() => {
      expect(screen.getByTestId("exposure-thumb-1")).toBeInTheDocument();
      expect(screen.getByTestId("exposure-thumb-2")).toBeInTheDocument();
      expect(screen.getByTestId("exposure-thumb-3")).toBeInTheDocument();
    });
  });

  it("marks a rejected exposure thumbnail", async () => {
    mockFetch({
      "/api/samples/7/exposures": [
        makeExposure({ id: 1, sample_id: 7, status: "rejected" }),
      ],
    });
    renderRow(makeSample({ id: 7 }));
    await waitFor(() => {
      expect(screen.getByTestId("exposure-thumb-1")).toHaveAttribute(
        "data-rejected",
        "true",
      );
    });
  });

  it("shows kept / total with a dropped sub-label", async () => {
    mockFetch({
      "/api/samples/7/exposures": [
        makeExposure({ id: 1, sample_id: 7, status: null }),
        makeExposure({ id: 2, sample_id: 7, status: "accepted" }),
        makeExposure({ id: 3, sample_id: 7, status: "rejected" }),
      ],
    });
    renderRow(makeSample({ id: 7 }));
    const kept = await screen.findByTestId("kept-cell");
    await waitFor(() => {
      expect(kept).toHaveTextContent("2");
      expect(kept).toHaveTextContent("3");
      expect(kept).toHaveTextContent("1 dropped");
    });
  });

  it("omits the dropped sub-label when nothing is dropped", async () => {
    mockFetch({
      "/api/samples/7/exposures": [makeExposure({ id: 1, sample_id: 7 })],
    });
    renderRow(makeSample({ id: 7 }));
    const kept = await screen.findByTestId("kept-cell");
    await waitFor(() => expect(kept).toHaveTextContent("1"));
    expect(kept).not.toHaveTextContent("dropped");
  });

  it("renders sample tags as chips", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(
      makeSample({
        id: 7,
        tags: [
          { id: 1, key: "lipid", value: "DOPC", source: "manual" },
          { id: 2, key: "", value: "LL37", source: "manual" },
        ],
      }),
    );
    const tags = await screen.findByTestId("tags-cell");
    expect(tags).toHaveTextContent("DOPC");
    expect(tags).toHaveTextContent("LL37");
  });

  it("renders an inert (disabled) tag-add button", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(makeSample({ id: 7 }));
    expect(await screen.findByTestId("tag-add")).toBeDisabled();
  });

  it("renders a 'Not indexed' status", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(makeSample({ id: 7 }));
    expect(await screen.findByTestId("status-cell")).toHaveTextContent(
      "Not indexed",
    );
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npm test -- contact-sheet`
Expected: FAIL — `Failed to resolve import "../src/components/ContactSheetRow"`.

- [ ] **Step 3: Create the `ContactSheetRow` component**

Create `src/components/ContactSheetRow.tsx`:

```tsx
import { useExposures } from "../queries";
import type { CorpusSample, Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";

/**
 * A sample's index/phase status column. #160 ships only "not-indexed";
 * a later issue wires the real phase call into this typed seam.
 */
export type SampleStatus = "not-indexed";

/**
 * Shared CSS grid template for the contact sheet — the column header in
 * SamplesPage and every ContactSheetRow use this so the columns align.
 */
export const CONTACT_SHEET_COLS =
  "grid grid-cols-[16rem_1fr_7rem_14rem_8rem] gap-4 items-center";

interface Props {
  sample: CorpusSample;
}

/** One exposure thumbnail — inert in #160 (culling wiring is #162). */
function ExposureThumb({ exposure }: { exposure: Exposure }): JSX.Element {
  const isRejected = exposure.status === "rejected";
  const isRepresentative = exposure.selected;
  return (
    <div
      data-testid={`exposure-thumb-${exposure.id}`}
      data-rejected={isRejected ? "true" : undefined}
      data-representative={isRepresentative ? "true" : undefined}
      className={[
        "relative w-12 shrink-0 aspect-[3/4] overflow-hidden rounded",
        "ring-1 ring-hair",
        isRejected ? "opacity-40 grayscale" : "",
      ].join(" ")}
    >
      <DetectorImage
        exposureId={exposure.id}
        imagePath={exposure.image_path}
        imageVersion={exposure.image_version}
        size="thumb"
        className="h-full w-full"
      />
      {isRepresentative && (
        <span
          className="absolute left-0.5 top-0.5 text-[10px] text-print-accent"
          title="representative exposure"
        >
          ⊙
        </span>
      )}
      {isRejected && (
        <span className="absolute inset-0 flex items-center justify-center text-print-accent">
          ✕
        </span>
      )}
    </div>
  );
}

/**
 * ContactSheetRow — one sample row of the corpus contact sheet (#160).
 *
 * Owns its own useExposures query (per-sample fan-out) so the table fills
 * in row-by-row. The same queryKeys.exposures(sampleId) cache entry is
 * reused by culling (#162) and the loupe (#161).
 *
 * Inert affordances: the thumbnails carry no onClick and the tag-add
 * button is disabled — selection (#162) and tag mutation (#159) wire in
 * separately.
 */
export function ContactSheetRow({ sample }: Props): JSX.Element {
  const exposuresQuery = useExposures(sample.id);
  const exposures = exposuresQuery.data ?? [];

  const total = exposures.length;
  const kept = exposures.filter((e) => e.status !== "rejected").length;
  const dropped = total - kept;

  const name = sample.display_name ?? sample.name ?? `#${sample.id}`;

  return (
    <div
      data-testid={`sample-row-${sample.id}`}
      className={`${CONTACT_SHEET_COLS} border-b border-hair px-4 py-3`}
    >
      {/* Sample — identity only (no screened mark; that is #162). */}
      <div data-testid="sample-cell" className="flex flex-col">
        <span className="font-semibold text-ink">{name}</span>
        <span className="text-xs text-ink-faint">#{sample.id}</span>
      </div>

      {/* Exposures — thumbnail strip. */}
      <div
        data-testid="exposures-cell"
        className="flex h-16 flex-row gap-2 overflow-x-auto"
      >
        {exposuresQuery.isLoading ? (
          <span className="self-center text-xs text-ink-faint">
            Loading frames…
          </span>
        ) : (
          exposures.map((e) => <ExposureThumb key={e.id} exposure={e} />)
        )}
      </div>

      {/* Kept — kept / total, plus an "N dropped" sub-label. */}
      <div data-testid="kept-cell" className="flex flex-col text-sm">
        {exposuresQuery.isLoading ? (
          <span className="text-ink-faint">—</span>
        ) : (
          <>
            <span className="text-ink">
              {kept}
              <span className="text-ink-faint"> / {total}</span>
            </span>
            {dropped > 0 && (
              <span className="text-xs text-print-accent">
                {dropped} dropped
              </span>
            )}
          </>
        )}
      </div>

      {/* Tags — read-only chips + inert add button (mutation is #159). */}
      <div data-testid="tags-cell" className="flex flex-wrap items-center gap-1">
        {sample.tags.map((t) => (
          <span
            key={t.id}
            title={t.key || undefined}
            className="rounded bg-paper-sunk px-1.5 py-0.5 text-xs text-ink-soft"
          >
            {t.value}
          </span>
        ))}
        <button
          type="button"
          data-testid="tag-add"
          disabled
          title="Add a tag (coming soon)"
          className="rounded border border-dashed border-hair-strong px-1.5
                     py-0.5 text-xs text-ink-faint"
        >
          + tag
        </button>
      </div>

      {/* Status — fixed placeholder behind the SampleStatus seam.
          TODO: wire the real phase call when an issue is scoped for it. */}
      <div data-testid="status-cell">
        <span className="text-xs text-ink-faint">Not indexed</span>
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npm test -- contact-sheet`
Expected: PASS — all 9 `ContactSheetRow` tests green.

- [ ] **Step 5: Commit**

```bash
git add src/components/ContactSheetRow.tsx test/contact-sheet.test.tsx
git commit -m "Add ContactSheetRow — one corpus contact-sheet row (#160)"
```

---

## Task 2: `SamplesPage` container

**Files:**
- Modify (rewrite): `packages/HimalayaUI/frontend/src/pages/SamplesPage.tsx`
- Modify: `packages/HimalayaUI/frontend/test/contact-sheet.test.tsx` (append)
- Modify: `packages/HimalayaUI/frontend/test/CorpusShell.test.tsx`

- [ ] **Step 1: Write the failing test**

Append to `test/contact-sheet.test.tsx` (the `mockFetch` / `makeSample` / `makeExposure` helpers from Task 1 are already in the file — reuse them; add the new imports at the top of the file):

Add to the import block at the top of the file:

```tsx
import { MemoryRouter } from "react-router-dom";
import { SamplesPage } from "../src/pages/SamplesPage";
```

Append these describe blocks at the end of the file:

```tsx
const EXPERIMENTS = [
  { id: 1, name: "SSRL Apr 2026", path: "/e1", data_dir: "/d1",
    analysis_dir: "/a1", manifest_path: null, created_at: "2026-04-01",
    q_units: "A-1" },
  { id: 2, name: "APS Jul 2026", path: "/e2", data_dir: "/d2",
    analysis_dir: "/a2", manifest_path: null, created_at: "2026-07-01",
    q_units: "A-1" },
];

/** Corpus of 3 samples: ids 10,11 in experiment 1; id 12 in experiment 2. */
const CORPUS = [
  makeSample({ id: 10, experiment_id: 1, name: "alpha" }),
  makeSample({ id: 11, experiment_id: 1, name: "beta" }),
  makeSample({ id: 12, experiment_id: 2, name: "gamma" }),
];

/** Route map covering the corpus query, experiments, and per-row fan-out. */
function corpusRoutes(): Record<string, unknown> {
  return {
    "/api/samples": CORPUS,
    "/api/experiments": EXPERIMENTS,
    "/api/samples/10/exposures": [makeExposure({ id: 1, sample_id: 10 })],
    "/api/samples/11/exposures": [makeExposure({ id: 2, sample_id: 11 })],
    "/api/samples/12/exposures": [makeExposure({ id: 3, sample_id: 12 })],
  };
}

function renderSamplesPage(initialPath = "/samples") {
  const client = makeClient();
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={[initialPath]}>
        <SamplesPage />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("SamplesPage", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("keeps the samples-page test id", async () => {
    mockFetch(corpusRoutes());
    renderSamplesPage();
    expect(await screen.findByTestId("samples-page")).toBeInTheDocument();
  });

  it("renders one row per corpus sample", async () => {
    mockFetch(corpusRoutes());
    renderSamplesPage();
    await waitFor(() => {
      expect(screen.getByTestId("sample-row-10")).toBeInTheDocument();
      expect(screen.getByTestId("sample-row-11")).toBeInTheDocument();
      expect(screen.getByTestId("sample-row-12")).toBeInTheDocument();
    });
  });

  it("shows the boneyard skeleton while the corpus query loads", async () => {
    // A fetch that never resolves keeps the query in isLoading.
    vi.spyOn(global, "fetch").mockReturnValue(new Promise(() => {}));
    renderSamplesPage();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
    expect(screen.queryByTestId("contact-sheet-rows")).toBeNull();
  });

  it("filters to one experiment when ?beamtime= is set", async () => {
    mockFetch(corpusRoutes());
    renderSamplesPage("/samples?beamtime=1");
    await waitFor(() =>
      expect(screen.getByTestId("sample-row-10")).toBeInTheDocument(),
    );
    expect(screen.getByTestId("sample-row-11")).toBeInTheDocument();
    expect(screen.queryByTestId("sample-row-12")).toBeNull();
  });

  it("names the active experiment in the header when filtered", async () => {
    mockFetch(corpusRoutes());
    renderSamplesPage("/samples?beamtime=2");
    expect(await screen.findByTestId("samples-scope")).toHaveTextContent(
      "APS Jul 2026",
    );
  });

  it("shows an error state when the corpus query fails", async () => {
    mockFetch({}); // every path → 404
    renderSamplesPage();
    expect(await screen.findByTestId("samples-error")).toBeInTheDocument();
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npm test -- contact-sheet`
Expected: FAIL — the new `SamplesPage` tests fail (the current placeholder renders neither rows, `samples-scope`, nor `samples-error`).

- [ ] **Step 3: Rewrite `SamplesPage`**

Replace the entire contents of `src/pages/SamplesPage.tsx`:

```tsx
import { useSearchParams } from "react-router-dom";
import { Skeleton } from "boneyard-js/react";
import { useCorpusSamples, useExperiments } from "../queries";
import {
  ContactSheetRow,
  CONTACT_SHEET_COLS,
} from "../components/ContactSheetRow";

/** Static skeleton shape for the boneyard headless capture (boneyard.md). */
const CONTACT_SHEET_FIXTURE = (
  <div>
    {[0, 1, 2, 3].map((i) => (
      <div
        key={i}
        className={`${CONTACT_SHEET_COLS} border-b border-hair px-4 py-3`}
      >
        <div className="h-8 rounded bg-paper-sunk" />
        <div className="h-16 rounded bg-paper-sunk" />
        <div className="h-8 rounded bg-paper-sunk" />
        <div className="h-8 rounded bg-paper-sunk" />
        <div className="h-8 rounded bg-paper-sunk" />
      </div>
    ))}
  </div>
);

/**
 * SamplesPage — the corpus contact-sheet table at /samples (#160 / I1.4).
 *
 * Owns the corpus query, the ?beamtime= URL filter, and the boneyard
 * skeleton. ?beamtime=<experiment_id> is read here (the topbar chip is the
 * writer — CorpusTopbar); the corpus is filtered client-side because
 * useCorpusSamples() fetches the whole corpus and exposes no filter.
 */
export function SamplesPage(): JSX.Element {
  const [searchParams] = useSearchParams();
  const beamtimeRaw = searchParams.get("beamtime");
  const beamtime =
    beamtimeRaw !== null && /^\d+$/.test(beamtimeRaw)
      ? Number(beamtimeRaw)
      : undefined;

  const corpusQuery = useCorpusSamples();
  const experimentsQuery = useExperiments();

  const samples = corpusQuery.data ?? [];
  const filtered =
    beamtime === undefined
      ? samples
      : samples.filter((s) => s.experiment_id === beamtime);

  const scopeName =
    beamtime === undefined
      ? "the whole corpus"
      : (experimentsQuery.data?.find((e) => e.id === beamtime)?.name ??
        `experiment ${beamtime}`);

  return (
    <div data-testid="samples-page" className="flex flex-col gap-4 p-6">
      <header className="flex flex-col gap-1">
        <div className="text-xs font-semibold uppercase tracking-wide text-print-accent">
          Contact sheet
        </div>
        <p data-testid="samples-scope" className="text-sm text-ink-faint">
          Showing {scopeName}.
        </p>
      </header>

      <div
        className={`${CONTACT_SHEET_COLS} border-b border-hair-strong px-4 pb-2
                    text-xs font-semibold uppercase tracking-wide text-ink-faint`}
      >
        <div>Sample</div>
        <div>Exposures</div>
        <div>Kept</div>
        <div>Tags</div>
        <div>Status</div>
      </div>

      {corpusQuery.isError ? (
        <div data-testid="samples-error" className="px-4 py-8 text-sm text-ink-soft">
          Could not load samples. Try reloading the page.
        </div>
      ) : (
        <Skeleton
          name="contact-sheet"
          className="w-full"
          loading={corpusQuery.isLoading}
          stagger={50}
          transition={200}
          fixture={CONTACT_SHEET_FIXTURE}
          fallback={
            <div data-testid="contact-sheet-skeleton" className="px-4 py-8" />
          }
        >
          <div data-testid="contact-sheet-rows">
            {filtered.map((s) => (
              <ContactSheetRow key={s.id} sample={s} />
            ))}
          </div>
        </Skeleton>
      )}

      <div className="flex gap-4 px-4 text-xs text-ink-faint">
        <span>click — select a frame</span>
        <span>X — drop the selected frames</span>
        <span>double-click — open the loupe</span>
        <span>Esc — clear</span>
      </div>
    </div>
  );
}
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npm test -- contact-sheet`
Expected: PASS — all `ContactSheetRow` and `SamplesPage` tests green.

- [ ] **Step 5: Fix `CorpusShell.test.tsx` for the now-querying `SamplesPage`**

`SamplesPage` now calls `useCorpusSamples()` / `useExperiments()`, which throw without a `QueryClientProvider`. Replace the entire contents of `test/CorpusShell.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { CorpusShell } from "../src/components/CorpusShell";
import { SamplesPage } from "../src/pages/SamplesPage";

/** SamplesPage and CorpusTopbar both issue queries — answer them all 200/[]. */
function mockEmptyApi(): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify([]), {
      status: 200,
      headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("CorpusShell", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("renders the topbar and the matched child route in its Outlet", () => {
    mockEmptyApi();
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/samples"]}>
          <Routes>
            <Route element={<CorpusShell />}>
              <Route path="/samples" element={<SamplesPage />} />
            </Route>
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );
    expect(screen.getByTestId("corpus-shell")).toBeInTheDocument();
    expect(screen.getByTestId("corpus-topbar")).toBeInTheDocument();
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
  });
});

describe("SamplesPage", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("renders the contact-sheet body", () => {
    mockEmptyApi();
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/samples"]}>
          <SamplesPage />
        </MemoryRouter>
      </QueryClientProvider>,
    );
    expect(screen.getByTestId("samples-page")).toBeInTheDocument();
  });
});
```

- [ ] **Step 6: Run the affected suites to verify they pass**

Run: `npm test -- contact-sheet CorpusShell`
Expected: PASS — both files green.

- [ ] **Step 7: Commit**

```bash
git add src/pages/SamplesPage.tsx test/contact-sheet.test.tsx test/CorpusShell.test.tsx
git commit -m "Build the SamplesPage corpus contact-sheet table (#160)"
```

---

## Task 3: Functional Beamtime chip in `CorpusTopbar`

**Files:**
- Modify: `packages/HimalayaUI/frontend/src/components/CorpusTopbar.tsx`
- Modify: `packages/HimalayaUI/frontend/test/CorpusTopbar.test.tsx`

- [ ] **Step 1: Write the failing test**

Replace the entire contents of `test/CorpusTopbar.test.tsx`:

```tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, useSearchParams } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { CorpusTopbar } from "../src/components/CorpusTopbar";

const EXPERIMENTS = [
  { id: 1, name: "SSRL Apr 2026", path: "/e1", data_dir: "/d1",
    analysis_dir: "/a1", manifest_path: null, created_at: "2026-04-01",
    q_units: "A-1" },
  { id: 2, name: "APS Jul 2026", path: "/e2", data_dir: "/d2",
    analysis_dir: "/a2", manifest_path: null, created_at: "2026-07-01",
    q_units: "A-1" },
];

function mockExperiments(): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(EXPERIMENTS), {
      status: 200,
      headers: { "Content-Type": "application/json" },
    }),
  );
}

/** Surfaces the live ?beamtime= value so tests can assert URL writes. */
function BeamtimeProbe(): JSX.Element {
  const [params] = useSearchParams();
  return <span data-testid="beamtime-probe">{params.get("beamtime") ?? ""}</span>;
}

function renderTopbar(initialPath = "/samples") {
  return render(
    <QueryClientProvider client={makeClient()}>
      <MemoryRouter initialEntries={[initialPath]}>
        <CorpusTopbar />
        <BeamtimeProbe />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("CorpusTopbar", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("shows the corpus wordmark", () => {
    mockExperiments();
    renderTopbar();
    const wordmark = screen.getByTestId("corpus-wordmark");
    expect(wordmark).toHaveTextContent("Himalaya");
    expect(wordmark).toHaveTextContent("SAXS");
  });

  it("renders three stage-tabs; Samples is active and links to /samples", () => {
    mockExperiments();
    renderTopbar();
    expect(screen.getByTestId("stage-tab-samples")).toHaveAttribute(
      "href",
      "/samples",
    );
    expect(screen.getByTestId("stage-tab-index")).toBeDisabled();
    expect(screen.getByTestId("stage-tab-series")).toBeDisabled();
  });

  it("lists experiments in the Beamtime chip once they load", async () => {
    mockExperiments();
    renderTopbar();
    await waitFor(() =>
      expect(screen.getByRole("option", { name: "SSRL Apr 2026" })).toBeInTheDocument(),
    );
    expect(screen.getByRole("option", { name: "APS Jul 2026" })).toBeInTheDocument();
  });

  it("writes ?beamtime= to the URL when an experiment is picked", async () => {
    mockExperiments();
    renderTopbar();
    await screen.findByRole("option", { name: "APS Jul 2026" });
    fireEvent.change(screen.getByTestId("beamtime-chip"), {
      target: { value: "2" },
    });
    expect(screen.getByTestId("beamtime-probe")).toHaveTextContent("2");
  });

  it("clears ?beamtime= when 'all' is picked", async () => {
    mockExperiments();
    renderTopbar("/samples?beamtime=2");
    await screen.findByRole("option", { name: "APS Jul 2026" });
    expect(screen.getByTestId("beamtime-probe")).toHaveTextContent("2");
    fireEvent.change(screen.getByTestId("beamtime-chip"), {
      target: { value: "" },
    });
    expect(screen.getByTestId("beamtime-probe")).toHaveTextContent("");
  });

  it("reflects the current ?beamtime= as the selected option", async () => {
    mockExperiments();
    renderTopbar("/samples?beamtime=1");
    await screen.findByRole("option", { name: "SSRL Apr 2026" });
    expect(screen.getByTestId<HTMLSelectElement>("beamtime-chip").value).toBe("1");
  });
});
```

- [ ] **Step 2: Run the test to verify it fails**

Run: `npm test -- CorpusTopbar`
Expected: FAIL — the chip is still an inert `<button>`; no `<option>` roles, no URL writes.

- [ ] **Step 3: Make the Beamtime chip functional**

In `src/components/CorpusTopbar.tsx`, change the import on line 1 from:

```tsx
import { Link } from "react-router-dom";
```

to:

```tsx
import { Link, useSearchParams } from "react-router-dom";
import { useExperiments } from "../queries";
```

Update the component doc comment (the block ending `…owned by the /samples route (#160).`) to:

```tsx
/**
 * CorpusTopbar — the topbar for the redesigned corpus-scoped shell: a
 * corpus-level wordmark, the three workflow stage-tabs, and the Beamtime
 * facet chip.
 *
 * The Beamtime chip is an experiment picker: it reads and writes the
 * `?beamtime=<experiment_id>` URL query that the /samples contact sheet
 * (#160) filters on. The URL is the only channel — no prop coupling.
 */
```

Inside `CorpusTopbar`, add at the top of the function body (before the `return`):

```tsx
  const [searchParams, setSearchParams] = useSearchParams();
  const experimentsQuery = useExperiments();
  const beamtime = searchParams.get("beamtime") ?? "";

  function handlePick(event: React.ChangeEvent<HTMLSelectElement>): void {
    const value = event.target.value;
    setSearchParams((prev) => {
      const next = new URLSearchParams(prev);
      if (value === "") next.delete("beamtime");
      else next.set("beamtime", value);
      return next;
    });
  }
```

Replace the entire `<button … data-testid="beamtime-chip" …>…</button>` element with:

```tsx
      <select
        data-testid="beamtime-chip"
        aria-label="Filter to a beamtime"
        value={beamtime}
        onChange={handlePick}
        className="rounded-full border border-hair-strong bg-plate px-2.5 py-1
                   text-xs font-semibold text-ink"
      >
        <option value="">Beamtime — all experiments</option>
        {(experimentsQuery.data ?? []).map((exp) => (
          <option key={exp.id} value={exp.id}>
            {exp.name ?? `Experiment ${exp.id}`}
          </option>
        ))}
      </select>
```

- [ ] **Step 4: Run the test to verify it passes**

Run: `npm test -- CorpusTopbar`
Expected: PASS — all 6 `CorpusTopbar` tests green.

- [ ] **Step 5: Commit**

```bash
git add src/components/CorpusTopbar.tsx test/CorpusTopbar.test.tsx
git commit -m "Make the CorpusTopbar Beamtime chip a functional experiment picker (#160)"
```

---

## Task 4: `?beamtime=` round-trip integration test + full verification

**Files:**
- Modify: `packages/HimalayaUI/frontend/test/contact-sheet.test.tsx` (append)

- [ ] **Step 1: Write the failing integration test**

Append to `test/contact-sheet.test.tsx`. Add to the import block at the top:

```tsx
import { Routes, Route } from "react-router-dom";
import { CorpusShell } from "../src/components/CorpusShell";
```

Append this describe block at the end of the file:

```tsx
describe("contact sheet — ?beamtime= round-trip", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("picking a beamtime in the topbar filters the table", async () => {
    mockFetch(corpusRoutes());
    render(
      <QueryClientProvider client={makeClient()}>
        <MemoryRouter initialEntries={["/samples"]}>
          <Routes>
            <Route element={<CorpusShell />}>
              <Route path="/samples" element={<SamplesPage />} />
            </Route>
          </Routes>
        </MemoryRouter>
      </QueryClientProvider>,
    );

    // All three corpus rows render unfiltered.
    await waitFor(() =>
      expect(screen.getByTestId("sample-row-12")).toBeInTheDocument(),
    );
    expect(screen.getByTestId("sample-row-10")).toBeInTheDocument();

    // Pick experiment 2 in the topbar chip.
    fireEvent.change(screen.getByTestId("beamtime-chip"), {
      target: { value: "2" },
    });

    // The table now shows only experiment 2's sample.
    await waitFor(() =>
      expect(screen.queryByTestId("sample-row-10")).toBeNull(),
    );
    expect(screen.queryByTestId("sample-row-11")).toBeNull();
    expect(screen.getByTestId("sample-row-12")).toBeInTheDocument();
  });
});
```

`fireEvent` is already imported in Task 1's import block — confirm `import { render, screen, waitFor, fireEvent } from "@testing-library/react";` at the top of the file; add `fireEvent` if missing.

- [ ] **Step 2: Run the test to verify it fails**

Run: `npm test -- contact-sheet`
Expected: the new round-trip test FAILS only if a wiring bug exists; if Tasks 1–3 are correct it may PASS immediately. Either way, confirm the test executes and the assertions run.

- [ ] **Step 3: No new implementation**

This task adds no production code — Tasks 1–3 already implement the round-trip. If Step 2 failed, debug the wiring (the `?beamtime=` write in `CorpusTopbar` and read in `SamplesPage` must share the one `MemoryRouter`).

- [ ] **Step 4: Run the full affected suites**

Run: `npm test -- contact-sheet CorpusShell CorpusTopbar AppRoutes`
Expected: PASS — all four files green. `AppRoutes.test.tsx` mounts `SamplesPage` at `/samples`; it provides a `QueryClientProvider`, so the now-querying `SamplesPage` does not throw (unmocked fetches resolve to the query error state, and the `samples-page` root test id still renders).

- [ ] **Step 5: Run the whole frontend suite**

Run: `npm test`
Expected: PASS — the full Vitest suite is green. If any unrelated file fails, confirm it fails on `main` too before treating it as a regression.

- [ ] **Step 6: Build**

Run: `npm run build`
Expected: PASS — `tsc --noEmit` reports no type errors and `vite build` completes.

- [ ] **Step 7: Commit**

```bash
git add test/contact-sheet.test.tsx
git commit -m "Add the ?beamtime= contact-sheet round-trip integration test (#160)"
```

---

## Self-Review — spec coverage

Every spec acceptance criterion maps to a task:

| Spec acceptance criterion | Task |
|---|---|
| `/samples` renders the whole corpus as a contact-sheet table | Task 2 |
| Each row shows identity, exposure strip, Kept count, Tags chips, Status | Task 1 |
| Kept count = `kept / total` with non-rejected kept, "N dropped" sub-label | Task 1 (steps 1, 3) |
| `?beamtime=<experiment_id>` filters and round-trips through the URL | Task 2 (read) + Task 3 (write) + Task 4 (round-trip) |
| A boneyard skeleton shows while the corpus query is loading | Task 2 (skeleton test + `<Skeleton>` gate) |
| Thumbnail-select and the Tags `+` button render inert | Task 1 (`ExposureThumb` has no `onClick`; `tag-add` is `disabled`) |
| Status renders "Not indexed" behind a typed seam | Task 1 (`SampleStatus` type + `status-cell`) |
| Vitest covers table, Kept count, exposure strip, `?beamtime=`, skeleton | Tasks 1, 2, 4 |
| `npm run build` passes | Task 4 (step 6) |

Spec scope items deliberately **not** implemented (out of scope per the spec): the screened mark / progress bar (#162), culling interactions (#162), tag mutation (#159), the loupe (#161), the real phase call, the Inspect cutover (#163), row virtualization. No task touches these.
