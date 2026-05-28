import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, waitFor, fireEvent } from "@testing-library/react";
import { QueryClientProvider } from "@tanstack/react-query";
import type { ReactNode } from "react";
import { MemoryRouter } from "react-router-dom";
import { Routes, Route } from "react-router-dom";
import { makeClient } from "./test-utils";
import { ContactSheetRow } from "../src/components/ContactSheetRow";
import type { CorpusSample, Exposure } from "../src/api";
import { SamplesPage } from "../src/pages/SamplesPage";
import { CorpusShell } from "../src/components/CorpusShell";

/** Render a row inside a router so loupe navigation can be asserted on. */
function renderRowRouted(sample: CorpusSample, initialPath = "/samples") {
  const client = makeClient();
  return render(
    <QueryClientProvider client={client}>
      <MemoryRouter initialEntries={[initialPath]}>
        <Routes>
          <Route
            path="/samples"
            element={<ContactSheetRow sample={sample} />}
          />
          <Route
            path="/samples/loupe/:sampleId"
            element={<div data-testid="loupe-stub" />}
          />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

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
    <QueryClientProvider client={client}>
      <MemoryRouter>{children}</MemoryRouter>
    </QueryClientProvider>
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

  it("renders a 'Not indexed' status with a leading hollow dot", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(makeSample({ id: 7 }));
    const cell = await screen.findByTestId("status-cell");
    expect(cell).toHaveTextContent("Not indexed");
    // M-6: the unset status carries a hollow dot, not just bare text.
    expect(cell.querySelector('[data-testid="status-dot"]')).toBeTruthy();
  });

  it("renders a phase chip when the sample carries a phase (M-6)", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(makeSample({ id: 7, phase: "Pn3m" }));
    const cell = await screen.findByTestId("status-cell");
    const chip = cell.querySelector('[data-testid="phase-chip"]');
    expect(chip).toBeTruthy();
    expect(chip).toHaveTextContent("Pn3m");
    expect(cell).not.toHaveTextContent("Not indexed");
  });

  // M-2: per-sample screened mark + unscreened-row tint.
  it("renders a filled screened mark for a screened sample", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(makeSample({ id: 7, screened: true }));
    const mark = await screen.findByTestId("screened-mark");
    expect(mark).toHaveAttribute("data-screened", "true");
  });

  it("renders a hollow screened mark and tints an unscreened row", async () => {
    mockFetch({ "/api/samples/7/exposures": [] });
    renderRow(makeSample({ id: 7, screened: false }));
    const mark = await screen.findByTestId("screened-mark");
    expect(mark).not.toHaveAttribute("data-screened", "true");
    expect(screen.getByTestId("sample-row-7")).toHaveAttribute(
      "data-unscreened",
      "true",
    );
  });

  it("derives screened from all exposures being triaged when no flag is set", async () => {
    // No explicit `screened` flag: a sample whose every exposure has a
    // non-null status reads as screened (M-2 derivation pending #162 backend).
    mockFetch({
      "/api/samples/7/exposures": [
        makeExposure({ id: 1, sample_id: 7, status: "accepted" }),
        makeExposure({ id: 2, sample_id: 7, status: "rejected" }),
      ],
    });
    renderRow(makeSample({ id: 7 }));
    await waitFor(() =>
      expect(screen.getByTestId("screened-mark")).toHaveAttribute(
        "data-screened",
        "true",
      ),
    );
  });

  // M-7: zero-padded frame-number badge on each thumbnail.
  it("renders a zero-padded frame number badge on thumbnails", async () => {
    mockFetch({
      "/api/samples/7/exposures": [
        makeExposure({ id: 1, sample_id: 7 }),
        makeExposure({ id: 2, sample_id: 7 }),
      ],
    });
    renderRow(makeSample({ id: 7 }));
    const badge = await screen.findByTestId("frame-no-1");
    expect(badge).toHaveTextContent("01");
    expect(await screen.findByTestId("frame-no-2")).toHaveTextContent("02");
  });

  // M-10: grease-pencil X mark renders on a rejected frame.
  it("renders the grease-pencil X mark on a rejected thumbnail", async () => {
    mockFetch({
      "/api/samples/7/exposures": [
        makeExposure({ id: 1, sample_id: 7, status: "rejected" }),
      ],
    });
    renderRow(makeSample({ id: 7 }));
    const thumb = await screen.findByTestId("exposure-thumb-1");
    expect(
      thumb.querySelector('[data-testid="reject-xmark"]'),
    ).toBeTruthy();
  });
});

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
    await waitFor(() =>
      expect(screen.getByTestId("samples-scope")).toHaveTextContent(
        "APS Jul 2026",
      ),
    );
  });

  it("shows an error state when the corpus query fails", async () => {
    mockFetch({}); // every path → 404
    renderSamplesPage();
    expect(await screen.findByTestId("samples-error")).toBeInTheDocument();
  });

  it("shows an empty state when the beamtime filter matches no samples", async () => {
    mockFetch(corpusRoutes());
    renderSamplesPage("/samples?beamtime=999");
    expect(await screen.findByTestId("samples-empty")).toBeInTheDocument();
    expect(screen.queryByTestId("contact-sheet-rows")).toBeNull();
  });

  // L-3: beamtime serif h1 + descriptive sub.
  it("renders a beamtime title and a descriptive sub", async () => {
    mockFetch(corpusRoutes());
    renderSamplesPage("/samples?beamtime=2");
    await waitFor(() =>
      expect(screen.getByTestId("samples-title")).toHaveTextContent(
        "APS Jul 2026",
      ),
    );
    expect(screen.getByTestId("samples-sub")).toBeInTheDocument();
  });

  // M-1: "N / M screened" progress block + terracotta bar.
  it("renders a screened-progress count and bar", async () => {
    mockFetch(corpusRoutes());
    renderSamplesPage();
    const progress = await screen.findByTestId("screened-progress");
    // 3 corpus samples, each with exposures that have no status → 0 screened.
    await waitFor(() =>
      expect(progress).toHaveTextContent("/ 3"),
    );
    expect(screen.getByTestId("screened-progress-bar")).toBeInTheDocument();
  });

  // L-5/L-8/L-11: keycap-chip legend with the 5 hints (incl. shift-click range).
  it("renders the keyboard legend with a shift-click range hint", async () => {
    mockFetch(corpusRoutes());
    renderSamplesPage();
    const legend = await screen.findByTestId("kb-legend");
    expect(legend).toHaveTextContent("extend the range");
    expect(legend.querySelectorAll("[data-kb-key]").length).toBe(5);
  });
});

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

// I1.6 (#162) — culling wiring. These tests prove the optimistic cache patch
// plus the outbound HTTP request through the real useQueueMutation → mutator →
// api → fetch path. They do NOT exercise the full queue lifecycle (no SSE
// confirmation / clearDeferred / replay) — that is covered by the queue suite.
describe("ContactSheetRow — culling", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("rejects a single exposure (optimistic patch + PATCH status request)", async () => {
    const patched: Array<{ url: string; body: unknown }> = [];
    vi.spyOn(global, "fetch").mockImplementation(
      (input: RequestInfo | URL, init?: RequestInit) => {
        const url = typeof input === "string" ? input : (input as Request).url;
        if (url === "/api/samples/7/exposures") {
          return Promise.resolve(
            new Response(
              JSON.stringify([
                makeExposure({ id: 1, sample_id: 7, status: "accepted" }),
              ]),
              { status: 200, headers: { "Content-Type": "application/json" } },
            ),
          );
        }
        if (/\/api\/exposures\/1\/status$/.test(url)) {
          patched.push({
            url,
            body: init?.body ? JSON.parse(String(init.body)) : null,
          });
          return Promise.resolve(
            new Response(JSON.stringify({ id: 1, status: "rejected" }), {
              status: 200,
              headers: { "Content-Type": "application/json" },
            }),
          );
        }
        return Promise.resolve(new Response("not found", { status: 404 }));
      },
    );

    renderRow(makeSample({ id: 7 }));
    const reject = await screen.findByTestId("exposure-reject-1");
    fireEvent.click(reject);

    await waitFor(() =>
      expect(
        patched.some((p) => /\/api\/exposures\/1\/status$/.test(p.url)),
      ).toBe(true),
    );
    expect(patched[0].body).toMatchObject({ status: "rejected" });
    await waitFor(() =>
      expect(screen.getByTestId("exposure-thumb-1")).toHaveAttribute(
        "data-rejected",
        "true",
      ),
    );
  });

  it("batch-rejects multiple selected exposures (N independent status requests)", async () => {
    const patched: string[] = [];
    vi.spyOn(global, "fetch").mockImplementation((input: RequestInfo | URL) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url === "/api/samples/7/exposures") {
        return Promise.resolve(
          new Response(
            JSON.stringify([
              makeExposure({ id: 1, sample_id: 7, status: "accepted" }),
              makeExposure({ id: 2, sample_id: 7, status: "accepted" }),
              makeExposure({ id: 3, sample_id: 7, status: "accepted" }),
            ]),
            { status: 200, headers: { "Content-Type": "application/json" } },
          ),
        );
      }
      const m = url.match(/\/api\/exposures\/(\d+)\/status$/);
      if (m) {
        patched.push(m[1]);
        return Promise.resolve(
          new Response(JSON.stringify({ id: Number(m[1]), status: "rejected" }), {
            status: 200,
            headers: { "Content-Type": "application/json" },
          }),
        );
      }
      return Promise.resolve(new Response("not found", { status: 404 }));
    });

    renderRow(makeSample({ id: 7 }));
    // Select exposures 1 and 3 for batch action.
    fireEvent.click(await screen.findByTestId("exposure-select-1"));
    fireEvent.click(screen.getByTestId("exposure-select-3"));
    // Action bar appears; reject the selection.
    fireEvent.click(screen.getByTestId("batch-reject"));

    await waitFor(() => expect(patched.slice().sort()).toEqual(["1", "3"]));
    // Exposure 2 was NOT selected → no op.
    expect(patched).not.toContain("2");
    // Selection clears after the batch.
    await waitFor(() =>
      expect(screen.queryByTestId("batch-reject")).toBeNull(),
    );
  });

  it("picks a representative exposure (optimistic mutual exclusion + PATCH select)", async () => {
    let repUrl = "";
    vi.spyOn(global, "fetch").mockImplementation((input: RequestInfo | URL) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url === "/api/samples/7/exposures") {
        return Promise.resolve(
          new Response(
            JSON.stringify([
              makeExposure({ id: 1, sample_id: 7, selected: true }),
              makeExposure({ id: 2, sample_id: 7, selected: false }),
            ]),
            { status: 200, headers: { "Content-Type": "application/json" } },
          ),
        );
      }
      if (/\/api\/exposures\/2\/select$/.test(url)) {
        repUrl = url;
        return Promise.resolve(
          new Response(JSON.stringify({ id: 2, selected: true }), {
            status: 200,
            headers: { "Content-Type": "application/json" },
          }),
        );
      }
      return Promise.resolve(new Response("not found", { status: 404 }));
    });

    renderRow(makeSample({ id: 7 }));
    fireEvent.click(await screen.findByTestId("exposure-represent-2"));

    await waitFor(() =>
      expect(repUrl).toMatch(/\/api\/exposures\/2\/select$/),
    );
    // Optimistic mutual exclusion: 2 becomes representative, 1 stops being it.
    await waitFor(() =>
      expect(screen.getByTestId("exposure-thumb-2")).toHaveAttribute(
        "data-representative",
        "true",
      ),
    );
    expect(screen.getByTestId("exposure-thumb-1")).not.toHaveAttribute(
      "data-representative",
    );
  });
});

// R1 round 2 — the keyboard/pointer affordances the footer legend + CullBar
// advertise must actually be wired (no silent no-ops). Selection stays
// per-sample (a keydown/range-select never crosses rows).
describe("ContactSheetRow — advertised affordances", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
  });

  it("selects a frame by clicking the thumbnail body (not just the checkbox)", async () => {
    mockFetch({
      "/api/samples/7/exposures": [makeExposure({ id: 1, sample_id: 7 })],
    });
    renderRowRouted(makeSample({ id: 7 }));
    const thumb = await screen.findByTestId("exposure-thumb-1");
    fireEvent.click(thumb);
    expect(thumb).toHaveAttribute("data-batch-selected", "true");
    expect(screen.getByTestId("cull-bar")).toBeInTheDocument();
  });

  it("clears the selection on Esc", async () => {
    mockFetch({
      "/api/samples/7/exposures": [makeExposure({ id: 1, sample_id: 7 })],
    });
    renderRowRouted(makeSample({ id: 7 }));
    fireEvent.click(await screen.findByTestId("exposure-select-1"));
    expect(screen.getByTestId("cull-bar")).toBeInTheDocument();
    fireEvent.keyDown(window, { key: "Escape" });
    await waitFor(() => expect(screen.queryByTestId("cull-bar")).toBeNull());
  });

  it("batch-rejects the selection on X", async () => {
    const patched: string[] = [];
    vi.spyOn(global, "fetch").mockImplementation((input: RequestInfo | URL) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url === "/api/samples/7/exposures") {
        return Promise.resolve(
          new Response(
            JSON.stringify([
              makeExposure({ id: 1, sample_id: 7, status: "accepted" }),
              makeExposure({ id: 2, sample_id: 7, status: "accepted" }),
            ]),
            { status: 200, headers: { "Content-Type": "application/json" } },
          ),
        );
      }
      const m = url.match(/\/api\/exposures\/(\d+)\/status$/);
      if (m) {
        patched.push(m[1]);
        return Promise.resolve(
          new Response(JSON.stringify({ id: Number(m[1]), status: "rejected" }), {
            status: 200,
            headers: { "Content-Type": "application/json" },
          }),
        );
      }
      return Promise.resolve(new Response("not found", { status: 404 }));
    });

    renderRowRouted(makeSample({ id: 7 }));
    fireEvent.click(await screen.findByTestId("exposure-select-1"));
    fireEvent.keyDown(window, { key: "x" });

    await waitFor(() => expect(patched).toEqual(["1"]));
    // Selection clears after the batch (cull bar gone).
    await waitFor(() => expect(screen.queryByTestId("cull-bar")).toBeNull());
  });

  it("does not bind X / Esc when there is no selection", async () => {
    const patched: string[] = [];
    vi.spyOn(global, "fetch").mockImplementation((input: RequestInfo | URL) => {
      const url = typeof input === "string" ? input : (input as Request).url;
      if (url === "/api/samples/7/exposures") {
        return Promise.resolve(
          new Response(
            JSON.stringify([
              makeExposure({ id: 1, sample_id: 7, status: "accepted" }),
            ]),
            { status: 200, headers: { "Content-Type": "application/json" } },
          ),
        );
      }
      const m = url.match(/\/api\/exposures\/(\d+)\/status$/);
      if (m) {
        patched.push(m[1]);
        return Promise.resolve(
          new Response(JSON.stringify({ id: Number(m[1]), status: "rejected" }), {
            status: 200,
            headers: { "Content-Type": "application/json" },
          }),
        );
      }
      return Promise.resolve(new Response("not found", { status: 404 }));
    });

    renderRowRouted(makeSample({ id: 7 }));
    await screen.findByTestId("exposure-thumb-1");
    // No selection: X is inert (no batch reject fires).
    fireEvent.keyDown(window, { key: "x" });
    await new Promise((r) => setTimeout(r, 0));
    expect(patched).toEqual([]);
  });

  it("extends a contiguous range on shift-click", async () => {
    mockFetch({
      "/api/samples/7/exposures": [
        makeExposure({ id: 1, sample_id: 7 }),
        makeExposure({ id: 2, sample_id: 7 }),
        makeExposure({ id: 3, sample_id: 7 }),
        makeExposure({ id: 4, sample_id: 7 }),
      ],
    });
    renderRowRouted(makeSample({ id: 7 }));
    // Anchor on frame 1, then shift-click frame 3 → 1,2,3 selected; 4 not.
    fireEvent.click(await screen.findByTestId("exposure-thumb-1"));
    fireEvent.click(screen.getByTestId("exposure-thumb-3"), { shiftKey: true });
    await waitFor(() => {
      expect(screen.getByTestId("exposure-thumb-1")).toHaveAttribute(
        "data-batch-selected",
        "true",
      );
      expect(screen.getByTestId("exposure-thumb-2")).toHaveAttribute(
        "data-batch-selected",
        "true",
      );
      expect(screen.getByTestId("exposure-thumb-3")).toHaveAttribute(
        "data-batch-selected",
        "true",
      );
    });
    expect(screen.getByTestId("exposure-thumb-4")).not.toHaveAttribute(
      "data-batch-selected",
    );
  });

  it("opens the loupe on double-click of a thumbnail", async () => {
    mockFetch({
      "/api/samples/7/exposures": [makeExposure({ id: 1, sample_id: 7 })],
    });
    renderRowRouted(makeSample({ id: 7 }));
    const thumb = await screen.findByTestId("exposure-thumb-1");
    fireEvent.doubleClick(thumb);
    expect(await screen.findByTestId("loupe-stub")).toBeInTheDocument();
  });
});
