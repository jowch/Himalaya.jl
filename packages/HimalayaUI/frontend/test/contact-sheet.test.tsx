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
});
