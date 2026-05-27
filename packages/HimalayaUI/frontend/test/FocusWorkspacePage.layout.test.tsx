import { describe, it, expect, beforeEach, vi } from "vitest";
import { render, screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { MemoryRouter } from "react-router-dom";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { AppRoutes } from "../src/components/AppRoutes";
import { useAppState } from "../src/state";

// Corpus list (keyed by sampleId alone) — how the layout learns experiment_id.
const CORPUS = [
  { id: 1, experiment_id: 1, name: "JC001", display_name: "Sample 1",
    notes: "a note", tags: [], q_units: "A-1" },
];
// Experiment-scoped list — the cache `updateSampleMutator` patches and the
// one the notes textarea reads from. Same row, different query key.
const EXP_SAMPLES = [
  { id: 1, experiment_id: 1, name: "JC001", display_name: "Sample 1",
    notes: "a note", tags: [] },
];
const EXPOSURES = [
  { id: 5, sample_id: 1, filename: "JC001-005.dat", kind: "file",
    selected: true, status: "accepted", image_path: "/img/5.png",
    image_version: "v1", tags: [], sources: [],
    trace_hash: null, analysis_inputs_hash: null },
];
// Trace shape is { q, I, sigma } (api.ts:175) — NOT `intensity`.
const TRACE = { q: [0.04, 0.1, 0.3], I: [10, 80, 20], sigma: [1, 1, 1] };

beforeEach(() => {
  localStorage.clear();
  // Seed NOTHING but username/theme: the route shim seeds activeSampleId,
  // and activeExperimentId is intentionally left undefined to mirror the real
  // /sample/:sampleId route (the I4.1 shim does not seed it).
  useAppState.setState({
    activeSampleId: undefined, activeExposureId: undefined,
    activeExperimentId: undefined, username: "tester", theme: "dark",
  });
  vi.stubGlobal("ResizeObserver", class {
    observe() {} unobserve() {} disconnect() {}
  });
  vi.stubGlobal("fetch", vi.fn(async (url: string, init?: RequestInit) => {
    const u = String(url);
    const json = (b: unknown) => new Response(JSON.stringify(b), {
      status: 200, headers: { "content-type": "application/json" } });
    // Order matters: the experiment-scoped route (/api/experiments/:id/samples)
    // must be matched BEFORE the corpus route (/api/samples), and PATCH writes
    // echo the patched row back.
    if (/\/api\/experiments\/\d+\/samples/.test(u)) return json(EXP_SAMPLES);
    if (/\/api\/samples\/\d+$/.test(u) && init?.method === "PATCH") {
      const patch = JSON.parse(String(init.body));
      return json({ ...EXP_SAMPLES[0], ...patch });
    }
    if (u.includes("/api/samples")) return json(CORPUS);
    // The sub-resources of /api/exposures/:id/* must be matched BEFORE the
    // bare exposures-list route (which the loose /exposures check below
    // would otherwise swallow, handing PlotCard the wrong shape).
    if (u.includes("/trace")) return json(TRACE);
    if (u.includes("/peaks")) return json([]);   // no peaks: TraceViewer renders empty
    if (u.includes("/indices")) return json([]); // no indices: PhasePanel empty state
    if (u.includes("/groups")) return json([]);  // no active group
    if (/\/exposures(\?|$)/.test(u)) return json(EXPOSURES); // GET /api/samples/:id/exposures
    return json([]);
  }));
});

function renderAt(path: string) {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={[path]}>
        <AppRoutes />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("focus workspace layout", () => {
  it("composes the trace hero, detector, rail, and notes for a sample", async () => {
    renderAt("/sample/1");
    // findBy* (async, retry-aware) — accounts for the shim -> auto-pick ->
    // corpus -> experiment-scoped settling chain before asserting presence.
    expect(await screen.findByTestId("focus-workspace-page")).toBeInTheDocument();
    expect(await screen.findByTestId("plot-card")).toBeInTheDocument();
    expect(await screen.findByTestId("indices-card")).toBeInTheDocument();
    expect(await screen.findByTestId("focus-detector-panel")).toBeInTheDocument();
    // Notes margin renders only once the experiment-scoped sample resolves —
    // its presence proves the corpus->experiment_id->useSamples chain works.
    expect(await screen.findByTestId("focus-notes-margin")).toBeInTheDocument();
  });

  // BLOCKING-review regression: the notes textarea must read from the SAME
  // cache the save patches. If it read from useCorpusSamples (which the
  // mutator never patches), the on-blur focus-gate would re-sync the stale
  // corpus value and silently revert the edit. This asserts the round-trip:
  // edit -> blur(save) -> the new value STICKS.
  it("retains an edited note after save (read-source == write-target)", async () => {
    const user = userEvent.setup();
    renderAt("/sample/1");
    const ta = await screen.findByTestId("focus-notes-input") as HTMLTextAreaElement;
    await user.click(ta);
    await user.clear(ta);
    await user.type(ta, "edited in the focus workspace");
    await user.tab(); // blur -> updateSample.mutate({ notes })
    // The optimistic mutator patch updates queryKeys.samples(1); useSamples(1)
    // re-delivers the patched row; the focus-gate (now unfocused) syncs it —
    // and because read-source == write-target, the synced value is the EDIT,
    // not a stale corpus value.
    await waitFor(() =>
      expect(ta.value).toBe("edited in the focus workspace"),
    );
  });
});
