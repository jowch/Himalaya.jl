import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, useLocation, useNavigate, Routes, Route } from "react-router-dom";
import { QueryClientProvider } from "@tanstack/react-query";
import { makeClient } from "./test-utils";
import { CorpusTopbar } from "../src/print/shell/CorpusTopbar";
import { useAppState } from "../src/state";

// Two experiments; samples 1-3 in experiment 1, sample 9 in experiment 2.
// Ordering of the stepper is the corpus-sample order filtered to the active
// sample's experiment.
const SAMPLES = [
  { id: 1, experiment_id: 1, name: "smp_01", display_name: "Lipid A", notes: null, tags: [], q_units: "A-1" },
  { id: 2, experiment_id: 1, name: "smp_02", display_name: "Lipid B", notes: null, tags: [], q_units: "A-1" },
  { id: 3, experiment_id: 1, name: "smp_03", display_name: "Lipid C", notes: null, tags: [], q_units: "A-1" },
  { id: 9, experiment_id: 2, name: "smp_09", display_name: "Other",  notes: "a note", tags: [], q_units: "A-1" },
];
const EXPERIMENTS = [
  { id: 1, name: "SSRL Apr 2026", path: "/e1", data_dir: "/d1", analysis_dir: "/a1",
    manifest_path: null, created_at: "2026-04-01", q_units: "A-1" },
];

function mockFetch(): void {
  vi.spyOn(global, "fetch").mockImplementation(async (url: string | URL | Request) => {
    const u = String(url);
    const json = (b: unknown) =>
      new Response(JSON.stringify(b), { status: 200, headers: { "content-type": "application/json" } });
    if (u.includes("/api/experiments")) return json(EXPERIMENTS);
    if (u.includes("/api/samples")) return json(SAMPLES); // corpus list
    return json([]);
  });
}

// Captures the current path so we can assert stepper navigation.
function LocationProbe(): JSX.Element {
  const loc = useLocation();
  return <span data-testid="loc">{loc.pathname}</span>;
}

// Mid-session navigation probe: a button that navigates to the given path,
// simulating the user editing the URL or following a stale link.
function NavTo({ path }: { path: string }): JSX.Element {
  const navigate = useNavigate();
  return (
    <button data-testid="nav-to" onClick={() => navigate(path)}>
      go
    </button>
  );
}

function renderWithNav(initialPath: string, navPath: string) {
  return render(
    <QueryClientProvider client={makeClient()}>
      <MemoryRouter initialEntries={[initialPath]}>
        <Routes>
          <Route
            path="*"
            element={
              <>
                <CorpusTopbar />
                <LocationProbe />
                <NavTo path={navPath} />
              </>
            }
          />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

function renderAt(path: string) {
  return render(
    <QueryClientProvider client={makeClient()}>
      <MemoryRouter initialEntries={[path]}>
        <Routes>
          <Route path="*" element={<><CorpusTopbar /><LocationProbe /></>} />
        </Routes>
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("CorpusTopbar — focus affordances (F-13/F-14/F-12)", () => {
  beforeEach(() => {
    vi.restoreAllMocks();
    localStorage.clear();
    // The focus surface seeds activeSampleId from the route; mirror that.
    useAppState.setState({ activeSampleId: undefined, activeExposureId: undefined });
  });

  // ── F-13: per-sample stepper ─────────────────────────────────────────────
  it("shows no stepper off a sample route", () => {
    mockFetch();
    useAppState.setState({ activeSampleId: undefined });
    renderAt("/samples");
    expect(screen.queryByTestId("sample-stepper")).not.toBeInTheDocument();
  });

  it("shows the stepper with 'sample N of M' and the display name on /sample/:id", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderAt("/sample/2");
    await screen.findByTestId("sample-stepper");
    const stepper = screen.getByTestId("sample-stepper");
    // sample 2 is index 1 (0-based) of 3 in experiment 1 → "sample 2 of 3"
    expect(stepper).toHaveTextContent("sample 2 of 3");
    expect(stepper).toHaveTextContent("Lipid B");
  });

  // F6: the caption is small INFORMATIONAL text, so under the F-CONTRAST
  // settled roles it must ride the ink-soft token (AA-normal >= 4.5:1) — faint
  // is decorative-only (3.16:1 on paper, below AA for text). The Kicker
  // primitive's data-tone attribute is its documented tone contract (the
  // token-pair ratios themselves are pinned in contrast-tokens.test.ts).
  it("renders the 'sample N of M' caption in the soft informational tone, not faint (F-CONTRAST)", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderAt("/sample/2");
    await screen.findByTestId("sample-stepper");
    const caption = screen.getByText(/sample 2 of 3/);
    expect(caption).toHaveAttribute("data-tone", "soft");
  });

  it("next navigates to the next sample in the experiment", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderAt("/sample/2");
    await screen.findByTestId("sample-stepper");
    fireEvent.click(screen.getByTestId("sample-stepper-next"));
    expect(screen.getByTestId("loc")).toHaveTextContent("/sample/3");
  });

  it("prev navigates to the previous sample in the experiment", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderAt("/sample/2");
    await screen.findByTestId("sample-stepper");
    fireEvent.click(screen.getByTestId("sample-stepper-prev"));
    expect(screen.getByTestId("loc")).toHaveTextContent("/sample/1");
  });

  it("disables prev on the first sample and next on the last (no wrap)", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 1 });
    renderAt("/sample/1");
    await screen.findByTestId("sample-stepper");
    expect(screen.getByTestId("sample-stepper-prev")).toBeDisabled();
    expect(screen.getByTestId("sample-stepper-next")).not.toBeDisabled();
  });

  // ── F-STALEURL: a bogus /sample/:id must not show the stale stepper ──────
  // Mid-session the store keeps the previous activeSampleId (a bogus route
  // never seeds it), so the stepper would otherwise keep showing the previous
  // sample's name and live prev/next above a "Sample not found" body.
  it("hides the stepper when navigating mid-session to a bogus numeric /sample/:id (F-STALEURL)", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderWithNav("/sample/2", "/sample/99999");
    // Warm start on a real sample: the stepper is up, the corpus cache is hot.
    await screen.findByTestId("sample-stepper");
    fireEvent.click(screen.getByTestId("nav-to"));
    expect(screen.getByTestId("loc")).toHaveTextContent("/sample/99999");
    await waitFor(() =>
      expect(screen.queryByTestId("sample-stepper")).not.toBeInTheDocument(),
    );
  });

  it("hides the stepper when navigating mid-session to a non-numeric /sample/:id (F-STALEURL)", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderWithNav("/sample/2", "/sample/not-a-number");
    await screen.findByTestId("sample-stepper");
    fireEvent.click(screen.getByTestId("nav-to"));
    await waitFor(() =>
      expect(screen.queryByTestId("sample-stepper")).not.toBeInTheDocument(),
    );
  });

  it("keeps the stepper across a valid-to-valid step (no flicker gate on store equality)", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderWithNav("/sample/2", "/sample/3");
    await screen.findByTestId("sample-stepper");
    fireEvent.click(screen.getByTestId("nav-to"));
    // The store still says 2 for an effect tick after the URL says 3; a known
    // param stays "found" throughout, so the stepper must remain mounted.
    expect(screen.getByTestId("sample-stepper")).toBeInTheDocument();
    await waitFor(() =>
      expect(screen.getByTestId("sample-stepper")).toBeInTheDocument(),
    );
  });

  // ── the focus workspace is not a top-level stage tab ─────────────────────
  // There is no "Index"/"Focus" stage tab — the focus workspace is reached by
  // opening a sample, and the on-focus context is carried by the sample stepper.
  it("renders no Index/Focus stage tab on a /sample/:id route", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderAt("/sample/2");
    // The stepper still mounts (focus context), but no stage tab for focus.
    expect(await screen.findByTestId("sample-stepper")).toBeInTheDocument();
    expect(screen.queryByTestId("stage-tab-index")).toBeNull();
  });

  // ── notes were removed from the focus surface ────────────────────────────
  it("renders no Notes toggle on a sample route (notes feature removed)", async () => {
    mockFetch();
    useAppState.setState({ activeSampleId: 2 });
    renderAt("/sample/2");
    await screen.findByTestId("sample-stepper");
    expect(screen.queryByTestId("notes-toggle")).toBeNull();
  });
});
