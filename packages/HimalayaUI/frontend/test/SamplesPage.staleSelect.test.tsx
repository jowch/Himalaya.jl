// test/SamplesPage.staleSelect.test.tsx
//
// SA-STALESELECT: a beamtime-filter change must NOT carry a working selection
// from the prior scope. Selection that survives the filter is invisible (no row
// to uncheck) and makes the floating bars + the cull toast lie about what's
// actionable. The page clears both selection grains whenever `beamtime` changes.
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter, useSearchParams } from "react-router-dom";
import { SamplesPage } from "../src/print/pages/SamplesPage";
import * as queries from "../src/queries";

// Two samples in DISTINCT experiments, so a ?beamtime=2 filter drops SAMPLE_A
// out of the visible scope (and the corpus vouches for id 2 → a real filter).
const SAMPLE_A = {
  id: 10, experiment_id: 1, name: "A", display_name: null,
  notes: null, tags: [], q_units: "A-1",
};
const SAMPLE_B = {
  id: 11, experiment_id: 2, name: "B", display_name: null,
  notes: null, tags: [], q_units: "A-1",
};

// A sibling control that re-scopes the page the same way the CorpusTopbar
// select does: by writing ?beamtime to the shared URL. (The select itself lives
// in the shell, outside SamplesPage's tree.)
function BeamtimeSwitch(): JSX.Element {
  const [, setSp] = useSearchParams();
  return (
    <button data-testid="set-beamtime" onClick={() => setSp({ experiment: "2" })}>
      filter to exp 2
    </button>
  );
}

function setup() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  vi.spyOn(queries, "useCorpusSamples").mockReturnValue({
    data: [SAMPLE_A, SAMPLE_B], isLoading: false, isError: false,
  } as unknown as ReturnType<typeof queries.useCorpusSamples>);
  vi.spyOn(queries, "useCorpusExposures").mockReturnValue(
    { byId: new Map(), isLoading: false });
  vi.spyOn(queries, "useScreenedProgress").mockReturnValue({ screened: 0, total: 2 });
  vi.spyOn(queries, "useExperiments").mockReturnValue(
    { data: [], isLoading: false } as unknown as ReturnType<typeof queries.useExperiments>);
  vi.spyOn(queries, "useSetExposureStatusBatch").mockReturnValue(
    { mutate: vi.fn() } as unknown as ReturnType<typeof queries.useSetExposureStatusBatch>);

  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter initialEntries={["/samples"]}>
        <BeamtimeSwitch />
        <SamplesPage />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("SamplesPage — SA-STALESELECT (selection clears on filter change)", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("clears the sample-grain compose selection when the beamtime filter changes", async () => {
    const user = userEvent.setup();
    setup();

    // Check SAMPLE_A → the compose bar carries one pick.
    const checkboxes = screen.getAllByRole("checkbox");
    await user.click(checkboxes[0]!);
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");
    expect(screen.getByTestId("compose-bar")).toHaveTextContent("1");

    // Re-scope to experiment 2 (SAMPLE_A leaves the visible list). The carry
    // must NOT survive: no row remains to uncheck it, so it would be a stale
    // "+ New series" carry pointing at an off-filter sample.
    await user.click(screen.getByTestId("set-beamtime"));
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
  });
});
