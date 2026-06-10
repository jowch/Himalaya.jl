// test/SamplesPage.samplePicker.test.tsx
import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, within } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { QueryClient, QueryClientProvider } from "@tanstack/react-query";
import { MemoryRouter } from "react-router-dom";
import { SamplesPage } from "../src/print/pages/SamplesPage";
import * as queries from "../src/queries";

// Minimal corpus: two samples, no exposures.
const SAMPLE_A = {
  id: 10, experiment_id: 1, name: "A", display_name: null,
  notes: null, tags: [], q_units: "A-1",
};
const SAMPLE_B = {
  id: 11, experiment_id: 1, name: "B", display_name: null,
  notes: null, tags: [], q_units: "A-1",
};

function setup() {
  const qc = new QueryClient({ defaultOptions: { queries: { retry: false } } });
  vi.spyOn(queries, "useCorpusSamples").mockReturnValue({
    data: [SAMPLE_A, SAMPLE_B], isLoading: false, isError: false,
  } as ReturnType<typeof queries.useCorpusSamples>);
  vi.spyOn(queries, "useCorpusExposures").mockReturnValue(
    { byId: new Map(), isLoading: false });
  vi.spyOn(queries, "useScreenedProgress").mockReturnValue({ screened: 0, total: 2 });
  vi.spyOn(queries, "useExperiments").mockReturnValue(
    { data: [], isLoading: false } as unknown as ReturnType<typeof queries.useExperiments>);
  vi.spyOn(queries, "useSetExposureStatusBatch").mockReturnValue(
    { mutate: vi.fn() } as unknown as ReturnType<typeof queries.useSetExposureStatusBatch>);

  return render(
    <QueryClientProvider client={qc}>
      <MemoryRouter>
        <SamplesPage />
      </MemoryRouter>
    </QueryClientProvider>,
  );
}

describe("SamplesPage — sample-grain picker", () => {
  beforeEach(() => vi.restoreAllMocks());

  it("compose bar is hidden when no samples are checked", () => {
    setup();
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
  });

  it("compose bar becomes visible when a sample is checked", async () => {
    const user = userEvent.setup();
    setup();
    const checkboxes = screen.getAllByRole("checkbox");
    await user.click(checkboxes[0]!);
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "true");
    expect(screen.getByTestId("compose-bar")).toHaveTextContent("1");
  });

  it("Clear resets selection and hides the compose bar", async () => {
    const user = userEvent.setup();
    setup();
    await user.click(screen.getAllByRole("checkbox")[0]!);
    // Scope to the compose bar: the hidden CullBar also has a Clear button,
    // and JSDOM does not implement inert a11y-tree exclusion (real browsers
    // do — pinned by the corpus-culling e2e).
    await user.click(
      within(screen.getByTestId("compose-bar")).getByRole("button", { name: /clear/i }),
    );
    expect(screen.getByTestId("compose-bar")).toHaveAttribute("data-show", "false");
  });

  it("cull bar remains independent of sample-grain selection", async () => {
    const user = userEvent.setup();
    setup();
    await user.click(screen.getAllByRole("checkbox")[0]!);
    // Frame-level CullBar should still be hidden (no exposure selected).
    expect(screen.getByTestId("cull-bar")).toHaveAttribute("data-show", "false");
  });
});
