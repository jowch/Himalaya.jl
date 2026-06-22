// test/ScanFailedPage.test.tsx
import { describe, it, expect, vi, beforeEach, afterEach } from "vitest";
import { render, screen, fireEvent, waitFor } from "@testing-library/react";
import { MemoryRouter, Routes, Route } from "react-router-dom";
import { ScanFailedPage } from "../src/print/pages/ScanFailedPage";
import * as api from "../src/api";
import type { ManifestUnmatched } from "../src/api";

// Module-level navigate mock — canonical pattern.
const navigate = vi.fn();
vi.mock("react-router-dom", async (orig) => {
  const m = await orig<typeof import("react-router-dom")>();
  return { ...m, useNavigate: () => navigate };
});

function renderScanFailed({
  unmatched,
  parsedCount,
  dataDir,
  analysisDir,
  patterns,
}: {
  unmatched: ManifestUnmatched[];
  parsedCount: number;
  dataDir?: string;
  analysisDir?: string;
  patterns?: { image?: string; metadata?: string; integration?: string };
}) {
  render(
    <MemoryRouter initialEntries={["/experiments/7/corpus"]}>
      <Routes>
        <Route
          path="/experiments/:id/corpus"
          element={
            <ScanFailedPage
              experimentId={7}
              unmatched={unmatched}
              parsedCount={parsedCount}
              {...(dataDir != null ? { dataDir } : {})}
              {...(analysisDir != null ? { analysisDir } : {})}
              {...(patterns != null ? { patterns } : {})}
            />
          }
        />
      </Routes>
    </MemoryRouter>,
  );
}

describe("ScanFailedPage", () => {
  beforeEach(() => {
    navigate.mockClear();
  });
  afterEach(() => vi.restoreAllMocks());

  it("groups misses, offers per-type test + ingest-N confirm", () => {
    renderScanFailed({
      unmatched: [
        { file: "a", miss: "metadata" },
        { file: "b", miss: "integration" },
      ],
      parsedCount: 5,
    });

    // Open configuration primary action (the bottom action-bar CTA).
    expect(
      screen.getByRole("button", { name: /open configuration/i }),
    ).toBeInTheDocument();

    // Unmatched stem shown.
    expect(screen.getByText("a")).toBeInTheDocument();

    // Miss-type section heading shown (scoped via stable testid — the summary
    // sentence also mentions "metadata", so getByText(/metadata/i) is ambiguous).
    expect(screen.getByTestId("miss-heading-metadata")).toBeInTheDocument();

    // One textbox per affected type (2 types: metadata + integration). The live
    // ✓ readout must NOT add a textbox.
    expect(screen.getAllByRole("textbox").length).toBe(2);

    // Two-stage confirm: first click shows Confirm.
    fireEvent.click(screen.getByRole("button", { name: /ingest 5/i }));
    expect(screen.getByText(/confirm/i)).toBeInTheDocument();
  });

  it("renders only one textbox when only one miss type is present", () => {
    renderScanFailed({
      unmatched: [
        { file: "x", miss: "metadata" },
        { file: "y", miss: "metadata" },
      ],
      parsedCount: 3,
    });
    expect(screen.getAllByRole("textbox").length).toBe(1);
  });

  it("navigates to config tab when Open configuration is clicked", () => {
    renderScanFailed({ unmatched: [], parsedCount: 0 });
    fireEvent.click(screen.getByRole("button", { name: /open configuration/i }));
    expect(navigate).toHaveBeenCalledWith("/experiments/7/config");
  });

  it("calls onIngestParsed after two-stage confirm when prop is provided", () => {
    const onIngestParsed = vi.fn();
    render(
      <MemoryRouter initialEntries={["/experiments/7/corpus"]}>
        <Routes>
          <Route
            path="/experiments/:id/corpus"
            element={
              <ScanFailedPage
                experimentId={7}
                unmatched={[{ file: "a.prp", miss: "metadata" }]}
                parsedCount={2}
                onIngestParsed={onIngestParsed}
              />
            }
          />
        </Routes>
      </MemoryRouter>,
    );
    // Arm the confirm.
    fireEvent.click(screen.getByRole("button", { name: /ingest 2/i }));
    // Confirm executes the callback.
    fireEvent.click(screen.getByTestId("ingest-confirm-yes"));
    expect(onIngestParsed).toHaveBeenCalledTimes(1);
  });

  it("shows a live ✓ N / N readout when a trial pattern clears (dataDir present)", async () => {
    vi.spyOn(api, "fetchManifest").mockResolvedValue({
      total: 19,
      matched: { image: 19, metadata: 0, integration: 19 },
      unmatched: [],
    });
    renderScanFailed({
      unmatched: [{ file: "HA_5_044_S1991", miss: "integration" }],
      parsedCount: 19,
      dataDir: "/data/exp1/raw",
      analysisDir: "/data/exp1/analysis",
      patterns: { image: "*.tif", integration: "*.dat" },
    });

    const input = screen.getByRole("textbox", { name: /integration pattern/i });
    fireEvent.change(input, { target: { value: "*_total.dat" } });

    // After the debounced fetch resolves, a cleared ✓ readout appears.
    await waitFor(() =>
      expect(screen.getByText(/✓ 19 \/ 19/)).toBeInTheDocument(),
    );
    expect(api.fetchManifest).toHaveBeenCalledWith(
      "/data/exp1/raw",
      expect.objectContaining({ integration: "*_total.dat" }),
      undefined,
      "/data/exp1/analysis",
    );
  });
});
