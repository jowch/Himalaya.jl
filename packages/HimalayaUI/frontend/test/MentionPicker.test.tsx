import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { MentionPicker } from "../src/components/MentionPicker";
import { useAppState } from "../src/state";
import * as api from "../src/api";

const INDICES: api.IndexEntry[] = [
  { id: 17, exposure_id: 1, phase: "Pn3m", basis: 1.0, score: 0.91, r_squared: 0.998,
    lattice_d: 5.14, ngc: -1.51, status: "candidate", peaks: [], predicted_q: [1.223, 1.414] },
];
const PEAKS: api.Peak[] = [
  { id: 42, exposure_id: 1, q: 1.223, intensity: 841, prominence: 4.2,
    sharpness: 0.3, source: "auto", excluded: false },
];
const EXPOSURES: api.Exposure[] = [
  { id: 1, sample_id: 3, filename: "JC001-120.dat", kind: "file",
    selected: true, status: "accepted", image_path: null, image_version: "", tags: [], sources: [] },
];

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ activeSampleId: 3, activeExposureId: 1, activeExperimentId: 1 });
  vi.spyOn(api, "listIndices").mockResolvedValue(INDICES);
  vi.spyOn(api, "listPeaks").mockResolvedValue(PEAKS);
  vi.spyOn(api, "listExposures").mockResolvedValue(EXPOSURES);
  vi.spyOn(api, "listSamples").mockResolvedValue([]);
});

describe("<MentionPicker>", () => {
  it("shows indices section when query matches phase", async () => {
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/Pn3m/)).toBeInTheDocument();
  });

  it("shows peak section when query matches q value", async () => {
    renderWithProviders(
      <MentionPicker query="1.22" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/1\.223/)).toBeInTheDocument();
  });

  it("calls onSelect with [[type:id]] token when row is clicked", async () => {
    const onSelect = vi.fn();
    const user = userEvent.setup();
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={onSelect} onDismiss={vi.fn()} />
    );
    await user.click(await screen.findByText(/Pn3m/));
    expect(onSelect).toHaveBeenCalledWith("[[index:17]]");
  });

  it("calls onDismiss on Escape key", async () => {
    const onDismiss = vi.fn();
    const user = userEvent.setup();
    renderWithProviders(
      <MentionPicker query="Pn3" onSelect={vi.fn()} onDismiss={onDismiss} />
    );
    await user.keyboard("{Escape}");
    expect(onDismiss).toHaveBeenCalled();
  });

  it("shows empty state when no results match", async () => {
    renderWithProviders(
      <MentionPicker query="zzznothing" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/no results/i)).toBeInTheDocument();
  });

  it("renders rows in a scrollable container with a footer hint when many matches", async () => {
    // 12 peaks + 1 default index + 1 default exposure = 14 total — over the 5-row cap.
    const manyPeaks: api.Peak[] = Array.from({ length: 12 }, (_, i) => ({
      id: 100 + i, exposure_id: 1, q: 1.0 + i * 0.01, intensity: 100,
      prominence: 1.0, sharpness: 0.1, source: "auto", excluded: false,
    }));
    vi.spyOn(api, "listPeaks").mockResolvedValue(manyPeaks);
    const { container } = renderWithProviders(
      <MentionPicker query="" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    expect(await screen.findByText(/14 matches.*scroll or refine/i)).toBeInTheDocument();
    const scroller = container.querySelector("[role='listbox'] .overflow-y-auto");
    expect(scroller).not.toBeNull();
  });

  it("does NOT show the footer hint when results fit (≤5)", async () => {
    // 1 index + 1 peak + 1 exposure = 3 total — under the cap.
    const { queryByText } = renderWithProviders(
      <MentionPicker query="" onSelect={vi.fn()} onDismiss={vi.fn()} />
    );
    // findByText would wait; queryByText returns immediately.
    await screen.findByText(/Pn3m/);
    expect(queryByText(/scroll or refine/i)).toBeNull();
  });

  // ─── Comparisons group (Phase 10) ──────────────────────────────────────
  describe("comparison rows", () => {
    const COMPARISONS: api.ComparisonSummary[] = [
      {
        id: 7, title: "DOPE titration",
        description: "5 candidates",
        content_hash: "a1b2c3d4e5f60718",
        created_by: 1,
        created_at: "2026-05-02 10:00:00",
        updated_at: "2026-05-02 10:00:00",
        forked_from_id: null,
        forked_at_hash: null,
      },
      {
        id: 8, title: "MO37 vs Im3m",
        description: null,
        content_hash: "deadbeef00112233",
        created_by: 1,
        created_at: "2026-05-03 09:00:00",
        updated_at: "2026-05-03 09:00:00",
        forked_from_id: null,
        forked_at_hash: null,
      },
    ];

    beforeEach(() => {
      vi.spyOn(api, "listExperimentComparisons").mockResolvedValue(COMPARISONS);
    });

    it("shows comparison rows when query matches a title", async () => {
      renderWithProviders(
        <MentionPicker query="DOPE" onSelect={vi.fn()} onDismiss={vi.fn()} />,
      );
      expect(await screen.findByText(/DOPE titration/)).toBeInTheDocument();
    });

    it("calls onSelect with eager-hash token for a comparison", async () => {
      const onSelect = vi.fn();
      const user = userEvent.setup();
      renderWithProviders(
        <MentionPicker query="DOPE" onSelect={onSelect} onDismiss={vi.fn()} />,
      );
      await user.click(await screen.findByText(/DOPE titration/));
      // hash8 = first 8 chars of content_hash, lowercased hex.
      expect(onSelect).toHaveBeenCalledWith("[[comparison:7@a1b2c3d4]]");
    });

    it("filters by `comparison:` type prefix", async () => {
      renderWithProviders(
        <MentionPicker
          query="comparison:Im3m"
          onSelect={vi.fn()}
          onDismiss={vi.fn()}
        />,
      );
      // Only matches MO37 vs Im3m (title contains "Im3m"); nothing else.
      expect(await screen.findByText(/MO37 vs Im3m/)).toBeInTheDocument();
      expect(screen.queryByText(/DOPE titration/)).toBeNull();
    });
  });
});
