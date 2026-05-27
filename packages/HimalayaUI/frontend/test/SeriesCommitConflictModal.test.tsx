import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen, fireEvent, act } from "@testing-library/react";
import { MemoryRouter } from "react-router-dom";
import { SeriesCommitConflictModal } from "../src/components/SeriesCommitConflictModal";
import { useAppState } from "../src/state";
import { ConflictError } from "../src/api";
import type { Series, Comparison } from "../src/api";

const h = vi.hoisted(() => ({
  commit: {
    mutate: vi.fn(),
    isPending: false,
    isSuccess: false,
    error: null as unknown,
    data: undefined as unknown,
  },
  navigate: vi.fn(),
}));

vi.mock("../src/queries", () => ({
  useCommitSeriesPlate: () => h.commit,
}));
vi.mock("react-router-dom", async (orig) => {
  const actual = await orig<typeof import("react-router-dom")>();
  return { ...actual, useNavigate: () => h.navigate };
});

function serverSeries(id: number, hash = "sha256:server"): Series {
  return {
    id, title: "Server series", description: "srv desc", content_hash: hash,
    created_by: 1, created_at: null, updated_at: "2026-05-02",
    forked_from_id: null, forked_at_hash: null, forked_from_title: null,
    view_grouping_mode: null, view_show_peak_ticks: null,
    view_show_peak_labels: null, ordering_variable: null, order_rule: "manual",
    state: "committed",
    members: [{
      id: 1, series_id: id, exposure_id: 101, display_order: 0,
      band_height: 1, y_offset: 0, normalization: "max",
      color_override: null, label_override: null, q_window_min: null,
      q_window_max: null, peak_display: null, snapshot: null,
      is_stale: false, created_by: 1, created_at: null,
    }],
    samples: [],
  };
}

function comparison(id: number): Comparison {
  return {
    id, title: "C", description: null, content_hash: "sha256:c",
    created_by: 1, created_at: null, updated_at: null, forked_from_id: null,
    forked_at_hash: null, forked_from_title: null, view_grouping_mode: null,
    view_show_peak_ticks: null, view_show_peak_labels: null,
    last_event_at: null, members: [],
  };
}

function renderModal() {
  return render(
    <MemoryRouter initialEntries={["/series/5"]}>
      <SeriesCommitConflictModal />
    </MemoryRouter>,
  );
}

describe("SeriesCommitConflictModal", () => {
  beforeEach(() => {
    h.commit.mutate.mockReset();
    h.navigate.mockReset();
    h.commit.isPending = false;
    h.commit.isSuccess = false;
    h.commit.error = null;
    h.commit.data = undefined;
    sessionStorage.clear();
    act(() => {
      useAppState.getState().setPendingConflict(null);
      useAppState.getState().discardSeriesDraft();
    });
  });

  it("renders nothing when the conflict is null", () => {
    const { container } = renderModal();
    expect(container.firstChild).toBeNull();
  });

  it("renders nothing for a COMPARISON conflict (wrong entity kind)", () => {
    act(() => {
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:c", comparison(9)),
      );
    });
    renderModal();
    expect(screen.queryByTestId("conflict-modal")).toBeNull();
  });

  it("opens for a SERIES conflict and shows server + draft panels", () => {
    act(() => {
      useAppState.getState().startSeriesDraftFromSeries(serverSeries(5, "sha256:old"));
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", serverSeries(5)),
      );
    });
    renderModal();
    expect(screen.getByTestId("conflict-modal")).toBeInTheDocument();
    expect(screen.getByText("Series changed while you were editing")).toBeInTheDocument();
    expect(screen.getByTestId("conflict-panel-server-title")).toHaveTextContent("Server series");
    // No Fork action for series.
    expect(screen.queryByTestId("conflict-fork")).toBeNull();
  });

  it("Discard clears the slot + draft and navigates to /series/:id", () => {
    act(() => {
      useAppState.getState().startSeriesDraftFromSeries(serverSeries(5, "sha256:old"));
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", serverSeries(5)),
      );
    });
    renderModal();
    fireEvent.click(screen.getByTestId("conflict-discard"));
    expect(h.navigate).toHaveBeenCalledWith("/series/5");
    expect(useAppState.getState().pendingConflict).toBeNull();
    expect(useAppState.getState().seriesDraft).toBeNull();
  });

  it("Overwrite re-commits with expected_content_hash = server current_hash", () => {
    act(() => {
      useAppState.getState().startSeriesDraftFromSeries(serverSeries(5, "sha256:old"));
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", serverSeries(5)),
      );
    });
    renderModal();
    fireEvent.click(screen.getByTestId("conflict-overwrite"));
    expect(h.commit.mutate).toHaveBeenCalledTimes(1);
    const arg = h.commit.mutate.mock.calls[0][0];
    expect(arg.id).toBe(5);
    expect(arg.expected_content_hash).toBe("sha256:server");
    // Plate rebuilt from the SERVER's current members (id-less, positional).
    expect(arg.members).toHaveLength(1);
    expect(arg.members[0].exposure_id).toBe(101);
    expect("id" in arg.members[0]).toBe(false);
  });

  it("double-click guard: a second synchronous click does not re-fire mutate", () => {
    act(() => {
      useAppState.getState().startSeriesDraftFromSeries(serverSeries(5, "sha256:old"));
      useAppState.getState().setPendingConflict(
        new ConflictError("sha256:server", serverSeries(5)),
      );
    });
    renderModal();
    const btn = screen.getByTestId("conflict-overwrite");
    fireEvent.click(btn);
    fireEvent.click(btn);
    expect(h.commit.mutate).toHaveBeenCalledTimes(1);
  });
});
