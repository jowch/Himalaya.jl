import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, waitFor } from "@testing-library/react";
import userEvent from "@testing-library/user-event";
import { renderWithProviders } from "./test-utils";
import { MentionCompose } from "../src/components/MentionCompose";
import { useAppState } from "../src/state";
import * as api from "../src/api";

beforeEach(() => {
  localStorage.clear();
  useAppState.setState({ activeSampleId: 3, activeExposureId: 1, activeExperimentId: 1 });
  vi.spyOn(api, "listIndices").mockResolvedValue([]);
  vi.spyOn(api, "listPeaks").mockResolvedValue([]);
  vi.spyOn(api, "listExposures").mockResolvedValue([]);
  vi.spyOn(api, "listSamples").mockResolvedValue([]);
});

describe("<MentionCompose>", () => {
  it("shows picker when @ is typed", async () => {
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={vi.fn()} />);
    const ta = screen.getByRole("textbox");
    await user.click(ta);
    await user.keyboard("@");
    expect(await screen.findByRole("listbox")).toBeInTheDocument();
  });

  it("dismisses picker on Escape", async () => {
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={vi.fn()} />);
    const ta = screen.getByRole("textbox");
    await user.click(ta);
    await user.keyboard("@foo");
    await screen.findByRole("listbox");
    await user.keyboard("{Escape}");
    await waitFor(() => {
      expect(screen.queryByRole("listbox")).toBeNull();
    });
  });

  it("inserts [[type:id]] token and closes picker on selection", async () => {
    vi.spyOn(api, "listIndices").mockResolvedValue([
      { id: 17, exposure_id: 1, phase: "Pn3m", basis: 1.0, score: 0.91,
        r_squared: 0.998, lattice_d: 5.14, ngc: -1.51, status: "candidate", peaks: [],
        predicted_q: [1.223] },
    ]);
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={vi.fn()} />);
    const ta = screen.getByRole("textbox") as HTMLTextAreaElement;
    await user.click(ta);
    await user.keyboard("see @Pn3");
    const row = await screen.findByText(/Pn3m/);
    await user.click(row);
    expect(ta.value).toBe("see [[index:17]]");
    expect(screen.queryByRole("listbox")).toBeNull();
  });

  it("submits the raw token body on Enter", async () => {
    const onSubmit = vi.fn();
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={onSubmit} />);
    const ta = screen.getByRole("textbox");
    await user.click(ta);
    await user.keyboard("hello");
    await user.keyboard("{Enter}");
    expect(onSubmit).toHaveBeenCalledWith("hello");
  });

  // Phase 10 — eager hash insertion for @comparison:N picks. The picker
  // tokenizes a comparison row as `[[comparison:N@<hash8>]]` (rowToken in
  // MentionPicker); MentionCompose only has to pass that string through.
  // This test is the contract that the compose path doesn't strip the
  // hash suffix.
  it("inserts eager-hash token for comparison selections", async () => {
    vi.spyOn(api, "listExperimentComparisons").mockResolvedValue([
      {
        id: 7, title: "DOPE titration",
        description: null,
        content_hash: "a1b2c3d4e5f60718",
        created_by: 1,
        created_at: "2026-05-02 10:00:00",
        updated_at: "2026-05-02 10:00:00",
        forked_from_id: null,
        forked_at_hash: null,
        view_grouping_mode: null,
        view_show_peak_ticks: null,
        view_show_peak_labels: null,
        last_event_at: null,
        author_username: null,
        member_count: 0,
        member_phases: [],
        member_phase_count: 0,
        has_stale_members: false,
      },
    ]);
    const user = userEvent.setup();
    renderWithProviders(<MentionCompose disabled={false} onSubmit={vi.fn()} />);
    const ta = screen.getByRole("textbox") as HTMLTextAreaElement;
    await user.click(ta);
    await user.keyboard("see @DOPE");
    const row = await screen.findByText(/DOPE titration/);
    await user.click(row);
    expect(ta.value).toBe("see [[comparison:7@a1b2c3d4]]");
  });
});
