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
});
