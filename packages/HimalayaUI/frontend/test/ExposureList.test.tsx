import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent, waitFor } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { ExposureList } from "../src/components/ExposureList";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({ activeSampleId: undefined, activeExposureId: undefined });
});

function mockExposures(items: unknown[]): void {
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response(JSON.stringify(items), {
      status: 200, headers: { "Content-Type": "application/json" },
    }),
  );
}

describe("<ExposureList>", () => {
  it("shows a hint when no sample is active", () => {
    renderWithProviders(<ExposureList />);
    expect(screen.getByText(/no sample selected/i)).toBeInTheDocument();
  });

  it("renders one row per exposure and marks the active one", async () => {
    useAppState.setState({ activeSampleId: 1, activeExposureId: 11 });
    mockExposures([
      { id: 10, sample_id: 1, filename: "a", kind: "file",
        selected: false, tags: [], sources: [] },
      { id: 11, sample_id: 1, filename: "b", kind: "file",
        selected: true,  tags: [], sources: [] },
    ]);
    renderWithProviders(<ExposureList />);
    await waitFor(() => expect(screen.getByText("a")).toBeInTheDocument());
    expect(screen.getByText("b")).toBeInTheDocument();
    const active = document.querySelector('[data-exposure-id="11"][data-active="true"]');
    expect(active).not.toBeNull();
  });

  it("clicking a row updates activeExposureId in Zustand", async () => {
    useAppState.setState({ activeSampleId: 1, activeExposureId: undefined });
    mockExposures([
      { id: 10, sample_id: 1, filename: "a", kind: "file",
        selected: false, tags: [], sources: [] },
    ]);
    renderWithProviders(<ExposureList />);
    await waitFor(() => expect(screen.getByText("a")).toBeInTheDocument());
    const row = document.querySelector<HTMLElement>('[data-exposure-id="10"]')!;
    fireEvent.click(row);
    expect(useAppState.getState().activeExposureId).toBe(10);
  });
});
