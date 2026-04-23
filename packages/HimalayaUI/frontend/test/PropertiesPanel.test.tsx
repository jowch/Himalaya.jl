import { describe, it, expect, vi, beforeEach } from "vitest";
import { screen, fireEvent } from "@testing-library/react";
import { renderWithProviders } from "./test-utils";
import { PropertiesPanel } from "../src/components/PropertiesPanel";
import { useAppState } from "../src/state";

beforeEach(() => {
  vi.restoreAllMocks();
  useAppState.setState({
    activeSampleId: undefined,
    activeExposureId: undefined,
  });
  vi.spyOn(global, "fetch").mockResolvedValue(
    new Response("[]", { status: 200, headers: { "Content-Type": "application/json" } }),
  );
});

describe("<PropertiesPanel>", () => {
  it("renders four tab buttons", () => {
    renderWithProviders(<PropertiesPanel />);
    expect(screen.getByRole("tab", { name: /exposures/i })).toBeInTheDocument();
    expect(screen.getByRole("tab", { name: /peaks/i })).toBeInTheDocument();
    expect(screen.getByRole("tab", { name: /tags/i })).toBeInTheDocument();
    expect(screen.getByRole("tab", { name: /notes/i })).toBeInTheDocument();
  });

  it("starts on Exposures tab", () => {
    renderWithProviders(<PropertiesPanel />);
    expect(screen.getByRole("tab", { name: /exposures/i }))
      .toHaveAttribute("aria-selected", "true");
  });

  it("switching tabs updates aria-selected", () => {
    renderWithProviders(<PropertiesPanel />);
    fireEvent.click(screen.getByRole("tab", { name: /peaks/i }));
    expect(screen.getByRole("tab", { name: /peaks/i }))
      .toHaveAttribute("aria-selected", "true");
    expect(screen.getByRole("tab", { name: /exposures/i }))
      .toHaveAttribute("aria-selected", "false");
  });

  it("tab panel has matching aria-labelledby for the active tab", () => {
    renderWithProviders(<PropertiesPanel />);
    fireEvent.click(screen.getByRole("tab", { name: /tags/i }));
    const panel = screen.getByRole("tabpanel");
    const activeTab = screen.getByRole("tab", { name: /tags/i });
    expect(panel.getAttribute("aria-labelledby")).toBe(activeTab.id);
  });
});
