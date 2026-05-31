import { render, screen } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ToggleSwitch } from "../../src/print/ui/ToggleSwitch";

describe("<ToggleSwitch>", () => {
  it("has role=switch with aria-checked reflecting checked=false", () => {
    render(<ToggleSwitch checked={false} onChange={() => {}} label="Heatmap" />);
    const sw = screen.getByRole("switch");
    expect(sw.getAttribute("aria-checked")).toBe("false");
  });

  it("has role=switch with aria-checked reflecting checked=true", () => {
    render(<ToggleSwitch checked={true} onChange={() => {}} label="Heatmap" />);
    const sw = screen.getByRole("switch");
    expect(sw.getAttribute("aria-checked")).toBe("true");
  });

  it("fires onChange(true) when clicked while off", () => {
    const onChange = vi.fn();
    render(<ToggleSwitch checked={false} onChange={onChange} label="Heatmap" />);
    screen.getByRole("switch").click();
    expect(onChange).toHaveBeenCalledWith(true);
  });

  it("fires onChange(false) when clicked while on", () => {
    const onChange = vi.fn();
    render(<ToggleSwitch checked={true} onChange={onChange} label="Heatmap" />);
    screen.getByRole("switch").click();
    expect(onChange).toHaveBeenCalledWith(false);
  });

  it("renders the label text as the switch accessible name", () => {
    render(<ToggleSwitch checked={false} onChange={() => {}} label="Heatmap" />);
    expect(screen.getByRole("switch", { name: "Heatmap" })).toBeInTheDocument();
    expect(screen.getByText("Heatmap")).toBeInTheDocument();
  });

  it("keeps the accessible name but visually hides the label text when hideLabel", () => {
    render(<ToggleSwitch checked={false} onChange={() => {}} label="Heatmap" hideLabel />);
    // accessible name preserved (aria-label on the switch)
    expect(screen.getByRole("switch", { name: "Heatmap" })).toBeInTheDocument();
    // the visible label span carries a data-role we can find, and is marked hidden
    const labelSpan = screen.getByTestId("toggle-switch").querySelector('[data-role="switch-label"]');
    expect(labelSpan).not.toBeNull();
    expect(labelSpan?.getAttribute("data-hidden")).toBe("true");
  });

  it("shows the label text by default", () => {
    render(<ToggleSwitch checked={false} onChange={() => {}} label="Heatmap" />);
    const labelSpan = screen.getByTestId("toggle-switch").querySelector('[data-role="switch-label"]');
    expect(labelSpan?.getAttribute("data-hidden")).toBe(null);
  });
});
