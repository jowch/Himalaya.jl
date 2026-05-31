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
});
