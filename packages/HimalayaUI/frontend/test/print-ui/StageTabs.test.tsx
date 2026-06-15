import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { StageTabs } from "../../src/print/ui/StageTabs";

describe("<StageTabs>", () => {
  it("renders all three tabs", () => {
    render(<StageTabs active="samples" onChange={() => {}} />);
    expect(screen.getByRole("tab", { name: /samples/i })).toBeInTheDocument();
    expect(screen.getByRole("tab", { name: /index/i })).toBeInTheDocument();
    expect(screen.getByRole("tab", { name: /series/i })).toBeInTheDocument();
  });

  it("marks the active tab via data-active and aria-selected", () => {
    render(<StageTabs active="index" onChange={() => {}} />);
    const indexTab = screen.getByRole("tab", { name: /index/i });
    expect(indexTab.getAttribute("data-active")).toBe("true");
    expect(indexTab.getAttribute("aria-selected")).toBe("true");
    const samplesTab = screen.getByRole("tab", { name: /samples/i });
    expect(samplesTab.getAttribute("data-active")).toBe("false");
  });

  it("calls onChange with the clicked stage key", () => {
    const onChange = vi.fn();
    render(<StageTabs active="samples" onChange={onChange} />);
    fireEvent.click(screen.getByRole("tab", { name: /series/i }));
    expect(onChange).toHaveBeenCalledWith("series");
  });

  it("tones the active tab's dot accent and inactive dots neutral", () => {
    render(<StageTabs active="index" onChange={() => {}} />);
    const indexTab = screen.getByRole("tab", { name: /index/i });
    expect(indexTab.querySelector("[data-tone]")?.getAttribute("data-tone")).toBe("accent");
    const samplesTab = screen.getByRole("tab", { name: /samples/i });
    expect(samplesTab.querySelector("[data-tone]")?.getAttribute("data-tone")).toBe("neutral");
  });
});
