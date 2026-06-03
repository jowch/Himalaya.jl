import { describe, it, expect } from "vitest";
import { render, screen } from "@testing-library/react";
import { BigFrame } from "../../src/print/components/BigFrame";

describe("BigFrame", () => {
  it("renders the frame with a detector-image placeholder when src is null", () => {
    render(<BigFrame src={null} />);
    expect(screen.getByTestId("big-frame")).toBeInTheDocument();
    expect(screen.getByTestId("detector-image-placeholder")).toBeInTheDocument();
  });

  it("renders the mono caption when given", () => {
    render(<BigFrame src={null} caption="frame 65 · 0.40 s" />);
    const cap = screen.getByText("frame 65 · 0.40 s");
    expect(cap).toHaveAttribute("data-role", "frame-caption");
  });

  it("kept (default): no reject overlay, no dropped tag, data-rejected absent", () => {
    render(<BigFrame src={null} caption="f1" />);
    expect(screen.queryByTestId("reject-overlay")).not.toBeInTheDocument();
    expect(screen.queryByText("Dropped")).not.toBeInTheDocument();
    expect(screen.getByTestId("big-frame")).not.toHaveAttribute("data-rejected", "true");
  });

  it("rejected: shows the dropped tag and the reject overlay, flags data-rejected", () => {
    render(<BigFrame src={null} caption="f1" rejected />);
    expect(screen.getByTestId("reject-overlay")).toBeInTheDocument();
    const tag = screen.getByText("Dropped");
    expect(tag).toHaveAttribute("data-role", "dropped-tag");
    expect(screen.getByTestId("big-frame")).toHaveAttribute("data-rejected", "true");
  });
});
