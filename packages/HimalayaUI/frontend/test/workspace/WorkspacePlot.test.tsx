import { describe, it, expect, vi } from "vitest";
import { render, fireEvent, screen } from "@testing-library/react";
import { WorkspacePlot } from "../../src/workspace/ui/WorkspacePlot";

const trace = { q: [0.05, 0.1, 0.15, 0.2, 0.25], I: [10, 50, 8, 6, 5], sigma: [0, 0, 0, 0, 0] };

describe("WorkspacePlot", () => {
  it("double-click does not reset the zoom", () => {
    const onXDomain = vi.fn();
    render(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} xDomain={[0.1, 0.2]} onXDomain={onXDomain} />);
    fireEvent.doubleClick(screen.getByTestId("workspace-trace"));
    expect(onXDomain).not.toHaveBeenCalled();
  });

  it("skips non-positive comparison points on a log scale but draws them on linear", () => {
    const compare = { id: 3, name: "buf", q: [0.05, 0.1, 0.15], I: [1, 0, 2], dash: "6 4", emphasis: false };
    const { rerender } = render(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} compares={[compare]} />);
    const d = (): string => screen.getByTestId("workspace-compares").querySelector("path")!.getAttribute("d")!;
    expect(d().match(/[ML]/g)).toHaveLength(2);
    rerender(<WorkspacePlot trace={trace} peaks={[]} phase={null} combs={[]} compares={[compare]} logY={false} />);
    expect(d().match(/[ML]/g)).toHaveLength(3);
  });
});
