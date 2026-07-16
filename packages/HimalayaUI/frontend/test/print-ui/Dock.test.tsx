import { describe, it, expect, vi, beforeEach } from "vitest";
import { render, screen } from "@testing-library/react";
import { Dock } from "../../src/print/ui/Dock";

// Mock floatingDock so the lane-claim effect doesn't touch Zustand.
const setCenterLaneOccupied = vi.fn();
vi.mock("../../src/print/shell/floatingDock", () => ({
  useFloatingDock: (sel: (s: { setCenterLaneOccupied: () => void }) => unknown) =>
    sel({ setCenterLaneOccupied }),
}));

beforeEach(() => {
  vi.clearAllMocks();
});

describe("Dock primitive", () => {
  it("renders children inside a fixed bar", () => {
    render(<Dock><button>seg</button></Dock>);
    expect(screen.getByRole("button", { name: "seg" })).toBeInTheDocument();
    expect(screen.getByTestId("dock")).toBeInTheDocument();
  });

  it("claims the center lane on mount", () => {
    render(<Dock><span>x</span></Dock>);
    expect(setCenterLaneOccupied).toHaveBeenCalledWith(true);
  });

  it("releases the center lane on unmount", () => {
    const { unmount } = render(<Dock><span>x</span></Dock>);
    unmount();
    expect(setCenterLaneOccupied).toHaveBeenCalledWith(false);
  });

  it("applies placement-only className after the base classes", () => {
    render(<Dock className="mt-4"><span>x</span></Dock>);
    expect(screen.getByTestId("dock").className).toContain("mt-4");
  });
});
