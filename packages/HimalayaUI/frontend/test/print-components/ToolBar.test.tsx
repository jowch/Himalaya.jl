import { render, screen } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { ToolBar } from "../../src/print/components/ToolBar";

describe("<ToolBar>", () => {
  it("renders a role=toolbar element", () => {
    render(<ToolBar><button>Auto-fit</button></ToolBar>);
    expect(screen.getByRole("toolbar")).toBeInTheDocument();
  });

  it("renders children inside the toolbar", () => {
    render(
      <ToolBar>
        <button>Auto-fit</button>
        <button>+ Peak</button>
      </ToolBar>,
    );
    expect(screen.getByRole("button", { name: /auto-fit/i })).toBeInTheDocument();
    expect(screen.getByRole("button", { name: /\+ peak/i })).toBeInTheDocument();
  });

  it("applies a placement className to the toolbar element", () => {
    render(<ToolBar className="ml-auto"><button>Auto-fit</button></ToolBar>);
    const toolbar = screen.getByRole("toolbar");
    expect(toolbar.className).toContain("ml-auto");
  });
});
