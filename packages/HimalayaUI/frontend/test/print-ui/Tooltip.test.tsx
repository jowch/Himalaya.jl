import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect } from "vitest";
import { Tooltip } from "../../src/print/ui/Tooltip";

describe("<Tooltip>", () => {
  it("is not in the DOM at rest", () => {
    render(
      <Tooltip label="Reanalyze">
        <button>Run</button>
      </Tooltip>,
    );
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("appears on mouseEnter of the wrapper/trigger", () => {
    render(
      <Tooltip label="Reanalyze">
        <button>Run</button>
      </Tooltip>,
    );
    fireEvent.mouseEnter(screen.getByText("Run"));
    const tip = screen.getByRole("tooltip");
    expect(tip).toHaveTextContent("Reanalyze");
  });

  it("appears on focus", () => {
    render(
      <Tooltip label="Reanalyze">
        <button>Run</button>
      </Tooltip>,
    );
    fireEvent.focus(screen.getByText("Run"));
    expect(screen.getByRole("tooltip")).toBeInTheDocument();
  });

  it("hides on mouseLeave", () => {
    render(
      <Tooltip label="Reanalyze">
        <button>Run</button>
      </Tooltip>,
    );
    fireEvent.mouseEnter(screen.getByText("Run"));
    expect(screen.getByRole("tooltip")).toBeInTheDocument();
    fireEvent.mouseLeave(screen.getByText("Run"));
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("hides on blur", () => {
    render(
      <Tooltip label="Reanalyze">
        <button>Run</button>
      </Tooltip>,
    );
    fireEvent.focus(screen.getByText("Run"));
    expect(screen.getByRole("tooltip")).toBeInTheDocument();
    fireEvent.blur(screen.getByText("Run"));
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("has role=tooltip and renders the label", () => {
    render(
      <Tooltip label="Caption text">
        <button>Run</button>
      </Tooltip>,
    );
    fireEvent.mouseEnter(screen.getByText("Run"));
    const tip = screen.getByRole("tooltip");
    expect(tip.getAttribute("data-testid")).toBe("tooltip");
    expect(tip).toHaveTextContent("Caption text");
  });

  it("carries the side as data-side (default top)", () => {
    render(
      <Tooltip label="x" side="bottom">
        <button>Run</button>
      </Tooltip>,
    );
    fireEvent.mouseEnter(screen.getByText("Run"));
    expect(screen.getByRole("tooltip").getAttribute("data-side")).toBe("bottom");
  });

  it("wires aria-describedby on the trigger to the tooltip id when open", () => {
    render(
      <Tooltip label="Reanalyze">
        <button>Run</button>
      </Tooltip>,
    );
    const trigger = screen.getByText("Run");
    expect(trigger.getAttribute("aria-describedby")).toBeNull();
    fireEvent.mouseEnter(trigger);
    const tip = screen.getByRole("tooltip");
    expect(trigger.getAttribute("aria-describedby")).toBe(tip.id);
    expect(tip.id).toBeTruthy();
  });
});
