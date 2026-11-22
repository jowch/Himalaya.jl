import { act, fireEvent, render, screen } from "@testing-library/react";
import { forwardRef, createRef } from "react";
import { afterEach, beforeEach, describe, expect, it, vi } from "vitest";
import { placeTooltip, resetTooltipWarmth, Tooltip, TOOLTIP_OPEN_DELAY_MS, TOOLTIP_WARM_WINDOW_MS } from "../../src/ui/Tooltip";

function hoverAndWait(el: HTMLElement): void {
  fireEvent.mouseEnter(el);
  act(() => {
    vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
  });
}

describe("<Tooltip>", () => {
  beforeEach(() => {
    vi.useFakeTimers();
    resetTooltipWarmth();
  });
  afterEach(() => {
    vi.useRealTimers();
  });

  it("is not in the DOM at rest", () => {
    render(<Tooltip content="Reanalyze"><button>Run</button></Tooltip>);
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("waits for the open delay before showing on hover", () => {
    render(<Tooltip content="Reanalyze"><button>Run</button></Tooltip>);
    fireEvent.mouseEnter(screen.getByText("Run"));
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS - 1);
    });
    expect(screen.queryByRole("tooltip")).toBeNull();
    act(() => {
      vi.advanceTimersByTime(1);
    });
    expect(screen.getByRole("tooltip")).toHaveTextContent("Reanalyze");
  });

  it("never opens when the pointer leaves before the delay", () => {
    render(<Tooltip content="Reanalyze"><button>Run</button></Tooltip>);
    fireEvent.mouseEnter(screen.getByText("Run"));
    fireEvent.mouseLeave(screen.getByText("Run"));
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS * 2);
    });
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("opens a neighbour instantly while warm, and waits again once cold", () => {
    render(
      <>
        <Tooltip content="First"><button>A</button></Tooltip>
        <Tooltip content="Second"><button>B</button></Tooltip>
      </>,
    );
    hoverAndWait(screen.getByText("A"));
    fireEvent.mouseLeave(screen.getByText("A"));
    fireEvent.mouseEnter(screen.getByText("B"));
    expect(screen.getByRole("tooltip")).toHaveTextContent("Second");
    fireEvent.mouseLeave(screen.getByText("B"));
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_WARM_WINDOW_MS);
    });
    fireEvent.mouseEnter(screen.getByText("A"));
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("shows at most one tooltip at a time", () => {
    render(
      <>
        <Tooltip content="First"><button>A</button></Tooltip>
        <Tooltip content="Second"><button>B</button></Tooltip>
      </>,
    );
    hoverAndWait(screen.getByText("A"));
    hoverAndWait(screen.getByText("B"));
    expect(screen.getAllByRole("tooltip")).toHaveLength(1);
    expect(screen.getByRole("tooltip")).toHaveTextContent("Second");
  });

  it("opens immediately on keyboard focus and closes on blur", () => {
    render(<Tooltip content="Reanalyze"><button>Run</button></Tooltip>);
    act(() => {
      screen.getByText("Run").focus();
    });
    expect(screen.getByRole("tooltip")).toBeInTheDocument();
    act(() => {
      screen.getByText("Run").blur();
    });
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("closes on Escape and stays closed until the pointer leaves", () => {
    render(<Tooltip content="Reanalyze"><button>Run</button></Tooltip>);
    const run = screen.getByText("Run");
    hoverAndWait(run);
    fireEvent.keyDown(window, { key: "Escape" });
    expect(screen.queryByRole("tooltip")).toBeNull();
    fireEvent.mouseEnter(run);
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_OPEN_DELAY_MS);
    });
    expect(screen.queryByRole("tooltip")).toBeNull();
    fireEvent.mouseLeave(run);
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_WARM_WINDOW_MS);
    });
    hoverAndWait(run);
    expect(screen.getByRole("tooltip")).toBeInTheDocument();
  });

  it("closes on pointer down and on scroll", () => {
    render(<Tooltip content="Reanalyze"><button>Run</button></Tooltip>);
    const run = screen.getByText("Run");
    hoverAndWait(run);
    fireEvent.pointerDown(run);
    expect(screen.queryByRole("tooltip")).toBeNull();
    fireEvent.mouseLeave(run);
    act(() => {
      vi.advanceTimersByTime(TOOLTIP_WARM_WINDOW_MS);
    });
    hoverAndWait(run);
    fireEvent.scroll(window);
    expect(screen.queryByRole("tooltip")).toBeNull();
  });

  it("adds no wrapper around the trigger and portals the tooltip to the body", () => {
    const { container } = render(<div data-testid="parent"><Tooltip content="Reanalyze"><button>Run</button></Tooltip></div>);
    const run = screen.getByText("Run");
    expect(run.parentElement).toBe(screen.getByTestId("parent"));
    hoverAndWait(run);
    const tip = screen.getByRole("tooltip");
    expect(container.contains(tip)).toBe(false);
    expect(tip.parentElement).toBe(document.body);
  });

  it("keeps the trigger's own handlers, ref and aria-describedby", () => {
    const onMouseEnter = vi.fn();
    const onFocus = vi.fn();
    const ref = createRef<HTMLButtonElement>();
    const Btn = forwardRef<HTMLButtonElement, React.ButtonHTMLAttributes<HTMLButtonElement>>((p, r) => <button ref={r} {...p} />);
    render(
      <>
        <span id="extra">extra</span>
        <Tooltip content="Reanalyze"><Btn ref={ref} aria-describedby="extra" onMouseEnter={onMouseEnter} onFocus={onFocus}>Run</Btn></Tooltip>
      </>,
    );
    const run = screen.getByText("Run");
    expect(ref.current).toBe(run);
    hoverAndWait(run);
    fireEvent.focus(run);
    expect(onMouseEnter).toHaveBeenCalledTimes(1);
    expect(onFocus).toHaveBeenCalledTimes(1);
    const tip = screen.getByRole("tooltip");
    expect(run.getAttribute("aria-describedby")).toBe(`extra ${tip.id}`);
  });

  it("only wires aria-describedby while open", () => {
    render(<Tooltip content="Reanalyze"><button>Run</button></Tooltip>);
    const run = screen.getByText("Run");
    expect(run.getAttribute("aria-describedby")).toBeNull();
    hoverAndWait(run);
    expect(run.getAttribute("aria-describedby")).toBe(screen.getByRole("tooltip").id);
  });

  it("renders a keycap for the shortcut", () => {
    render(<Tooltip content="Keyboard shortcuts" shortcut="?"><button>Help</button></Tooltip>);
    hoverAndWait(screen.getByText("Help"));
    const tip = screen.getByRole("tooltip");
    expect(tip).toHaveTextContent("Keyboard shortcuts");
    expect(tip.querySelector("kbd")).toHaveTextContent("?");
  });

  it("with overflowOnly, shows only when the trigger's text is truncated", () => {
    render(<Tooltip content="A very long sample name" overflowOnly><span>A very long sample name</span></Tooltip>);
    const name = screen.getByText("A very long sample name");
    hoverAndWait(name);
    expect(screen.queryByRole("tooltip")).toBeNull();
    fireEvent.mouseLeave(name);
    Object.defineProperty(name, "scrollWidth", { configurable: true, value: 200 });
    Object.defineProperty(name, "clientWidth", { configurable: true, value: 80 });
    hoverAndWait(name);
    expect(screen.getByRole("tooltip")).toBeInTheDocument();
  });

  it("renders the bare trigger when there is no content", () => {
    render(<Tooltip content={undefined}><button>Run</button></Tooltip>);
    hoverAndWait(screen.getByText("Run"));
    expect(screen.queryByRole("tooltip")).toBeNull();
  });
});

describe("placeTooltip", () => {
  const viewport = { width: 1000, height: 800 };
  const tip = { width: 120, height: 24 };

  it("centres above the trigger by default", () => {
    const p = placeTooltip({ top: 400, bottom: 420, left: 500, width: 40 }, tip, viewport, "top");
    expect(p).toEqual({ side: "top", top: 400 - 6 - 24, left: 520 - 60 });
  });

  it("flips below when there is no room above", () => {
    const p = placeTooltip({ top: 10, bottom: 30, left: 500, width: 40 }, tip, viewport, "top");
    expect(p.side).toBe("bottom");
    expect(p.top).toBe(36);
  });

  it("flips above when bottom is preferred but there is no room below", () => {
    const p = placeTooltip({ top: 770, bottom: 790, left: 500, width: 40 }, tip, viewport, "bottom");
    expect(p.side).toBe("top");
  });

  it("clamps to the viewport edges", () => {
    expect(placeTooltip({ top: 400, bottom: 420, left: 0, width: 20 }, tip, viewport, "top").left).toBe(8);
    expect(placeTooltip({ top: 400, bottom: 420, left: 990, width: 10 }, tip, viewport, "top").left).toBe(1000 - 8 - 120);
  });
});
