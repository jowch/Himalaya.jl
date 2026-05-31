import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { Chip } from "../../src/print/ui/Chip";

describe("<Chip> variants", () => {
  it('static is the default variant and sets data-variant="static"', () => {
    render(<Chip>Pn3m</Chip>);
    const el = screen.getByTestId("chip");
    expect(el.getAttribute("data-variant")).toBe("static");
    expect(el.textContent).toContain("Pn3m");
  });

  it('removable sets data-variant="removable"', () => {
    render(<Chip variant="removable">Pn3m</Chip>);
    expect(screen.getByTestId("chip").getAttribute("data-variant")).toBe(
      "removable",
    );
  });

  it('add sets data-variant="add"', () => {
    render(<Chip variant="add">+ add</Chip>);
    expect(screen.getByTestId("chip").getAttribute("data-variant")).toBe("add");
  });

  it('toggle sets data-variant="toggle"', () => {
    render(<Chip variant="toggle">Kept</Chip>);
    expect(screen.getByTestId("chip").getAttribute("data-variant")).toBe(
      "toggle",
    );
  });

  it('trigger sets data-variant="trigger"', () => {
    render(<Chip variant="trigger">Beamtime</Chip>);
    expect(screen.getByTestId("chip").getAttribute("data-variant")).toBe(
      "trigger",
    );
  });
});

describe("<Chip> testId override", () => {
  it('defaults data-testid to "chip"', () => {
    render(<Chip>x</Chip>);
    expect(screen.getByTestId("chip")).toBeTruthy();
  });

  it("overrides data-testid via the testId prop", () => {
    render(
      <Chip variant="toggle" testId="filter-chip">
        x
      </Chip>,
    );
    expect(screen.getByTestId("filter-chip")).toBeTruthy();
    expect(screen.queryByTestId("chip")).toBeNull();
  });

  it("override applies on every variant root (trigger)", () => {
    render(
      <Chip variant="trigger" testId="facet-chip">
        x
      </Chip>,
    );
    expect(screen.getByTestId("facet-chip")).toBeTruthy();
  });
});

describe("<Chip> removable", () => {
  it('exposes a remove control with aria-label="Remove"', () => {
    render(
      <Chip variant="removable" onRemove={() => {}}>
        Pn3m
      </Chip>,
    );
    expect(screen.getByRole("button", { name: "Remove" })).toBeTruthy();
  });

  it("fires onRemove when the × is clicked", () => {
    const onRemove = vi.fn();
    render(
      <Chip variant="removable" onRemove={onRemove}>
        Pn3m
      </Chip>,
    );
    fireEvent.click(screen.getByRole("button", { name: "Remove" }));
    expect(onRemove).toHaveBeenCalledTimes(1);
  });

  it("the × renders within the chip root", () => {
    render(
      <Chip variant="removable" onRemove={() => {}}>
        Pn3m
      </Chip>,
    );
    const root = screen.getByTestId("chip");
    const remove = screen.getByRole("button", { name: "Remove" });
    expect(root.contains(remove)).toBe(true);
  });
});

describe("<Chip> toggle", () => {
  it("reflects active on aria-pressed and data-active", () => {
    const { rerender } = render(
      <Chip variant="toggle" active={false}>
        Kept
      </Chip>,
    );
    let el = screen.getByTestId("chip");
    expect(el.getAttribute("aria-pressed")).toBe("false");
    expect(el.getAttribute("data-active")).toBe("false");

    rerender(
      <Chip variant="toggle" active>
        Kept
      </Chip>,
    );
    el = screen.getByTestId("chip");
    expect(el.getAttribute("aria-pressed")).toBe("true");
    expect(el.getAttribute("data-active")).toBe("true");
  });

  it("fires onClick when toggled", () => {
    const onClick = vi.fn();
    render(
      <Chip variant="toggle" onClick={onClick}>
        Kept
      </Chip>,
    );
    fireEvent.click(screen.getByTestId("chip"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });
});

describe("<Chip> trigger", () => {
  it('declares aria-haspopup="menu"', () => {
    render(<Chip variant="trigger">Beamtime</Chip>);
    expect(screen.getByTestId("chip").getAttribute("aria-haspopup")).toBe(
      "menu",
    );
  });

  it("renders an aria-hidden chevron", () => {
    render(<Chip variant="trigger">Beamtime</Chip>);
    const chevron = screen
      .getByTestId("chip")
      .querySelector('[aria-hidden="true"]');
    expect(chevron).toBeTruthy();
    expect(chevron?.textContent).toBe("▾");
  });

  it("fires onClick when clicked", () => {
    const onClick = vi.fn();
    render(
      <Chip variant="trigger" onClick={onClick}>
        Beamtime
      </Chip>,
    );
    fireEvent.click(screen.getByTestId("chip"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });
});

describe("<Chip> add", () => {
  it("fires onClick when clicked", () => {
    const onClick = vi.fn();
    render(
      <Chip variant="add" onClick={onClick}>
        + add
      </Chip>,
    );
    fireEvent.click(screen.getByTestId("chip"));
    expect(onClick).toHaveBeenCalledTimes(1);
  });

  it("defaults its content to + add", () => {
    render(<Chip variant="add" />);
    expect(screen.getByTestId("chip").textContent).toContain("+ add");
  });
});
