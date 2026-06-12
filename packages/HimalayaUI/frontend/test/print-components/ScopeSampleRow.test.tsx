import { render, screen, fireEvent } from "@testing-library/react";
import { describe, it, expect, vi } from "vitest";
import { ScopeSampleRow } from "../../src/print/components/ScopeSampleRow";

const TRACE = { q: [0.02, 0.04, 0.08, 0.16], I: [10, 80, 30, 5] };

describe("<ScopeSampleRow>", () => {
  it("renders name, id, sparkline, grip, and the value", () => {
    render(<ScopeSampleRow name="Lipid 1-2 + LL37 1:0.25" sampleId="smp_07" trace={TRACE} phase="Pn3m" value="1 : 0.25" />);
    const row = screen.getByTestId("scope-sample-row");
    expect(row).toHaveTextContent("Lipid 1-2 + LL37 1:0.25");
    expect(row).toHaveTextContent("smp_07");
    expect(screen.getByTestId("sparkline")).toBeInTheDocument();
    expect(screen.getByTestId("grip-handle")).toBeInTheDocument();
    expect(screen.getByTestId("flag-button")).toHaveTextContent("1 : 0.25");
  });
  it("marks the flagged state on the row and the value", () => {
    render(<ScopeSampleRow name="Lipid 1-2, no LL37" sampleId="smp_04" trace={TRACE} phase="Pn3m" value="1 : 0" flagged />);
    expect(screen.getByTestId("scope-sample-row").getAttribute("data-flagged")).toBe("true");
    expect(screen.getByTestId("flag-button").getAttribute("data-state")).toBe("flagged");
  });
  it("forwards the flag toggle", () => {
    const onToggleFlag = vi.fn();
    render(<ScopeSampleRow name="x" sampleId="smp_04" trace={TRACE} value="1 : 0" flagged onToggleFlag={onToggleFlag} />);
    fireEvent.click(screen.getByTestId("flag-button"));
    expect(onToggleFlag).toHaveBeenCalledOnce();
  });

  describe("keyboard reorder contract (SC-KBD)", () => {
    it("with onMoveBy: the grip is a real button named 'Reorder {name}' (anchored)", () => {
      render(
        <ScopeSampleRow name="Lipid A" sampleId="smp_01" trace={TRACE} value="1 : 0" onMoveBy={() => {}} />,
      );
      expect(screen.getByRole("button", { name: /^reorder Lipid A$/i })).toBeInTheDocument();
      // The glyph itself is unchanged (still the aria-hidden visual handle).
      expect(screen.getByTestId("grip-handle")).toBeInTheDocument();
    });

    it("ArrowUp calls onMoveBy(-1) and prevents default", () => {
      const onMoveBy = vi.fn();
      render(
        <ScopeSampleRow name="Lipid A" sampleId="smp_01" trace={TRACE} value="1 : 0" onMoveBy={onMoveBy} />,
      );
      const btn = screen.getByRole("button", { name: /^reorder Lipid A$/i });
      // fireEvent returns false when a handler called preventDefault on a
      // cancelable event — the page must not scroll on a handled arrow key.
      const notPrevented = fireEvent.keyDown(btn, { key: "ArrowUp" });
      expect(onMoveBy).toHaveBeenCalledTimes(1);
      expect(onMoveBy).toHaveBeenCalledWith(-1);
      expect(notPrevented).toBe(false);
    });

    it("ArrowDown calls onMoveBy(1) and prevents default", () => {
      const onMoveBy = vi.fn();
      render(
        <ScopeSampleRow name="Lipid A" sampleId="smp_01" trace={TRACE} value="1 : 0" onMoveBy={onMoveBy} />,
      );
      const btn = screen.getByRole("button", { name: /^reorder Lipid A$/i });
      const notPrevented = fireEvent.keyDown(btn, { key: "ArrowDown" });
      expect(onMoveBy).toHaveBeenCalledTimes(1);
      expect(onMoveBy).toHaveBeenCalledWith(1);
      expect(notPrevented).toBe(false);
    });

    it("modified arrows stay native (no move, no preventDefault)", () => {
      // Cmd+ArrowDown is macOS scroll-to-end; Alt/Ctrl arrows are nav chords.
      // Hijacking them from a focused control would be its own keyboard trap.
      const onMoveBy = vi.fn();
      render(
        <ScopeSampleRow name="Lipid A" sampleId="smp_01" trace={TRACE} value="1 : 0" onMoveBy={onMoveBy} />,
      );
      const btn = screen.getByRole("button", { name: /^reorder Lipid A$/i });
      for (const mod of [{ metaKey: true }, { ctrlKey: true }, { altKey: true }]) {
        const notPrevented = fireEvent.keyDown(btn, { key: "ArrowDown", ...mod });
        expect(notPrevented).toBe(true);
      }
      expect(onMoveBy).not.toHaveBeenCalled();
    });

    it("other keys neither call onMoveBy nor prevent default", () => {
      const onMoveBy = vi.fn();
      render(
        <ScopeSampleRow name="Lipid A" sampleId="smp_01" trace={TRACE} value="1 : 0" onMoveBy={onMoveBy} />,
      );
      const btn = screen.getByRole("button", { name: /^reorder Lipid A$/i });
      const notPrevented = fireEvent.keyDown(btn, { key: "ArrowLeft" });
      expect(onMoveBy).not.toHaveBeenCalled();
      expect(notPrevented).toBe(true);
    });

    it("without onMoveBy: no reorder button renders, the visual grip remains", () => {
      render(<ScopeSampleRow name="Lipid A" sampleId="smp_01" trace={TRACE} value="1 : 0" />);
      expect(screen.queryByRole("button", { name: /reorder/i })).toBeNull();
      expect(screen.getByTestId("grip-handle")).toBeInTheDocument();
    });

    it("the reorder grip declares a >=24px hit area (WCAG 2.5.8)", () => {
      render(
        <ScopeSampleRow name="Lipid A" sampleId="smp_01" trace={TRACE} value="1 : 0" onMoveBy={() => {}} />,
      );
      // jsdom cannot compute Tailwind box sizes, so the padded-hit-area contract
      // is carried as a data attribute (Checkbox precedent).
      expect(
        screen.getByRole("button", { name: /^reorder Lipid A$/i }),
      ).toHaveAttribute("data-hit-area", "24");
    });

    it("surfaces the arrow-key reorder shortcut to sighted keyboard users (tooltip)", () => {
      render(
        <ScopeSampleRow name="Lipid A" sampleId="smp_01" trace={TRACE} value="1 : 0" onMoveBy={() => {}} />,
      );
      // The shortcut is otherwise only in aria-keyshortcuts (AT-only); a visible
      // title makes it discoverable.
      const grip = screen.getByRole("button", { name: /^reorder Lipid A$/i });
      expect(grip.getAttribute("title")).toMatch(/↑ ↓/);
      expect(grip).toHaveAttribute("aria-keyshortcuts", "ArrowUp ArrowDown");
    });
  });
});
