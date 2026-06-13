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

  describe("value correction (SC-VALUECORRECT)", () => {
    it("without onEditValue: no edit affordance, value stays a skip toggle only", () => {
      render(<ScopeSampleRow name="HEPES Only" sampleId="smp_42" trace={TRACE} value="1 : 0.25" />);
      expect(screen.queryByRole("button", { name: /edit value/i })).toBeNull();
      expect(screen.getByTestId("flag-button")).toHaveTextContent("1 : 0.25");
    });

    it("with onEditValue: a pencil 'Edit value for {name}' control renders beside the value", () => {
      render(
        <ScopeSampleRow
          name="HEPES Only"
          sampleId="smp_42"
          trace={TRACE}
          value="1 : 0.25"
          onEditValue={() => {}}
        />,
      );
      // Skip toggle and the new edit affordance are DISTINCT controls.
      expect(screen.getByTestId("flag-button")).toBeInTheDocument();
      expect(screen.getByRole("button", { name: /^edit value for HEPES Only$/i })).toBeInTheDocument();
    });

    it("clicking the pencil swaps the value into a focused input and hides the skip toggle", () => {
      render(
        <ScopeSampleRow
          name="HEPES Only"
          sampleId="smp_42"
          trace={TRACE}
          value="1 : 0.25"
          onEditValue={() => {}}
        />,
      );
      fireEvent.click(screen.getByRole("button", { name: /^edit value for HEPES Only$/i }));
      const input = screen.getByRole("textbox", { name: /^corrected value for HEPES Only$/i });
      expect(input).toHaveValue("1 : 0.25");
      expect(input).toHaveFocus();
      // While editing, the skip toggle is replaced by the correction field.
      expect(screen.queryByTestId("flag-button")).toBeNull();
    });

    it("Enter commits the corrected value and returns to the skip-toggle display", () => {
      const onEditValue = vi.fn();
      render(
        <ScopeSampleRow
          name="HEPES Only"
          sampleId="smp_42"
          trace={TRACE}
          value="10.25"
          onEditValue={onEditValue}
        />,
      );
      fireEvent.click(screen.getByRole("button", { name: /^edit value/i }));
      const input = screen.getByRole("textbox", { name: /corrected value/i });
      fireEvent.change(input, { target: { value: "1 : 0.25" } });
      fireEvent.keyDown(input, { key: "Enter" });
      expect(onEditValue).toHaveBeenCalledTimes(1);
      expect(onEditValue).toHaveBeenCalledWith("1 : 0.25");
      // Back to display: the flag toggle returns, the input is gone.
      expect(screen.getByTestId("flag-button")).toBeInTheDocument();
      expect(screen.queryByRole("textbox", { name: /corrected value/i })).toBeNull();
    });

    it("blur commits the corrected value", () => {
      const onEditValue = vi.fn();
      render(
        <ScopeSampleRow name="x" sampleId="smp_1" trace={TRACE} value="10.25" onEditValue={onEditValue} />,
      );
      fireEvent.click(screen.getByRole("button", { name: /^edit value/i }));
      const input = screen.getByRole("textbox", { name: /corrected value/i });
      fireEvent.change(input, { target: { value: "1 : 0.25" } });
      fireEvent.blur(input);
      expect(onEditValue).toHaveBeenCalledTimes(1);
      expect(onEditValue).toHaveBeenCalledWith("1 : 0.25");
    });

    it("Escape cancels: onEditValue is not called and the original value is restored", () => {
      const onEditValue = vi.fn();
      render(
        <ScopeSampleRow name="x" sampleId="smp_1" trace={TRACE} value="10.25" onEditValue={onEditValue} />,
      );
      fireEvent.click(screen.getByRole("button", { name: /^edit value/i }));
      const input = screen.getByRole("textbox", { name: /corrected value/i });
      fireEvent.change(input, { target: { value: "999" } });
      fireEvent.keyDown(input, { key: "Escape" });
      expect(onEditValue).not.toHaveBeenCalled();
      expect(screen.getByTestId("flag-button")).toHaveTextContent("10.25");
    });

    it("committing an unchanged value does not fire onEditValue (no spurious history)", () => {
      const onEditValue = vi.fn();
      render(
        <ScopeSampleRow name="x" sampleId="smp_1" trace={TRACE} value="1 : 0.25" onEditValue={onEditValue} />,
      );
      fireEvent.click(screen.getByRole("button", { name: /^edit value/i }));
      fireEvent.keyDown(screen.getByRole("textbox", { name: /corrected value/i }), { key: "Enter" });
      expect(onEditValue).not.toHaveBeenCalled();
    });

    it("an emptied value is not committed (the write must never carry value:'')", () => {
      const onEditValue = vi.fn();
      render(
        <ScopeSampleRow name="x" sampleId="smp_1" trace={TRACE} value="1 : 0.25" onEditValue={onEditValue} />,
      );
      fireEvent.click(screen.getByRole("button", { name: /^edit value/i }));
      const input = screen.getByRole("textbox", { name: /corrected value/i });
      fireEvent.change(input, { target: { value: "   " } });
      fireEvent.keyDown(input, { key: "Enter" });
      expect(onEditValue).not.toHaveBeenCalled();
      // Reverts to the original display.
      expect(screen.getByTestId("flag-button")).toHaveTextContent("1 : 0.25");
    });
  });
});
