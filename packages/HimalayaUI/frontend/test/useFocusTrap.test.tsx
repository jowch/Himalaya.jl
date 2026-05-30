import { describe, it, expect } from "vitest";
import { render, fireEvent } from "@testing-library/react";
import { useRef } from "react";
import { useFocusTrap } from "../src/hooks/useFocusTrap";

function TrapContainer({ active = true }: { active?: boolean }): JSX.Element {
  const ref = useRef<HTMLDivElement>(null);
  useFocusTrap(ref, active);
  return (
    <div ref={ref}>
      <button data-testid="btn-a">A</button>
      <button data-testid="btn-b">B</button>
      <button data-testid="btn-c">C</button>
    </div>
  );
}

function TextareaTrap({ active = true }: { active?: boolean }): JSX.Element {
  const ref = useRef<HTMLDivElement>(null);
  useFocusTrap(ref, active);
  return (
    <div ref={ref}>
      <textarea data-testid="ta-only" defaultValue="" />
    </div>
  );
}

describe("useFocusTrap", () => {
  it("Tab on the last element wraps focus to the first", () => {
    const { getByTestId } = render(<TrapContainer />);
    const a = getByTestId("btn-a");
    const c = getByTestId("btn-c");

    c.focus();
    expect(document.activeElement).toBe(c);

    fireEvent.keyDown(c, { key: "Tab", bubbles: true });
    expect(document.activeElement).toBe(a);
  });

  it("Shift+Tab on the first element wraps focus to the last", () => {
    const { getByTestId } = render(<TrapContainer />);
    const a = getByTestId("btn-a");
    const c = getByTestId("btn-c");

    a.focus();
    expect(document.activeElement).toBe(a);

    fireEvent.keyDown(a, { key: "Tab", shiftKey: true, bubbles: true });
    expect(document.activeElement).toBe(c);
  });

  it("Tab on a non-boundary element does not move focus", () => {
    const { getByTestId } = render(<TrapContainer />);
    const b = getByTestId("btn-b");

    b.focus();
    expect(document.activeElement).toBe(b);

    // Keydown fires but handler only acts on first/last — focus stays on b.
    fireEvent.keyDown(b, { key: "Tab", bubbles: true });
    expect(document.activeElement).toBe(b);
  });

  it("does not intercept Tab when inactive", () => {
    const { getByTestId } = render(<TrapContainer active={false} />);
    const c = getByTestId("btn-c");

    c.focus();
    const evt = fireEvent.keyDown(c, { key: "Tab", bubbles: true });
    // No handler registered — event should not be prevented.
    expect(evt).toBe(true); // fireEvent returns true when defaultPrevented is false
    expect(document.activeElement).toBe(c);
  });

  it("restores the previously-focused element on cleanup", () => {
    const trigger = document.createElement("button");
    document.body.appendChild(trigger);
    trigger.focus();
    expect(document.activeElement).toBe(trigger);

    const { getByTestId, unmount } = render(<TrapContainer />);
    getByTestId("btn-a").focus();

    unmount();
    expect(document.activeElement).toBe(trigger);

    document.body.removeChild(trigger);
  });

  it("traps a textarea-only container (Tab on the sole textarea stays put)", () => {
    const { getByTestId } = render(<TextareaTrap />);
    const ta = getByTestId("ta-only");
    ta.focus();
    expect(document.activeElement).toBe(ta);
    // Sole focusable is both first and last → Tab wraps to itself, focus held.
    const evt = fireEvent.keyDown(ta, { key: "Tab", bubbles: true });
    expect(evt).toBe(false); // defaultPrevented: trap acted (textarea is now in FOCUSABLE)
    expect(document.activeElement).toBe(ta);
  });

  it("includes anchors with href in the focusable set", () => {
    function AnchorTrap(): JSX.Element {
      const ref = useRef<HTMLDivElement>(null);
      useFocusTrap(ref, true);
      return (
        <div ref={ref}>
          <a href="#x" data-testid="lnk">link</a>
          <button data-testid="btn-z">Z</button>
        </div>
      );
    }
    const { getByTestId } = render(<AnchorTrap />);
    const lnk = getByTestId("lnk");
    const z = getByTestId("btn-z");
    z.focus();
    fireEvent.keyDown(z, { key: "Tab", bubbles: true }); // last → wraps to first (anchor)
    expect(document.activeElement).toBe(lnk);
  });
});
