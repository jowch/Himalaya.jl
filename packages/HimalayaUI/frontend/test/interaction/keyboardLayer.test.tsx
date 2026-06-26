import { describe, it, expect, vi, afterEach } from "vitest";
import { render, fireEvent, cleanup } from "@testing-library/react";
import { useKeyboardLayer } from "../../src/print/interaction/useKeyboardLayer";
import { useInteraction } from "../../src/print/interaction/registry";
import { core, page } from "../../src/print/interaction/core";

function Harness(): JSX.Element {
  useKeyboardLayer();
  return <div data-interaction-scope tabIndex={-1} data-testid="scope" />;
}

afterEach(() => { useInteraction.getState().clearPage(); cleanup(); });

describe("useKeyboardLayer", () => {
  it("fires a matched action's run() and preventDefaults", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [page("cull", { label: "Cull", keys: ["x"], group: "Act", run })]);
    render(<Harness />);
    document.querySelector<HTMLElement>('[data-interaction-scope]')!.focus();
    const evt = fireEvent.keyDown(document.activeElement!, { key: "x" });
    expect(run).toHaveBeenCalledTimes(1);
    expect(evt).toBe(false); // preventDefault was called (jsdom convention)
  });

  it("respects enabled() === false (inert, no preventDefault)", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [
      page("cull", { label: "Cull", keys: ["x"], group: "Act", enabled: () => false, run }),
    ]);
    render(<Harness />);
    document.querySelector<HTMLElement>('[data-interaction-scope]')!.focus();
    const evt = fireEvent.keyDown(document.activeElement!, { key: "x" });
    expect(run).not.toHaveBeenCalled();
    expect(evt).toBe(true); // not prevented
  });

  it("bails when the event was already defaultPrevented (rung 1/2 hand-off)", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("back", { run })]); // Escape
    render(<Harness />);
    // Simulate a modal that already consumed Escape:
    const handled = new KeyboardEvent("keydown", { key: "Escape", cancelable: true, bubbles: true });
    handled.preventDefault();
    window.dispatchEvent(handled);
    expect(run).not.toHaveBeenCalled();
  });

  it("ignores a bare key while typing in an input", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [page("cull", { label: "Cull", keys: ["x"], group: "Act", run })]);
    render(<><Harness /><input data-testid="field" /></>);
    const input = document.querySelector<HTMLInputElement>('[data-testid="field"]')!;
    input.focus();
    fireEvent.keyDown(input, { key: "x" });
    expect(run).not.toHaveBeenCalled();
  });

  it("still fires a Mod-chord while typing (undo)", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("undo", { run })]);
    render(<><Harness /><input data-testid="field" /></>);
    const input = document.querySelector<HTMLInputElement>('[data-testid="field"]')!;
    input.focus();
    fireEvent.keyDown(input, { key: "z", metaKey: true });
    expect(run).toHaveBeenCalledTimes(1);
  });

  // Regression: Enter on a native interactive target (button, input) must NOT
  // fire the page's openFocus action — the control owns that Enter natively.
  it("Enter on a native button target does NOT fire the registered Enter action", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("openFocus", { run })]);
    render(
      <>
        <Harness />
        <div data-interaction-scope>
          <button data-testid="btn">Activate</button>
        </div>
      </>
    );
    const btn = document.querySelector<HTMLButtonElement>('[data-testid="btn"]')!;
    btn.focus();
    fireEvent.keyDown(btn, { key: "Enter" });
    expect(run).not.toHaveBeenCalled();
  });

  it("Enter on a native input target does NOT fire the registered Enter action", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("openFocus", { run })]);
    render(
      <>
        <Harness />
        <div data-interaction-scope>
          <input data-testid="field" />
        </div>
      </>
    );
    const input = document.querySelector<HTMLInputElement>('[data-testid="field"]')!;
    input.focus();
    fireEvent.keyDown(input, { key: "Enter" });
    expect(run).not.toHaveBeenCalled();
  });

  // Preservation: Enter on a non-native row (div[role="row"]) must STILL fire —
  // this is the Corpus/Loupe Enter-on-row navigation path.
  it("Enter on a non-native row element (div[role=row]) DOES fire the registered Enter action", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("openFocus", { run })]);
    render(
      <>
        <Harness />
        <div data-interaction-scope>
          <div role="row" data-testid="row" tabIndex={0}>Sample row</div>
        </div>
      </>
    );
    const row = document.querySelector<HTMLDivElement>('[data-testid="row"]')!;
    row.focus();
    fireEvent.keyDown(row, { key: "Enter" });
    expect(run).toHaveBeenCalledTimes(1);
  });
});
