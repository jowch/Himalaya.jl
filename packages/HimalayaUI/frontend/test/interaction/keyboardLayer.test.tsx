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

  // Regression: Mod+z inside a text input must NOT fire the page's undo action —
  // the native text-undo must reach the input's own edit history. (Previously the
  // guard only suppressed bare keys while typing; now it suppresses ALL shortcuts
  // so chorded gestures like Mod+z also reach native handlers.)
  it("Mod+z while typing in an input does NOT fire the registered undo action", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("undo", { run })]);
    render(<><Harness /><input data-testid="field" /></>);
    const input = document.querySelector<HTMLInputElement>('[data-testid="field"]')!;
    input.focus();
    fireEvent.keyDown(input, { key: "z", metaKey: true });
    expect(run).not.toHaveBeenCalled();
  });

  // Preservation: Mod+z outside a text field (e.g. on the scope container) DOES
  // fire the registered undo action — the widening must not suppress global-scope
  // chorded shortcuts.
  it("Mod+z on a non-typing element (scope container) DOES fire the registered undo action", () => {
    const run = vi.fn();
    useInteraction.getState().setPage(null, [core("undo", { run })]);
    render(<Harness />);
    const scope = document.querySelector<HTMLElement>('[data-interaction-scope]')!;
    scope.focus();
    fireEvent.keyDown(scope, { key: "z", metaKey: true });
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

  // Regression: a focus-trapped overlay ([role="dialog"]) owns the keyboard —
  // scope-exempt arrows must NOT navigate the surface behind an open modal
  // (ModalShell only preventDefaults Escape, never arrows).
  it("an arrow inside an open [role=dialog] does NOT fire the page arrowHandler", () => {
    const arrow = vi.fn();
    useInteraction.getState().setPage(null, [], [], null, arrow);
    render(
      <>
        <Harness />
        <div role="dialog" aria-modal="true">
          <button data-testid="dlg-btn">Cancel</button>
        </div>
      </>
    );
    const btn = document.querySelector<HTMLButtonElement>('[data-testid="dlg-btn"]')!;
    btn.focus();
    fireEvent.keyDown(btn, { key: "ArrowDown" });
    expect(arrow).not.toHaveBeenCalled();
  });

  // Preservation: the same arrow OUTSIDE any overlay (on the scope) DOES fire.
  it("an arrow on the scope container DOES fire the page arrowHandler", () => {
    const arrow = vi.fn();
    useInteraction.getState().setPage(null, [], [], null, arrow);
    render(<Harness />);
    const scope = document.querySelector<HTMLElement>('[data-interaction-scope]')!;
    scope.focus();
    fireEvent.keyDown(scope, { key: "ArrowDown" });
    expect(arrow).toHaveBeenCalledTimes(1);
  });
});
