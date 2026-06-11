import { describe, it, expect, afterEach } from "vitest";
import { fireEvent } from "@testing-library/dom";
import { suppressGlobalKeys } from "../../src/lib/keys";

// Dispatch a real keydown from `target`, capture it at the window (where the
// page-level shortcut listeners live), and return what the predicate said.
function probe(target: Element | Window): boolean {
  let result: boolean | null = null;
  const listener = (e: Event): void => {
    result = suppressGlobalKeys(e as KeyboardEvent);
  };
  window.addEventListener("keydown", listener);
  fireEvent.keyDown(target, { key: "x" });
  window.removeEventListener("keydown", listener);
  if (result === null) throw new Error("keydown never reached the window listener");
  return result;
}

function mount<K extends keyof HTMLElementTagNameMap>(tag: K): HTMLElementTagNameMap[K] {
  const el = document.createElement(tag);
  document.body.appendChild(el);
  return el;
}

afterEach(() => {
  document.body.innerHTML = "";
});

describe("suppressGlobalKeys — typing/selection targets", () => {
  const editingTags = ["input", "textarea", "select"] as const;
  for (const tag of editingTags) {
    it(`suppresses when the target is a <${tag}>`, () => {
      expect(probe(mount(tag))).toBe(true);
    });
  }

  it("suppresses when the target is a contenteditable element", () => {
    const div = mount("div");
    div.setAttribute("contenteditable", "true");
    expect(probe(div)).toBe(true);
  });

  it("suppresses when the target sits INSIDE a contenteditable region", () => {
    const div = mount("div");
    div.setAttribute("contenteditable", "");
    const span = document.createElement("span");
    div.appendChild(span);
    expect(probe(span)).toBe(true);
  });

  it("does not suppress a contenteditable=\"false\" island", () => {
    const div = mount("div");
    div.setAttribute("contenteditable", "false");
    expect(probe(div)).toBe(false);
  });

  it("does not suppress a plain button target", () => {
    expect(probe(mount("button"))).toBe(false);
  });

  it("does not suppress the body / window target", () => {
    expect(probe(document.body)).toBe(false);
    expect(probe(window)).toBe(false);
  });
});

describe("suppressGlobalKeys — open modal dialog", () => {
  it("suppresses regardless of the event target while an aria-modal dialog is open", () => {
    const dialog = mount("div");
    dialog.setAttribute("role", "dialog");
    dialog.setAttribute("aria-modal", "true");
    // Target is a plain button OUTSIDE the dialog — aria-modal means everything
    // outside is inert, so the page shortcut must still be off.
    expect(probe(mount("button"))).toBe(true);
    // …and from the window itself (focus can transiently sit on body before
    // the focus trap engages).
    expect(probe(window)).toBe(true);
  });

  it("does NOT suppress for a non-modal role=dialog (no aria-modal)", () => {
    const dialog = mount("div");
    dialog.setAttribute("role", "dialog");
    expect(probe(mount("button"))).toBe(false);
  });
});

describe("suppressGlobalKeys — already-consumed event (defaultPrevented)", () => {
  // The real bug this guards: ModalShell's `document` Escape listener
  // preventDefault()s and closes the modal; React can unmount it before the
  // page's `window` listener runs, so the open-dialog check sees nothing.
  // `defaultPrevented` survives that race. Here a closer listener consumes the
  // key (no dialog present at all), and the window-level predicate must still
  // suppress.
  it("suppresses when an upstream listener has called preventDefault — even with no dialog", () => {
    const target = mount("button");
    let result: boolean | null = null;
    const consumer = (e: Event): void => e.preventDefault();
    const probeListener = (e: Event): void => {
      result = suppressGlobalKeys(e as KeyboardEvent);
    };
    // document consumes first (bubble: document before window), window observes.
    document.addEventListener("keydown", consumer);
    window.addEventListener("keydown", probeListener);
    fireEvent.keyDown(target, { key: "Escape" });
    window.removeEventListener("keydown", probeListener);
    document.removeEventListener("keydown", consumer);
    expect(result).toBe(true);
  });
});
