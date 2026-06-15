import { describe, it, expect } from "vitest";
import {
  SHORTCUTS,
  matchShortcut,
  keyComboLabel,
  shortcutLabel,
} from "../../src/print/shell/shortcuts";

/** Build a minimal KeyboardEvent-like object for the matcher. */
function ev(
  key: string,
  mods: { meta?: boolean; ctrl?: boolean; alt?: boolean; shift?: boolean } = {},
): KeyboardEvent {
  return {
    key,
    metaKey: !!mods.meta,
    ctrlKey: !!mods.ctrl,
    altKey: !!mods.alt,
    shiftKey: !!mods.shift,
  } as KeyboardEvent;
}

describe("shortcut registry", () => {
  it("every entry has a unique, non-empty key binding", () => {
    const seen = new Set<string>();
    for (const def of Object.values(SHORTCUTS)) {
      expect(def.keys.length).toBeGreaterThan(0);
      for (const k of def.keys) {
        expect(seen.has(k)).toBe(false); // no two actions share a binding
        seen.add(k);
      }
    }
  });

  it("the sample stepper is [ and ] only (', '.' retired)", () => {
    expect(SHORTCUTS.prevSample.keys).toEqual(["["]);
    expect(SHORTCUTS.nextSample.keys).toEqual(["]"]);
    // , and . must not appear anywhere in the registry
    const all = Object.values(SHORTCUTS).flatMap((d) => d.keys);
    expect(all).not.toContain(",");
    expect(all).not.toContain(".");
  });
});

describe("matchShortcut", () => {
  it("matches bare letter/bracket keys only when NO modifier is held", () => {
    expect(matchShortcut(ev("x"))).toBe("drop");
    expect(matchShortcut(ev("k"))).toBe("keep");
    expect(matchShortcut(ev("r"))).toBe("representative");
    expect(matchShortcut(ev("["))).toBe("prevSample");
    expect(matchShortcut(ev("]"))).toBe("nextSample");
    // a held Mod means it is NOT the page shortcut (e.g. Cmd+] = browser nav)
    expect(matchShortcut(ev("]", { meta: true }))).toBeNull();
    expect(matchShortcut(ev("x", { ctrl: true }))).toBeNull();
  });

  it("lowercases letters so CapsLock X still drops, but Shift+X does not", () => {
    expect(matchShortcut(ev("X"))).toBe("drop"); // capslock: shiftKey false
    expect(matchShortcut(ev("X", { shift: true }))).toBeNull(); // shift held = different combo
  });

  it("matches arrows for the Focus/Loupe sub-entity steps", () => {
    expect(matchShortcut(ev("ArrowLeft"))).toBe("prevExposure");
    expect(matchShortcut(ev("ArrowRight"))).toBe("nextExposure");
    expect(matchShortcut(ev("ArrowUp"))).toBe("prevCandidate");
    expect(matchShortcut(ev("ArrowDown"))).toBe("nextCandidate");
  });

  it("matches Mod+Z undo and Mod+Shift+Z redo cross-platform (meta OR ctrl)", () => {
    expect(matchShortcut(ev("z", { meta: true }))).toBe("undo");
    expect(matchShortcut(ev("z", { ctrl: true }))).toBe("undo");
    expect(matchShortcut(ev("z", { meta: true, shift: true }))).toBe("redo");
    // plain z is not undo
    expect(matchShortcut(ev("z"))).toBeNull();
  });

  it("matches Alt+Arrow reorder and the find combos", () => {
    expect(matchShortcut(ev("ArrowUp", { alt: true }))).toBe("reorderUp");
    expect(matchShortcut(ev("ArrowDown", { alt: true }))).toBe("reorderDown");
    expect(matchShortcut(ev("k", { meta: true }))).toBe("find");
    expect(matchShortcut(ev("/"))).toBe("find");
  });

  it("matches Escape to dismiss", () => {
    expect(matchShortcut(ev("Escape"))).toBe("dismiss");
  });
});

describe("keyComboLabel (display)", () => {
  it("renders mac glyphs", () => {
    expect(keyComboLabel("Mod+z", true)).toBe("⌘Z");
    expect(keyComboLabel("Mod+Shift+z", true)).toBe("⌘⇧Z");
    expect(keyComboLabel("Alt+ArrowUp", true)).toBe("⌥↑");
    expect(keyComboLabel("ArrowLeft", true)).toBe("←");
    expect(keyComboLabel("Escape", true)).toBe("Esc");
    expect(keyComboLabel("[", true)).toBe("[");
  });

  it("renders non-mac words", () => {
    expect(keyComboLabel("Mod+z", false)).toBe("Ctrl+Z");
    expect(keyComboLabel("Alt+ArrowUp", false)).toBe("Alt+↑");
  });

  it("shortcutLabel shows the first binding of an action", () => {
    expect(shortcutLabel("prevSample", true)).toBe("[");
    expect(shortcutLabel("undo", true)).toBe("⌘Z");
  });
});
