import { describe, it, expect } from "vitest";
import {
  SHORTCUTS,
  eventCombo,
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

  it("the sample axis is ↑/↓ and the detail axis is ←/→ (rev-2 axes)", () => {
    expect(SHORTCUTS.prevSample.keys).toEqual(["ArrowUp"]);
    expect(SHORTCUTS.nextSample.keys).toEqual(["ArrowDown"]);
    expect(SHORTCUTS.prevDetail.keys).toEqual(["ArrowLeft"]);
    expect(SHORTCUTS.nextDetail.keys).toEqual(["ArrowRight"]);
    // old [ ] bindings must not appear anywhere in the registry
    const all = Object.values(SHORTCUTS).flatMap((d) => d.keys);
    expect(all).not.toContain("[");
    expect(all).not.toContain("]");
    expect(all).not.toContain(",");
    expect(all).not.toContain(".");
  });
});

describe("matchShortcut", () => {
  it("matches bare letter keys only when NO modifier is held", () => {
    expect(matchShortcut(ev("x"))).toBe("drop");
    expect(matchShortcut(ev("k"))).toBe("keep");
    expect(matchShortcut(ev("r"))).toBe("representative");
    // old [ ] are unbound in rev-2; arrows now drive sample/detail nav
    expect(matchShortcut(ev("["))).toBeNull();
    expect(matchShortcut(ev("]"))).toBeNull();
    // a held Mod means it is NOT the page shortcut
    expect(matchShortcut(ev("x", { ctrl: true }))).toBeNull();
  });

  it("lowercases letters so CapsLock X still drops, but Shift+X does not", () => {
    expect(matchShortcut(ev("X"))).toBe("drop"); // capslock: shiftKey false
    expect(matchShortcut(ev("X", { shift: true }))).toBeNull(); // shift held = different combo
  });

  it("matches arrows for the rev-2 sample/detail axes", () => {
    expect(matchShortcut(ev("ArrowUp"))).toBe("prevSample");
    expect(matchShortcut(ev("ArrowDown"))).toBe("nextSample");
    expect(matchShortcut(ev("ArrowLeft"))).toBe("prevDetail");
    expect(matchShortcut(ev("ArrowRight"))).toBe("nextDetail");
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

  it("new verbs are bound: Enter→openFocus, Backspace→restore, Space→toggleSelect", () => {
    expect(matchShortcut(ev("Enter"))).toBe("openFocus");
    expect(matchShortcut(ev("Backspace"))).toBe("restore");
    expect(matchShortcut(ev(" "))).toBe("toggleSelect");
  });
});

describe("eventCombo (normalization)", () => {
  it("? normalizes to a stable token regardless of Shift being held", () => {
    expect(eventCombo(ev("?", { shift: true }))).toBe("?");
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
    expect(shortcutLabel("prevSample", true)).toBe("↑");
    expect(shortcutLabel("undo", true)).toBe("⌘Z");
  });
});
