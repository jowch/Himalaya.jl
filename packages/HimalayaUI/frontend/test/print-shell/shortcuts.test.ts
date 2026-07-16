import { describe, it, expect } from "vitest";
import {
  SHORTCUTS,
  keyComboLabel,
  shortcutLabel,
} from "../../src/print/shell/shortcuts";

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
