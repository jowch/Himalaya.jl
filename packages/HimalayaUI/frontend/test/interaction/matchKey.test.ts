import { describe, it, expect } from "vitest";
import { comboOf, matchesKeys, isTyping } from "../../src/print/interaction/matchKey";

function key(init: Partial<KeyboardEvent> & { key: string }): KeyboardEvent {
  return new KeyboardEvent("keydown", init);
}

describe("comboOf", () => {
  it("normalizes a bare letter to lowercase", () => {
    expect(comboOf(key({ key: "X" }))).toBe("x");
  });
  it("emits Mod for meta or ctrl", () => {
    expect(comboOf(key({ key: "z", metaKey: true }))).toBe("Mod+z");
    expect(comboOf(key({ key: "z", ctrlKey: true }))).toBe("Mod+z");
  });
  it("orders modifiers Mod+Alt+Shift then key", () => {
    expect(comboOf(key({ key: "ArrowUp", shiftKey: true }))).toBe("Shift+ArrowUp");
    expect(comboOf(key({ key: "z", metaKey: true, shiftKey: true }))).toBe("Mod+Shift+z");
  });
  it("maps space to Space", () => {
    expect(comboOf(key({ key: " " }))).toBe("Space");
  });
  it("passes ? through as bare '?' even though the real keystroke carries Shift", () => {
    // A real `?` is Shift+/, so the event has shiftKey:true and key:"?". The
    // shifted-punctuation char already encodes Shift → combo is bare "?", not
    // "Shift+?" (which would never match the help action's keys:["?"]).
    expect(comboOf(key({ key: "?", shiftKey: true }))).toBe("?");
  });
  it("keeps Shift as a modifier for alphanumerics (Mod+Shift+z redo)", () => {
    expect(comboOf(key({ key: "z", metaKey: true, shiftKey: true }))).toBe("Mod+Shift+z");
  });
});

describe("matchesKeys", () => {
  it("matches any combo in the list, case-insensitively for letters", () => {
    expect(matchesKeys(key({ key: "K" }), ["k"])).toBe(true);
    expect(matchesKeys(key({ key: "z", metaKey: true }), ["Mod+z"])).toBe(true);
    expect(matchesKeys(key({ key: "z" }), ["Mod+z"])).toBe(false);
  });
  it("matches a declared \" \" against the space bar (Space-normalized both sides)", () => {
    expect(matchesKeys(key({ key: " " }), [" "])).toBe(true);
    expect(matchesKeys(key({ key: " " }), ["Space"])).toBe(true);
  });
});

describe("isTyping", () => {
  it("detects inputs, textareas, and contenteditable", () => {
    const input = document.createElement("input");
    expect(isTyping(input)).toBe(true);
    const div = document.createElement("div");
    expect(isTyping(div)).toBe(false);
    div.setAttribute("contenteditable", "true");
    expect(isTyping(div)).toBe(true);
  });
});
