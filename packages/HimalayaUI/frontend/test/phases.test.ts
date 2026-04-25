import { describe, it, expect } from "vitest";
import { phaseColor, KNOWN_PHASES } from "../src/phases";

// Accepts CSS hex (#aabbcc), oklch(...), or any other valid CSS color string —
// we just want to confirm phaseColor returns a non-empty string.
const COLOR_RE = /^(#[0-9a-f]{3,8}|oklch\(.+\)|rgb\(.+\)|var\(--.+\))$/i;

describe("phases", () => {
  it("returns a distinct color for each known phase", () => {
    const colors = new Set(KNOWN_PHASES.map((p) => phaseColor(p)));
    expect(colors.size).toBe(KNOWN_PHASES.length);
    for (const c of colors) expect(c).toMatch(COLOR_RE);
  });

  it("returns a fallback color for unknown phases", () => {
    expect(phaseColor("Unknown")).toMatch(COLOR_RE);
  });
});
