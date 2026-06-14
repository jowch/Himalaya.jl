import { describe, it, expect } from "vitest";
import { peakGlyph } from "../src/print/ui/peakMark";

const PHASE = "oklch(0.570 0.150 58)"; // a resolved phase colour (Pn3m amber)

describe("peakGlyph — §5.1 encoding atoms", () => {
  it("auto → filled downward triangle", () => {
    const g = peakGlyph({ source: "auto", color: PHASE });
    expect(g.shape).toBe("triangle-down");
    expect(g.fill).toBe(PHASE);
    expect(g.stroke).toBe(PHASE);
    expect(g.interactive).toBe(true);
  });

  it("manual (indexed) → filled diamond in the phase colour", () => {
    const g = peakGlyph({ source: "manual", color: PHASE });
    expect(g.shape).toBe("diamond");
    expect(g.fill).toBe(PHASE);
  });

  it("manual unindexed → neutral-gray (caller passes the neutral colour)", () => {
    const NEUTRAL = "var(--color-ink-faint)";
    const g = peakGlyph({ source: "manual", color: NEUTRAL });
    expect(g.shape).toBe("diamond");
    expect(g.stroke).toBe(NEUTRAL);
    expect(g.fill).toBe(NEUTRAL);
  });

  it("excluded → ghosted (hollow, not struck)", () => {
    const g = peakGlyph({ source: "auto", color: PHASE, excluded: true });
    expect(g.fill).toBe("none");
    expect(g.stroke).toBe(PHASE);
  });

  it("optimistic (id<0) → outline-only + non-interactive", () => {
    const g = peakGlyph({ source: "manual", color: PHASE, optimistic: true });
    expect(g.fill).toBe("none");
    expect(g.interactive).toBe(false);
    // outline-only draws a slightly heavier stroke so the empty glyph reads.
    expect(g.strokeWidth).toBeGreaterThan(
      peakGlyph({ source: "manual", color: PHASE }).strokeWidth,
    );
  });

  it("predicted-but-absent → hollow caret", () => {
    const g = peakGlyph({ source: "auto", color: PHASE, predictedAbsent: true });
    expect(g.shape).toBe("caret");
    expect(g.fill).toBe("none");
  });

  it("hot (q-link) → sets the ring flag but does NOT grow or recolour the mark", () => {
    const base = peakGlyph({ source: "auto", color: PHASE, r: 4 });
    const hot = peakGlyph({ source: "auto", color: PHASE, r: 4, hot: true });
    // Greenfield q-link redesign: the mark is UNCHANGED when hot — emphasis is
    // carried by the q-line + q-readout, not by growing the glyph. `ring` is now
    // a pure flag (surfaced as data-hot on the SVG), no geometry change.
    expect(hot.r).toBe(base.r);
    expect(hot.ring).toBe(true);
    // hue is unchanged — provenance/state never rides on colour.
    expect(hot.fill).toBe(base.fill);
    expect(base.ring).toBe(false);
  });

  it("offsetPx defaults to 7 and is parameterizable (member layers use 5)", () => {
    expect(peakGlyph({ source: "auto", color: PHASE }).offsetPx).toBe(7);
    expect(peakGlyph({ source: "auto", color: PHASE, offsetPx: 5 }).offsetPx).toBe(5);
  });

  // No test asserts magenta — the manual atom is a diamond, full stop.
  it("manual is NOT magenta — geometry (diamond) carries the provenance", () => {
    const g = peakGlyph({ source: "manual", color: PHASE });
    expect(g.shape).toBe("diamond");
    expect(g.fill).not.toMatch(/340/); // the retired magenta hue
  });
});
