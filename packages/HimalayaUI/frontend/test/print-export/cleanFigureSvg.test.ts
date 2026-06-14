import { describe, it, expect } from "vitest";
import { buildCleanFigureSvg } from "../../src/print/export/cleanFigureSvg";
import type { WaterfallRow } from "../../src/print/waterfall/waterfallModel";

function trace(): { q: number[]; I: number[]; sigma: number[] } {
  const q = Array.from({ length: 40 }, (_, i) => 0.01 + i * 0.005);
  const I = q.map((qq, i) => 1000 / (1 + qq * 50) + i); // monotonic-ish, all > 0
  return { q, I, sigma: q.map(() => 1) };
}

function row(over: Partial<WaterfallRow> = {}): WaterfallRow {
  return {
    key: "1",
    label: "exp 1",
    trace: trace(),
    phase: "Im3m",
    state: "indexed",
    anchors: [
      { id: 10, q: 0.04, intensity: 200, phase: "Im3m" },
      { id: 11, q: 0.06, intensity: 150, phase: "Im3m" },
    ],
    bandHeight: 1,
    yOffset: 0,
    ...over,
  };
}

const BASE = {
  title: "LL37 titration",
  footer: "q normalized · intensity offset for clarity",
};

describe("buildCleanFigureSvg", () => {
  it("returns standalone, well-formed SVG markup with the title + axis titles", () => {
    const svg = buildCleanFigureSvg({ rows: [row()], ...BASE });
    expect(svg.startsWith("<svg")).toBe(true);
    expect(svg.trimEnd().endsWith("</svg>")).toBe(true);
    // Parses round-trip (no XML error).
    const doc = new DOMParser().parseFromString(svg, "image/svg+xml");
    expect(doc.querySelector("parsererror")).toBeNull();
    expect(svg).toContain("LL37 titration");
    expect(svg).toContain("q (Å⁻¹)");
    expect(svg).toContain("Intensity (a.u.) + offset");
    expect(svg).toContain(BASE.footer);
  });

  it("uses the clean export skin: white ground, Arial, literal hex (no Print OKLCH / var())", () => {
    const svg = buildCleanFigureSvg({ rows: [row()], ...BASE });
    expect(svg).toContain('fill="#ffffff"'); // white background rect
    expect(svg).toContain("Arial");
    expect(svg).not.toMatch(/var\(--/);
    expect(svg).not.toMatch(/oklch\(/);
  });

  it("colours each trace by its phase hex (Im3m → sage #007d53, unphased → grey)", () => {
    const svg = buildCleanFigureSvg({
      rows: [row({ key: "a", phase: "Im3m" }), row({ key: "b", phase: null, label: "exp 2", anchors: [] })],
      ...BASE,
    });
    expect(svg).toContain('stroke="#007d53"'); // Im3m sage
    expect(svg).toContain('stroke="#777777"'); // unphased fallback
  });

  it("draws peak glyphs only when showPeakTicks, and ordinal labels only when showPeakLabels", () => {
    const ticksOnly = buildCleanFigureSvg({ rows: [row()], showPeakTicks: true, showPeakLabels: false, ...BASE });
    // Two anchors → two glyph triangles (filled paths in the phase colour).
    expect((ticksOnly.match(/d="M[^"]*Z" fill="#007d53"/g) ?? []).length).toBe(2);
    // No ordinal label text yet.
    expect(ticksOnly).not.toMatch(/text-anchor="middle" font-size="9" fill="#007d53">1</);

    const withLabels = buildCleanFigureSvg({ rows: [row()], showPeakTicks: true, showPeakLabels: true, ...BASE });
    expect(withLabels).toContain(">1</text>");
    expect(withLabels).toContain(">2</text>");

    const noPeaks = buildCleanFigureSvg({ rows: [row()], showPeakTicks: false, ...BASE });
    expect(noPeaks).not.toMatch(/d="M[^"]*Z" fill="#007d53"/);
  });

  it("renders a phase legend (one entry per distinct phase, unphased spelled out)", () => {
    const svg = buildCleanFigureSvg({
      rows: [row({ phase: "Im3m" }), row({ key: "b", phase: null, anchors: [] })],
      ...BASE,
    });
    expect(svg).toContain(">Im3m</text>");
    expect(svg).toContain(">unphased / unbound</text>");
  });

  it("honours the q-domain and renders a row label per trace", () => {
    const svg = buildCleanFigureSvg({
      rows: [row({ label: "exp 20" }), row({ key: "b", label: "exp 6" })],
      ...BASE,
    });
    expect(svg).toContain(">exp 20</text>");
    expect(svg).toContain(">exp 6</text>");
  });
});
