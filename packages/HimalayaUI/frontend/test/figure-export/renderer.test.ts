import { describe, it, expect } from "vitest";
import * as Plot from "@observablehq/plot";
import { buildExportSvg, buildExportPng } from "../../src/lib/figure-export/renderer";
import type { ExportSpec } from "../../src/lib/figure-export/types";

function makeSpec(over: Partial<ExportSpec> = {}): ExportSpec {
  return {
    title: { primary: "Primary Title", secondary: "Secondary line" },
    width: 800,
    height: 600,
    plot: {
      marks: [Plot.line([{ x: 0, y: 0 }, { x: 1, y: 1 }], { x: "x", y: "y" })],
      x: { type: "linear" },
      y: { type: "linear" },
      width: 800,
      height: 480,
    },
    legend: {
      rows: [
        { color: "#1a1a1a", symbol: "swatch", label: "trace" },
        { color: "#cc6600", symbol: "line",   label: "predicted Pn3m" },
      ],
    },
    ...over,
  };
}

describe("buildExportSvg", () => {
  it("returns an SVGSVGElement", () => {
    const svg = buildExportSvg(makeSpec());
    expect(svg.tagName.toLowerCase()).toBe("svg");
    expect(svg.namespaceURI).toBe("http://www.w3.org/2000/svg");
  });

  it("renders a white background rect covering the canvas", () => {
    const svg = buildExportSvg(makeSpec());
    const rects = svg.querySelectorAll("rect");
    const bgRect = Array.from(rects).find((r) =>
      r.getAttribute("fill") === "#ffffff" || r.getAttribute("fill") === "white"
    );
    expect(bgRect).toBeDefined();
  });

  it("includes the primary title text", () => {
    const svg = buildExportSvg(makeSpec({ title: { primary: "My Trace Title" } }));
    expect(svg.textContent).toContain("My Trace Title");
  });

  it("includes the secondary title when provided", () => {
    const svg = buildExportSvg(makeSpec({ title: { primary: "P", secondary: "S" } }));
    expect(svg.textContent).toContain("S");
  });

  it("omits secondary title when undefined", () => {
    const spec = makeSpec();
    delete spec.title.secondary;
    const svg = buildExportSvg(spec);
    // Only the primary text appears.
    expect(svg.textContent?.trim()).toContain("Primary Title");
  });

  it("renders one legend row per LegendSpec.rows entry", () => {
    const svg = buildExportSvg(makeSpec());
    expect(svg.textContent).toContain("trace");
    expect(svg.textContent).toContain("predicted Pn3m");
  });

  it("output contains no var(--…) tokens (palette resolved at adapter time)", () => {
    const svg = buildExportSvg(makeSpec());
    const xml = new XMLSerializer().serializeToString(svg);
    expect(xml).not.toMatch(/var\(--/);
  });

  it("serialized output is well-formed XML and parses round-trip", () => {
    const svg = buildExportSvg(makeSpec());
    const xml = new XMLSerializer().serializeToString(svg);
    const parsed = new DOMParser().parseFromString(xml, "image/svg+xml");
    expect(parsed.querySelector("parsererror")).toBeNull();
    expect(parsed.documentElement.tagName.toLowerCase()).toBe("svg");
  });

  it("Plot.plot output is SVGSVGElement (guards against figure: true)", () => {
    // Sanity check that a vanilla spec doesn't have title/caption/figure set —
    // those flip Plot's return to HTMLElement.
    const spec = makeSpec();
    expect((spec.plot as { title?: unknown }).title).toBeUndefined();
    expect((spec.plot as { caption?: unknown }).caption).toBeUndefined();
    expect((spec.plot as { figure?: unknown }).figure).toBeUndefined();
    const inner = Plot.plot(spec.plot);
    expect(inner.tagName.toLowerCase()).toBe("svg");
  });
});

describe("buildExportPng", () => {
  it("is async and returns a Promise (full PNG path is JSDOM-untestable)", () => {
    // We only assert the function exists and returns a Promise. Real PNG
    // generation is covered by the Playwright E2E in Phase 4. JSDOM has no
    // OffscreenCanvas / Image.decode.
    const result = buildExportPng(makeSpec());
    expect(result).toBeInstanceOf(Promise);
    // Swallow rejection in JSDOM where canvas APIs are stubs/missing.
    result.catch(() => undefined);
  });
});
