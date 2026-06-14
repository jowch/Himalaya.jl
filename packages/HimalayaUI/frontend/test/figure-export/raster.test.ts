import { describe, it, expect } from "vitest";
import { svgScaledToPixels, EXPORT_PNG_DPI } from "../../src/lib/figure-export/raster";

describe("svgScaledToPixels", () => {
  it("rewrites root width/height to target px but preserves the point viewBox", () => {
    const svg = '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="300" viewBox="0 0 200 300"><rect/></svg>';
    const out = svgScaledToPixels(svg, 600, 900);
    expect(out).toContain('width="600"');
    expect(out).toContain('height="900"');
    expect(out).toContain('viewBox="0 0 200 300"'); // coordinate system unchanged
    expect(out).toContain("<rect/>"); // body untouched
  });

  it("only touches the root <svg> dims, not inner element width/height", () => {
    const svg = '<svg width="100" height="120" viewBox="0 0 100 120"><rect x="0" y="0" width="11" height="11"/></svg>';
    const out = svgScaledToPixels(svg, 300, 360);
    expect(out).toContain('<rect x="0" y="0" width="11" height="11"/>'); // swatch untouched
    expect(out).toMatch(/<svg[^>]*width="300"[^>]*height="360"/);
  });

  it("synthesizes a viewBox when the source has none", () => {
    const out = svgScaledToPixels('<svg width="50" height="40"></svg>', 150, 120);
    expect(out).toContain('viewBox="0 0 50 40"');
  });

  it("the export DPI is print-grade (300, ≥4× the 72pt authoring grid)", () => {
    expect(EXPORT_PNG_DPI).toBe(300);
    expect(EXPORT_PNG_DPI / 72).toBeGreaterThanOrEqual(3);
  });
});
