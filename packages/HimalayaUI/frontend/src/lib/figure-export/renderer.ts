// renderer.ts — turn an ExportSpec into SVG (string-serializable) or PNG.
import * as Plot from "@observablehq/plot";
import type { ExportSpec, LegendRow } from "./types";
import {
  EXPORT_MARGIN,
  LIGHT_PALETTE,
} from "./presets";

const SVG_NS = "http://www.w3.org/2000/svg";

/** Fallback sans stack when a spec doesn't pin a family (clean preset pins Arial). */
const EXPORT_SANS = "Arial, Helvetica, system-ui, sans-serif";

/**
 * Set the three SVG font presentation attributes explicitly. The CSS `font`
 * shorthand is NOT a valid SVG presentation attribute (it is silently ignored),
 * so text authored with `setAttribute("font", "600 16px …")` falls back to the
 * UA default size AND, with no family, the default serif — the export's old
 * "sloppy" look. Always set family/size/weight as their own attributes.
 */
function applyFont(el: SVGElement, family: string, sizePx: number, weight: number): void {
  el.setAttribute("font-family", family);
  el.setAttribute("font-size", String(sizePx));
  el.setAttribute("font-weight", String(weight));
}

/**
 * Layout the export SVG: white background, primary + secondary title,
 * Plot.plot() body, legend rows beneath. Returns a fresh SVGSVGElement
 * not attached to the document.
 */
export function buildExportSvg(spec: ExportSpec): SVGSVGElement {
  const svg = document.createElementNS(SVG_NS, "svg");
  // Note: `createElementNS` already declares the SVG namespace on serialization;
  // setting `xmlns` explicitly produces a duplicate-attribute parse error in
  // strict XML serializers (e.g. JSDOM). Browsers ship the namespace declaration
  // automatically when the element is in the SVG namespace.
  svg.setAttribute("width", String(spec.width));
  svg.setAttribute("height", String(spec.height));
  svg.setAttribute("viewBox", `0 0 ${spec.width} ${spec.height}`);

  // White background.
  const bg = document.createElementNS(SVG_NS, "rect");
  bg.setAttribute("x", "0");
  bg.setAttribute("y", "0");
  bg.setAttribute("width", String(spec.width));
  bg.setAttribute("height", String(spec.height));
  bg.setAttribute("fill", LIGHT_PALETTE.background);
  svg.appendChild(bg);

  // One resolved family for every hand-drawn text element (clean preset → Arial).
  const family = spec.fontFamily ?? EXPORT_SANS;

  // Title block (primary + secondary).
  const titleY = 24;
  const primary = document.createElementNS(SVG_NS, "text");
  primary.setAttribute("x", String(EXPORT_MARGIN.left));
  primary.setAttribute("y", String(titleY));
  applyFont(primary, family, 16, 600);
  primary.setAttribute("fill", LIGHT_PALETTE.text);
  primary.textContent = spec.title.primary;
  svg.appendChild(primary);

  if (spec.title.secondary) {
    const secondary = document.createElementNS(SVG_NS, "text");
    secondary.setAttribute("x", String(EXPORT_MARGIN.left));
    secondary.setAttribute("y", String(titleY + 18));
    applyFont(secondary, family, 12, 400);
    secondary.setAttribute("fill", LIGHT_PALETTE.textMuted);
    secondary.textContent = spec.title.secondary;
    svg.appendChild(secondary);
  }

  // Plot body — Plot.plot() returns SVGSVGElement when title/caption/figure
  // are NOT set (adapters must NOT set them).
  const plotEl = Plot.plot({
    style: {
      background: "transparent",
      color: LIGHT_PALETTE.text,
      // A real family, not a CSS `font` shorthand — Plot writes this straight to
      // the body's font-family (the clean preset pins Arial).
      fontFamily: spec.fontFamily ?? EXPORT_SANS,
    },
    ...spec.plot,
  });
  // Position the plot below the title block. NEST the Plot <svg> directly — do
  // NOT unwrap its children into a <g>. Observable Plot scopes its generated
  // <style> to a `:where(.plot-XXXX)` class on its ROOT <svg>; discarding that
  // root orphans every scoped rule (font-family included) and the body falls
  // back to the UA default serif. A nested <svg x y> keeps the class, the
  // <style> block, and the inline styles self-contained and correctly placed.
  plotEl.setAttribute("x", String(EXPORT_MARGIN.left));
  plotEl.setAttribute("y", String(EXPORT_MARGIN.top));
  // A nested <svg> clips to its viewport by default; Plot draws the x-axis tick
  // labels a few px BELOW its height box, so they would be cut. `overflow:visible`
  // renders them (like the old unwrapping <g>) while keeping the scoped styles.
  plotEl.setAttribute("overflow", "visible");
  svg.appendChild(plotEl);

  // Legend rows (beneath the plot, inside the bottom margin).
  if (spec.legend && spec.legend.rows.length > 0) {
    const legendY = spec.height - 28;
    const legendG = document.createElementNS(SVG_NS, "g");
    legendG.setAttribute("transform", `translate(${EXPORT_MARGIN.left}, ${legendY})`);
    let cursorX = 0;
    for (const row of spec.legend.rows) {
      const item = renderLegendItem(row, cursorX, family);
      legendG.appendChild(item.group);
      cursorX += item.width + 16;
    }
    svg.appendChild(legendG);
  }

  // Centered footnote line under the plot (clean preset, E-8).
  if (spec.footnote) {
    const foot = document.createElementNS(SVG_NS, "text");
    foot.setAttribute("x", String(spec.width / 2));
    foot.setAttribute("y", String(spec.height - 10));
    foot.setAttribute("text-anchor", "middle");
    applyFont(foot, family, 12, 400);
    foot.setAttribute("fill", LIGHT_PALETTE.textMuted);
    foot.textContent = spec.footnote;
    svg.appendChild(foot);
  }

  return svg;
}

function renderLegendItem(row: LegendRow, x: number, family: string): { group: SVGGElement; width: number } {
  const g = document.createElementNS(SVG_NS, "g");
  g.setAttribute("transform", `translate(${x}, 0)`);

  if (row.symbol === "swatch") {
    const sw = document.createElementNS(SVG_NS, "rect");
    sw.setAttribute("x", "0");
    sw.setAttribute("y", "-9");
    sw.setAttribute("width", "10");
    sw.setAttribute("height", "10");
    sw.setAttribute("fill", row.color);
    g.appendChild(sw);
  } else {
    const ln = document.createElementNS(SVG_NS, "line");
    ln.setAttribute("x1", "0");
    ln.setAttribute("y1", "-4");
    ln.setAttribute("x2", "0");
    ln.setAttribute("y2", "-12");
    ln.setAttribute("stroke", row.color);
    ln.setAttribute("stroke-width", "2");
    g.appendChild(ln);
  }

  const text = document.createElementNS(SVG_NS, "text");
  text.setAttribute("x", "16");
  text.setAttribute("y", "0");
  applyFont(text, family, 12, 400);
  text.setAttribute("fill", LIGHT_PALETTE.text);
  text.textContent = row.label;
  g.appendChild(text);

  // Approximate width — symbol + gap + len(label)*7. Good enough for layout.
  const width = 16 + row.label.length * 7;
  return { group: g, width };
}

/**
 * Pre-flight feature check for the PNG path. Mirrors `canCopyPngToClipboard()`
 * in clipboard.ts — the Save button uses this to disable itself on browsers
 * that lack the OffscreenCanvas/convertToBlob pipeline. Without this gate the
 * user clicks Save and only finds out the export failed when the error toast
 * fires; with it the Save button stays disabled with no surprise.
 */
export function canExportPng(): boolean {
  if (typeof OffscreenCanvas === "undefined") return false;
  // convertToBlob landed alongside OffscreenCanvas in evergreen browsers but
  // is missing on some Safari versions where OffscreenCanvas exists.
  const proto = (OffscreenCanvas as unknown as { prototype?: { convertToBlob?: unknown } }).prototype;
  if (!proto || typeof proto.convertToBlob !== "function") return false;
  if (typeof Image === "undefined") return false;
  return true;
}

/**
 * Render the export as a 2× DPI PNG blob. Pipeline: SVG → blob URL → Image
 * decode → OffscreenCanvas drawImage → convertToBlob.
 *
 * `URL.revokeObjectURL` runs in a `finally` so a thrown render still cleans
 * up (otherwise every export leaks an object URL). JSDOM doesn't have these
 * canvas APIs — full path covered by the Playwright E2E in Phase 4.
 */
export async function buildExportPng(spec: ExportSpec, scale = 2): Promise<Blob> {
  const svg = buildExportSvg(spec);
  const xml = new XMLSerializer().serializeToString(svg);
  const url = URL.createObjectURL(new Blob([xml], { type: "image/svg+xml" }));
  try {
    const img = new Image();
    img.src = url;
    await img.decode();
    const off = new OffscreenCanvas(spec.width * scale, spec.height * scale);
    const ctx = off.getContext("2d");
    if (!ctx) throw new Error("OffscreenCanvas 2d context unavailable");
    ctx.drawImage(img, 0, 0, off.width, off.height);
    return await off.convertToBlob({ type: "image/png" });
  } finally {
    URL.revokeObjectURL(url);
  }
}
