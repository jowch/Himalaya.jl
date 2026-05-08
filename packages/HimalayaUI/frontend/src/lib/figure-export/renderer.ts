// renderer.ts — turn an ExportSpec into SVG (string-serializable) or PNG.
import * as Plot from "@observablehq/plot";
import type { ExportSpec, LegendRow } from "./types";
import {
  EXPORT_MARGIN,
  TITLE_FONT_PRIMARY,
  TITLE_FONT_SECONDARY,
  BODY_FONT,
  LIGHT_PALETTE,
} from "./presets";

const SVG_NS = "http://www.w3.org/2000/svg";

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

  // Title block (primary + secondary).
  const titleY = 24;
  const primary = document.createElementNS(SVG_NS, "text");
  primary.setAttribute("x", String(EXPORT_MARGIN.left));
  primary.setAttribute("y", String(titleY));
  primary.setAttribute("font", TITLE_FONT_PRIMARY);
  primary.setAttribute("fill", LIGHT_PALETTE.text);
  primary.textContent = spec.title.primary;
  svg.appendChild(primary);

  if (spec.title.secondary) {
    const secondary = document.createElementNS(SVG_NS, "text");
    secondary.setAttribute("x", String(EXPORT_MARGIN.left));
    secondary.setAttribute("y", String(titleY + 18));
    secondary.setAttribute("font", TITLE_FONT_SECONDARY);
    secondary.setAttribute("fill", LIGHT_PALETTE.textMuted);
    secondary.textContent = spec.title.secondary;
    svg.appendChild(secondary);
  }

  // Plot body — Plot.plot() returns SVGSVGElement when title/caption/figure
  // are NOT set (adapters must NOT set them).
  const plotEl = Plot.plot({
    style: { background: "transparent", color: LIGHT_PALETTE.text, fontFamily: BODY_FONT },
    ...spec.plot,
  });
  const plotG = document.createElementNS(SVG_NS, "g");
  // Position the plot below the title block.
  const plotX = EXPORT_MARGIN.left;
  const plotY = EXPORT_MARGIN.top;
  plotG.setAttribute("transform", `translate(${plotX}, ${plotY})`);
  // Move all children of the inner SVG into the group, dropping the inner
  // <svg> wrapper itself (we are nesting inside our outer SVG).
  while (plotEl.firstChild) {
    plotG.appendChild(plotEl.firstChild);
  }
  svg.appendChild(plotG);

  // Legend rows (beneath the plot, inside the bottom margin).
  if (spec.legend && spec.legend.rows.length > 0) {
    const legendY = spec.height - 28;
    const legendG = document.createElementNS(SVG_NS, "g");
    legendG.setAttribute("transform", `translate(${EXPORT_MARGIN.left}, ${legendY})`);
    let cursorX = 0;
    for (const row of spec.legend.rows) {
      const item = renderLegendItem(row, cursorX);
      legendG.appendChild(item.group);
      cursorX += item.width + 16;
    }
    svg.appendChild(legendG);
  }

  return svg;
}

function renderLegendItem(row: LegendRow, x: number): { group: SVGGElement; width: number } {
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
  text.setAttribute("font", BODY_FONT);
  text.setAttribute("fill", LIGHT_PALETTE.text);
  text.textContent = row.label;
  g.appendChild(text);

  // Approximate width — symbol + gap + len(label)*7. Good enough for layout.
  const width = 16 + row.label.length * 7;
  return { group: g, width };
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
