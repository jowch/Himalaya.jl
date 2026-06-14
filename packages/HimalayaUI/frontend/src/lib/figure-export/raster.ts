// raster.ts — turn a standalone SVG markup string into a downloadable blob / PNG.
// The export engine now produces SVG markup directly (the greenfield CleanFigure
// builder), so the download path is engine-agnostic: blob it for SVG, or
// rasterize it through an OffscreenCanvas for a 2× PNG.
import { canExportPng } from "./renderer";

/** Whether the 2× PNG rasterization pipeline is available (OffscreenCanvas +
 *  convertToBlob + Image). A plain function (not a re-export) so it is spyable. */
export function canRasterizePng(): boolean {
  return canExportPng();
}

export function svgStringToBlob(svg: string): Blob {
  return new Blob([svg], { type: "image/svg+xml" });
}

function dim(svg: string, attr: "width" | "height", fallback: number): number {
  const m = new RegExp(`<svg[^>]*\\b${attr}="(\\d+(?:\\.\\d+)?)"`).exec(svg);
  return m ? Number(m[1]) : fallback;
}

/** Target raster density for the PNG export. The SVG authors its geometry in
 *  user units treated as points (72/inch), so the canvas scale is dpi/72 — 150
 *  DPI (scale ≈ 2.08) gives a crisp figure for slides and print. */
export const EXPORT_PNG_DPI = 150;

/**
 * Render an SVG string to a PNG blob at EXPORT_PNG_DPI: SVG → blob URL → Image
 * decode → OffscreenCanvas drawImage → convertToBlob. `URL.revokeObjectURL` runs
 * in a `finally` so a thrown render still cleans up. JSDOM lacks these canvas
 * APIs — the full path is covered by the Playwright E2E.
 */
export async function svgStringToPng(svg: string, dpi = EXPORT_PNG_DPI): Promise<Blob> {
  const scale = dpi / 72;
  const w = dim(svg, "width", 800);
  const h = dim(svg, "height", 600);
  const url = URL.createObjectURL(svgStringToBlob(svg));
  try {
    const img = new Image();
    img.src = url;
    await img.decode();
    const off = new OffscreenCanvas(w * scale, h * scale);
    const ctx = off.getContext("2d");
    if (!ctx) throw new Error("OffscreenCanvas 2d context unavailable");
    ctx.drawImage(img, 0, 0, off.width, off.height);
    return await off.convertToBlob({ type: "image/png" });
  } finally {
    URL.revokeObjectURL(url);
  }
}
