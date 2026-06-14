// raster.ts — turn a standalone SVG markup string into a downloadable blob / PNG.
// The export engine produces SVG markup directly (the greenfield cleanFigureSvg
// builder), so the download path is engine-agnostic: blob it for SVG, or
// rasterize it through an OffscreenCanvas for a high-DPI PNG.

/**
 * Whether the PNG rasterization pipeline is available (OffscreenCanvas +
 * convertToBlob + Image). The Save/Copy buttons use it to disable themselves on
 * browsers lacking the pipeline, so the user never clicks into a guaranteed
 * failure. Mirrors `canCopyPngToClipboard()` in clipboard.ts.
 */
export function canRasterizePng(): boolean {
  if (typeof OffscreenCanvas === "undefined") return false;
  // convertToBlob landed alongside OffscreenCanvas in evergreen browsers but is
  // missing on some Safari versions where OffscreenCanvas exists.
  const proto = (OffscreenCanvas as unknown as { prototype?: { convertToBlob?: unknown } }).prototype;
  if (!proto || typeof proto.convertToBlob !== "function") return false;
  if (typeof Image === "undefined") return false;
  return true;
}

export function svgStringToBlob(svg: string): Blob {
  return new Blob([svg], { type: "image/svg+xml" });
}

function dim(svg: string, attr: "width" | "height", fallback: number): number {
  const m = new RegExp(`<svg[^>]*\\b${attr}="(\\d+(?:\\.\\d+)?)"`).exec(svg);
  return m ? Number(m[1]) : fallback;
}

/** Target raster density for the PNG export. The SVG authors its geometry in
 *  points (72/inch) at the figure's true physical size, so the canvas scale is
 *  dpi/72 — 216 DPI (scale = 3×) upscales the compact point layout to a crisp
 *  raster for slides and print, while the SVG download stays the true size. */
export const EXPORT_PNG_DPI = 216;

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
