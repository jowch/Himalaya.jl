// measureText.ts — measure a label's rendered width so peak-label dodging only
// kicks in when labels would ACTUALLY overlap (vs a fixed worst-case width that
// dodges even well-spaced labels). Uses a memoized 2D canvas (the standard
// text-metrics path); falls back to a monospace character estimate where canvas
// is unavailable (jsdom under Vitest) so the dodge math stays deterministic in
// tests.

let _ctx: CanvasRenderingContext2D | null | undefined;

function ctx(): CanvasRenderingContext2D | null {
  if (_ctx !== undefined) return _ctx;
  _ctx =
    typeof document !== "undefined"
      ? document.createElement("canvas").getContext("2d")
      : null;
  return _ctx;
}

export interface TextFont {
  px: number;
  weight: number;
  /** A concrete font family (e.g. "monospace"); a CSS var won't resolve in a
   *  canvas font string. */
  family: string;
}

/**
 * Width in px of `text` rendered in `font`. Canvas `measureText` when available,
 * else a monospace estimate (the labels ARE a monospace font, so character
 * count × an em fraction is close). Never returns 0 for non-empty text.
 */
export function measureTextWidth(text: string, font: TextFont): number {
  if (text === "") return 0;
  const c = ctx();
  if (c) {
    c.font = `${font.weight} ${font.px}px ${font.family}`;
    const w = c.measureText(text).width;
    if (w > 0) return w;
  }
  // Monospace estimate: ~0.6em per character is a good bold-mono average.
  return text.length * font.px * 0.6;
}
