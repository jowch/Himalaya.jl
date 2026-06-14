// traceColors.ts — the export figure's distinguishable per-trace palette, as
// literal sRGB hex. This is the readback of lib/comparison/coloring.ts
// COMPARE_PALETTE (OKLCH), which is tuned to be mutually distinct AND held ≥13°
// off every phase hue — so a categorical trace colour can never be mistaken for
// a phase encoding. Kept separate from cleanFigureSvg (the SVG builder) so it's
// plain colour data the key adapter can read without pulling in the builder.

/** 12-colour categorical palette, walked by trace order. */
export const COMPARE_PALETTE_HEX: readonly string[] = [
  "#b44c38", "#a56400", "#287c39", "#007e68", "#007c84", "#007798",
  "#096cb5", "#5a5eb2", "#635bb0", "#a04e94", "#b24866", "#7e7400",
];

/** Neutral cool-gray for orphan / unphased traces (COMPARE coloring ORPHAN). */
export const TRACE_ORPHAN_HEX = "#8a8f9c";

/** Phase → literal sRGB hex, the readback of phases.ts PHASE_PALETTE (OKLCH),
 *  so a phase-coloured export (the single-trace Focus figure) matches the
 *  on-screen trace hue. The standalone SVG can't carry oklch()/var(). */
export const PHASE_HEX: Record<string, string> = {
  Pn3m: "#b65b00",
  Im3m: "#007d53",
  Ia3d: "#7555a8",
  Fm3m: "#b54952",
  Fd3m: "#855095",
  Hexagonal: "#a44b79",
  Lamellar: "#375fb9",
  Square: "#548126",
};
/** Neutral fallback for an unphased trace under phase colouring. */
export const PHASE_FALLBACK_HEX = "#777777";

/** Resolve a phase's export hex (fallback for null / unknown). */
export function phaseHex(phase: string | null): string {
  return (phase != null && PHASE_HEX[phase]) || PHASE_FALLBACK_HEX;
}
