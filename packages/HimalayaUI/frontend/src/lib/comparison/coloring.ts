/**
 * Coloring library for Compare-page traces (Plan §Phase 9, Task 9.1;
 * spec §Trace coloring).
 *
 * Three grouping modes:
 *   bySample  — color is a property of `sample_id`. Reordering doesn't
 *               shift colors. Same-sample multiple traces share a color
 *               (acceptable for v1; user breaks ties via `color_override`).
 *   byPhase   — color is `phaseColor(snapshot.confirmed_index.phase)`.
 *               Unindexed traces fall back to `ORPHAN_FALLBACK`.
 *   distinct  — palette walked by `display_order` (in `allMembers` array
 *               order). Cycles when N > palette.length.
 *
 * Resolution order: `member.color_override` → grouping-mode default →
 * `ORPHAN_FALLBACK`.
 *
 * The palette is a small categorical set tuned for dark backgrounds (chroma
 * 0.10–0.13, luminance 0.76–0.80). Hues are spread ≥30° apart to keep
 * adjacent palette entries visually distinct. Vendoring here rather than
 * reusing `phases.ts` because phase colors carry semantic meaning (phase A
 * is always color X) — categorical trace colors must not collide with that
 * mapping or `byPhase` and `bySample` would visually conflate.
 *
 * If a comparison surface ever needs cross-experiment color stability, we
 * can swap the palette indexer for a `hash(sample_id) % len` scheme without
 * changing call sites — the API takes a member + context and returns a
 * resolved color.
 */
import type { ComparisonMember } from "../../api";
import { phaseColor } from "../../phases";

export type GroupingMode = "bySample" | "byPhase" | "distinct";

/**
 * 12-color categorical palette (OKLCH). Hues spaced ~30° apart so two
 * adjacent palette entries are always distinguishable. Chroma + luminance
 * tuned to match the dark theme's text/background contrast band.
 */
export const COMPARE_PALETTE: readonly string[] = [
  "oklch(0.78 0.13  18)", // coral
  "oklch(0.80 0.13  62)", // amber
  "oklch(0.79 0.12 132)", // chartreuse
  "oklch(0.78 0.12 160)", // sage
  "oklch(0.79 0.12 185)", // seafoam teal
  "oklch(0.78 0.12 205)", // sky
  "oklch(0.80 0.10 248)", // periwinkle
  "oklch(0.76 0.12 280)", // indigo-violet
  "oklch(0.76 0.12 300)", // violet
  "oklch(0.76 0.12 318)", // rose-purple
  "oklch(0.78 0.12 348)", // pink
  "oklch(0.79 0.10  90)", // olive-yellow
];

/** Neutral cool-gray for orphans (NULL exposure, no phase, etc.). */
export const ORPHAN_FALLBACK = "oklch(0.65 0.02 270)";

/**
 * Context passed to `colorFor`. `allMembers` is needed for `bySample` and
 * `distinct` modes to derive a deterministic palette index. `sampleIdFor`
 * is supplied by the caller so `colorFor` itself stays pure — it doesn't
 * need to know about TanStack cache lookups.
 */
export interface ColorContext {
  /** All members in the same comparison, in display order. */
  allMembers: ReadonlyArray<ComparisonMember>;
  /**
   * Resolve a member's `sample_id`. Returns `null` when the exposure is
   * unknown (cache miss / orphan). Caller wires this against the TanStack
   * exposure cache or another lookup source.
   */
  sampleIdFor: (m: ComparisonMember) => number | null;
}

/**
 * Resolve the rendered color for one member under the given grouping mode.
 *
 * Always pure — the only side-effect is reading `context.sampleIdFor(member)`,
 * which is supplied by the caller. Suitable for use in a render path.
 */
export function colorFor(
  member: ComparisonMember,
  mode: GroupingMode,
  palette: readonly string[],
  context: ColorContext,
): string {
  // 1. Override always wins.
  if (member.color_override != null && member.color_override !== "") {
    return member.color_override;
  }

  // 2. Grouping-mode default.
  switch (mode) {
    case "bySample":
      return defaultBySample(member, palette, context);
    case "byPhase":
      return defaultByPhase(member);
    case "distinct":
      return defaultDistinct(member, palette, context);
    default:
      // Unknown mode (defensive — TS should prevent this at compile time).
      return ORPHAN_FALLBACK;
  }
}

function defaultBySample(
  member: ComparisonMember,
  palette: readonly string[],
  context: ColorContext,
): string {
  if (member.exposure_id === null) return ORPHAN_FALLBACK;
  const sampleId = context.sampleIdFor(member);
  if (sampleId === null) return ORPHAN_FALLBACK;

  // Per spec §Trace coloring: `palette[hash(sample_id) % palette.length]`.
  // Color is a property of the sample, NOT the slot — reordering must not
  // shift colors. A pure hash on sample_id achieves that without scanning
  // `allMembers`. (`allMembers` stays in the context for `distinct` mode.)
  const idx = hashSampleId(sampleId) % palette.length;
  return palette[idx] ?? ORPHAN_FALLBACK;
}

/**
 * Deterministic 32-bit FNV-1a hash of an integer sample id, returned as a
 * non-negative 31-bit integer for use as a palette index. Stable across
 * reloads and tab sessions; only depends on the input id, so the same
 * sample lands in the same palette bucket every time.
 */
function hashSampleId(id: number): number {
  // FNV-1a over the 4 little-endian bytes of the id (we treat ids as
  // unsigned for the hash). For palette sizes ≤16 this is plenty of
  // entropy to spread small id ranges (1, 2, 3, …) evenly.
  let h = 0x811c9dc5;
  let n = id >>> 0;
  for (let i = 0; i < 4; i++) {
    h ^= n & 0xff;
    h = Math.imul(h, 0x01000193) >>> 0;
    n >>>= 8;
  }
  // Drop the sign bit so `% palette.length` always returns a non-negative.
  return h & 0x7fffffff;
}

function defaultByPhase(member: ComparisonMember): string {
  const phase = member.snapshot?.confirmed_index?.phase;
  if (!phase) return ORPHAN_FALLBACK;
  return phaseColor(phase);
}

function defaultDistinct(
  member: ComparisonMember,
  palette: readonly string[],
  context: ColorContext,
): string {
  // Position in `allMembers` (caller already sorts by display_order).
  const idx = context.allMembers.findIndex((m) => m.id === member.id);
  if (idx < 0) return ORPHAN_FALLBACK;
  return palette[idx % palette.length] ?? ORPHAN_FALLBACK;
}
