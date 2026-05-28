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
 * The palette is a small categorical set paper-tuned for "The Print" (L
 * 0.52–0.56, chroma 0.12–0.14). Hues are spread ~30° apart to keep adjacent
 * palette entries visually distinct, and offset from the canonical phase
 * hues so `byPhase` and `bySample` never visually conflate. Vendoring here
 * rather than reusing `phases.ts` because phase colors carry semantic
 * meaning (phase A is always color X) — categorical trace colors must not
 * collide with that mapping.
 *
 * If a comparison surface ever needs cross-experiment color stability, we
 * can swap the palette indexer for a `hash(sample_id) % len` scheme without
 * changing call sites — the API takes a member + context and returns a
 * resolved color.
 */
import type { SeriesMember } from "../../api";
import { phaseColor } from "../../phases";

export type GroupingMode = "bySample" | "byPhase" | "distinct";

/**
 * 12-color categorical palette (OKLCH), paper-tuned for "The Print" (R8-N2,
 * round-2 finding) so bySample/distinct traces read at WCAG AA on `--plate`.
 *
 * Tuning: luminance L 0.52–0.56, chroma 0.12–0.14 — the same band the R0b
 * PHASE_PALETTE retune (#222) settled on. The retired dark-field tuning
 * (L 0.76–0.80) washed out on the warm paper surface — round-2 caught the
 * regression because the live builder defaults to `bySample` and most series
 * are unindexed, making this palette the dominant visual on the destination
 * surface.
 *
 * **Phase-offset invariant.** Each entry is held ≥13° from every hue in
 * `PHASE_PALETTE` (`phases.ts`), with nine of twelve sitting at exactly 15°.
 * At L≈0.55 C≈0.13 a 13° hue shift is ΔE2000 ≈ 2.5 — perceptually distinct,
 * the smallest gap the geometry allows while keeping twelve entries (the
 * eight phase hues pack the warm and purple sectors tightly enough that
 * ≥15° everywhere AND twelve mutually-distinct entries is infeasible — see
 * the proof in `test/coloring.test.ts`'s phase-offset block). The egregious
 * 3° collision the round-1 review caught (entry 315 vs Fd3m 318) has been
 * resolved by relocating that entry to 285.
 *
 * The offset matters because `byPhase` resolves through `PHASE_PALETTE` and
 * `bySample`/`distinct` walk this palette: a sample never indexed as Pn3m
 * must not render in something perceptually identical to Pn3m's amber.
 *
 * Tradeoff: hue 285 sits 6° from neighbour 279, so the two render as
 * "purple" and "violet-purple" in legend order — visually paired but
 * distinct. Acceptable; this palette is categorical, not a rainbow.
 *
 * Regression tests in `test/coloring.test.ts` pin the invariants: the
 * numeric phase-offset floor (≥13°), the L-band (0.50–0.58), and the AA-
 * on-`--plate` contrast floor mirroring `phases.test.ts`. Keep all three
 * green if you re-tune hues here.
 */
export const COMPARE_PALETTE: readonly string[] = [
  "oklch(0.55 0.14  33)", // warm coral
  "oklch(0.56 0.14  73)", // gold
  "oklch(0.52 0.13 147)", // moss
  "oklch(0.52 0.12 177)", // teal
  "oklch(0.52 0.12 200)", // cyan
  "oklch(0.52 0.12 220)", // azure
  "oklch(0.52 0.14 249)", // lavender
  "oklch(0.52 0.13 279)", // purple
  "oklch(0.52 0.13 285)", // violet-purple (was 315 = 3° from Fd3m 318; #251 r1)
  "oklch(0.55 0.14 333)", // raspberry
  "oklch(0.55 0.14   5)", // pink-red
  "oklch(0.55 0.13 105)", // chartreuse
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
  allMembers: ReadonlyArray<SeriesMember>;
  /**
   * Resolve a member's `sample_id`. Returns `null` when the exposure is
   * unknown (cache miss / orphan). Caller wires this against the TanStack
   * exposure cache or another lookup source.
   */
  sampleIdFor: (m: SeriesMember) => number | null;
}

/**
 * Resolve the rendered color for one member under the given grouping mode.
 *
 * Always pure — the only side-effect is reading `context.sampleIdFor(member)`,
 * which is supplied by the caller. Suitable for use in a render path.
 */
export function colorFor(
  member: SeriesMember,
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
  member: SeriesMember,
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

function defaultByPhase(member: SeriesMember): string {
  const phase = member.snapshot?.confirmed_index?.phase;
  if (!phase) return ORPHAN_FALLBACK;
  return phaseColor(phase);
}

function defaultDistinct(
  member: SeriesMember,
  palette: readonly string[],
  context: ColorContext,
): string {
  // Position in `allMembers` (caller already sorts by display_order).
  const idx = context.allMembers.findIndex((m) => m.id === member.id);
  if (idx < 0) return ORPHAN_FALLBACK;
  return palette[idx % palette.length] ?? ORPHAN_FALLBACK;
}
