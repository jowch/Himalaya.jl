// seriesFigureKeys.ts — turn the series' members into the export figure's
// per-trace KEY entries (one per member, parallel to the waterfall rows): a
// distinguishable categorical colour + the sample name + every assigned phase
// with its lattice (a/d) and, for bicontinuous cubics, κ. Coexistence lists all
// phases with equal billing — there is no "dominant" phase. Pure + testable;
// the page supplies the exposure→sample-name resolver and the q-units.
import type { SeriesMember, MemberSnapshotPhase } from "../../api";
import { CUBIC_PHASES } from "../../phases";
import { resolveState } from "../../lib/series/memberRead";
import { latticeUnitFromQUnits, inverseSquareUnits, formatKappaPretty } from "../../lib/units";
import { kappaForPhase } from "../../lib/curvature";
import { COMPARE_PALETTE_HEX } from "./traceColors";
import type { FigureTraceKey, FigureKeySegment } from "./cleanFigureSvg";

export interface SeriesFigureKeyContext {
  /** Resolve a member's exposure to its sample's display name. Returns null on
   *  a cache miss; the caller falls back to the member's own row label. */
  sampleNameForExposure: (exposureId: number | null) => string | null;
  /** Per-member fallback label (memberRowLabel) when the sample is unresolved. */
  fallbackLabel: (m: SeriesMember) => string;
  /** Experiment q-units, for the lattice unit (Å / nm) and κ unit (Å⁻² / nm⁻²). */
  qUnits: string | null | undefined;
}

/** Lattice symbol: cubics report the lattice parameter `a`; lamellar / hexagonal
 *  report the d-spacing `d` (mirrors seriesReading.latticeSymbol). */
function latticeSymbol(phase: string): "a" | "d" {
  return CUBIC_PHASES.has(phase) ? "a" : "d";
}

function fmtLattice(v: number | null): string {
  return v == null ? "—" : Number.isInteger(v) ? String(v) : v.toFixed(1);
}

/** The phases a member carries (confirmed_phases, else confirmed_index), honouring
 *  the 3-state — form-factor / null members carry none. */
function phasesOf(m: SeriesMember): MemberSnapshotPhase[] {
  const snap = m.snapshot;
  if (!snap) return [];
  if (resolveState(m) !== "indexed") return [];
  if (snap.confirmed_phases && snap.confirmed_phases.length > 0) return snap.confirmed_phases;
  if (snap.confirmed_index) {
    return [{ phase: snap.confirmed_index.phase, lattice_d: snap.confirmed_index.lattice_d }];
  }
  return [];
}

/** Build the `detail` string for one phase: "a = 252 Å · κ = 1.70×10⁻⁴ Å⁻²".
 *  κ is appended only for the bicontinuous cubics; no lattice → empty detail. */
function segmentDetail(phase: string, latticeD: number | null, qUnits: string | null | undefined): string {
  if (latticeD == null) return "";
  const lat = `${latticeSymbol(phase)} = ${fmtLattice(latticeD)} ${latticeUnitFromQUnits(qUnits)}`;
  const k = kappaForPhase(phase, latticeD);
  if (k == null) return lat;
  return `${lat} · κ = ${formatKappaPretty(k)} ${inverseSquareUnits(qUnits)}`;
}

/** One phase → its figure key segment (phase name + formatted lattice/κ detail).
 *  Shared by the series adapter and the single-trace Focus figure. */
export function phaseSegment(
  phase: string,
  latticeD: number | null,
  qUnits: string | null | undefined,
): FigureKeySegment {
  return { phase, detail: segmentDetail(phase, latticeD, qUnits) };
}

/** The phaseless note for a non-indexed member. */
function phaselessNote(m: SeriesMember): string {
  return resolveState(m) === "form_factor" ? "form factor (no Bragg)" : "unindexed";
}

/** One FigureTraceKey per member, in input order (parallel to the waterfall rows). */
export function buildSeriesFigureKeys(
  members: SeriesMember[],
  ctx: SeriesFigureKeyContext,
): FigureTraceKey[] {
  return members.map((m, i) => {
    const phases = phasesOf(m);
    const segments: FigureKeySegment[] = phases.map((p) =>
      phaseSegment(p.phase, p.lattice_d, ctx.qUnits),
    );
    const label = ctx.sampleNameForExposure(m.exposure_id) ?? ctx.fallbackLabel(m);
    const key: FigureTraceKey = {
      color: COMPARE_PALETTE_HEX[i % COMPARE_PALETTE_HEX.length]!,
      label,
      segments,
    };
    if (segments.length === 0) key.note = phaselessNote(m);
    return key;
  });
}
