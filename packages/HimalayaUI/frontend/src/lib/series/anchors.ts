/**
 * series/anchors — the `(phase, Miller-order)` anchor map driving the Series
 * migration-tracking layer (Plan E, Task E-2).
 *
 * The waterfall's job is "what changes across the variable?". The migration
 * track answers it by connecting the SAME reflection — keyed by `(phase,
 * order)` — across every member that carries the phase, so the connector reads
 * as a peak drifting in q as the ordering variable changes.
 *
 * Crucially, a member can carry a phase yet NOT observe a particular Miller
 * order (the √3 Pn3m reflection is weak and routinely absent). We must mark
 * that order as *predicted-but-absent* — a hollow ghost ring at the predicted
 * q — rather than silently skipping the member, so the track threads through a
 * visible "missing peak" instead of jumping over it.
 *
 * **Predicted-q physics (client-derived, no backend).** For order k the
 * predicted q is `2π·√(radicand[k]) / lattice_d` for cubics and `2π·n / d`
 * for lamellar — the SAME relation as `predictpeaks` / the Focus combs
 * (`basis × phaseratios(P)`, with `basis = 2π / lattice_d`). The radicands
 * come from the shared `RADICANDS` table (mirrors `src/phase.jl phaseratios`).
 *
 * Observed-vs-absent is decided by matching each predicted-q against the
 * member's `effective_peaks` within a relative tolerance: a predicted order
 * with a nearby observed peak is "present" (carries the observed q); one with
 * none is "absent" (q === null, predictedQ set).
 */
import type { SeriesMember } from "../../api";
import { RADICANDS } from "../seriesRatio";

const TWO_PI = 2 * Math.PI;

/** Number of Miller orders to track per phase — matches the cross-trace layer's
 *  MAX_MILLER_ORDERS cap (deeper orders crowd the track and rarely add insight). */
const MAX_ORDERS = 4;

/** Relative q tolerance for matching a predicted order to an observed peak.
 *  Loose enough to absorb lattice-fit residual, tight enough not to bind an
 *  adjacent order. */
const MATCH_REL_TOL = 0.02;

/**
 * Predicted q for Miller order `order` (0-based) of `phase` at lattice
 * parameter `latticeD`. Cubic: `2π·√(radicand) / a`; lamellar: `2π·n / d`
 * (radicand[k] = n² for lamellar, so √ is exact). Returns `null` for a phase
 * not in the radicand table or an out-of-range order.
 */
export function predictedQForOrder(
  phase: string,
  latticeD: number,
  order: number,
): number | null {
  const radicands = RADICANDS[phase];
  if (!radicands || order < 0 || order >= radicands.length) return null;
  if (!(latticeD > 0)) return null;
  return (TWO_PI * Math.sqrt(radicands[order]!)) / latticeD;
}

/** One vertex in a `(phase, order)` track. `q` is the OBSERVED q when the order
 *  is present, or `null` when the phase is present but the order is absent
 *  (predicted-but-not-observed). `predictedQ` is always the comb position so
 *  the layer can plant the ghost ring there. */
export interface AnchorVertex {
  /** Member position in the render order (aligned with the `members` array). */
  memberPos: number;
  /** Observed q, or `null` when the order is predicted-but-absent. */
  q: number | null;
  /** Predicted comb position (2π·√radicand / a) for this order, always set. */
  predictedQ: number;
  /** True when this order is predicted but not observed in this member. */
  absent: boolean;
}

/** `"phase:order"` → vertices across carrying members, in member-position order. */
export type AnchorMap = Map<string, AnchorVertex[]>;

/**
 * Build the `(phase, order)` anchor map from member snapshots. For each member
 * with a `confirmed_index`, derive the predicted-q comb from `lattice_d`, then
 * for each order decide present (a nearby observed peak) vs absent. Members
 * without a confirmed index (form-factor / null) contribute nothing.
 */
export function buildAnchorMap(members: SeriesMember[]): AnchorMap {
  const map: AnchorMap = new Map();

  for (let pos = 0; pos < members.length; pos++) {
    const m = members[pos]!;
    const ci = m.snapshot?.confirmed_index;
    if (!ci) continue; // form-factor / null / unindexed → no anchors
    const phase = ci.phase;
    const radicands = RADICANDS[phase];
    if (!radicands || ci.lattice_d == null || !(ci.lattice_d > 0)) continue;

    const observed = (m.snapshot?.effective_peaks ?? []).map((p) => p.q);
    const maxOrders = Math.min(MAX_ORDERS, radicands.length);

    for (let k = 0; k < maxOrders; k++) {
      const predictedQ = predictedQForOrder(phase, ci.lattice_d, k);
      if (predictedQ === null) continue;
      // Match against the nearest observed peak within relative tolerance.
      let observedQ: number | null = null;
      let best = Infinity;
      for (const oq of observed) {
        const rel = Math.abs(oq - predictedQ) / predictedQ;
        if (rel <= MATCH_REL_TOL && rel < best) {
          best = rel;
          observedQ = oq;
        }
      }
      const key = `${phase}:${k}`;
      if (!map.has(key)) map.set(key, []);
      map.get(key)!.push({
        memberPos: pos,
        q: observedQ,
        predictedQ,
        absent: observedQ === null,
      });
    }
  }

  return map;
}
