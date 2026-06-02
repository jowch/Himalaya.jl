import { positiveExtent } from "../plot/projection";

/** One predicted Bragg reflection in a phase's comb. */
export interface CombTooth {
  /** Predicted reflection q (Å⁻¹). */
  q: number;
  /** hkl or √N label shown above the tooth, e.g. "√6". */
  label: string;
  /** A claimed observed peak exists → solid tooth; otherwise a predicted-absent caret. */
  observed: boolean;
  /** Fractional residual (q_obs − q_pred)/q_pred for the claimed peak. Used by ResidualChart. */
  residual?: number;
}

/** One assigned (or hovered-preview) phase's comb row. `color` is a token ref
 *  (phaseColor(phase)) supplied by the composite. */
export interface CombSeries {
  phase: string;
  color: string;
  /** Comb gutter line 2, e.g. "a = 197 Å" (from lattice_d). */
  latticeLabel?: string;
  /** Residual gutter line 2 value, e.g. 0.998 (from r_squared). */
  rSquared?: number;
  teeth: CombTooth[];
}

/** A laid-out row. `assigned`/`preview` carry a series; `leftover` carries raw q-values. */
export type CombRow =
  | { kind: "assigned"; series: CombSeries }
  | { kind: "preview"; series: CombSeries }
  | { kind: "leftover"; qs: number[] };

/** Ordered rows: assigned phases, then the optional hovered-preview row, then the
 *  optional leftover (unindexed-peaks) row. ResidualChart passes `leftover: []`
 *  so it never gets a leftover row. */
export function assembleRows(
  assigned: CombSeries[],
  hovered: CombSeries | undefined,
  leftover: number[],
): CombRow[] {
  const rows: CombRow[] = assigned.map((series) => ({ kind: "assigned", series }));
  if (hovered) rows.push({ kind: "preview", series: hovered });
  if (leftover.length > 0) rows.push({ kind: "leftover", qs: leftover });
  return rows;
}

/** Padded log-q domain over every q in the rows (teeth + leftover). ~10% pad each
 *  side (mockup convention). Non-positive q are ignored by positiveExtent. */
export function combQDomain(rows: CombRow[]): [number, number] {
  const qs: number[] = [];
  for (const row of rows) {
    if (row.kind === "leftover") qs.push(...row.qs);
    else for (const t of row.series.teeth) qs.push(t.q);
  }
  const [lo, hi] = positiveExtent(qs, [0.01, 0.2]);
  return [lo * 0.9, hi * 1.1];
}
