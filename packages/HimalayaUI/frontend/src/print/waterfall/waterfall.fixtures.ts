import { toWaterfallRows, type WaterfallRow } from "./waterfallModel";
import {
  realMembers,
  transitionSeries,
  formFactorMember,
  unindexedMember,
} from "../fixtures/realSeriesMembers";
import { realTraces } from "../fixtures/realTraces";

/** The full real Sample-9 set (Ia3d → Im3m+Ia3d → Im3m + stress members). */
export const FULL: WaterfallRow[] = toWaterfallRows(realMembers, realTraces);

/** The three-member transition only — the canonical hero story. */
export const TRANSITION: WaterfallRow[] = toWaterfallRows(transitionSeries, realTraces);

/** Mixed states: an indexed member, a form-factor member, an unindexed member. */
// The synthetic form-factor / no-index members predate real trace capture
// (exposure_id null). Borrow a real measured trace so the story shows the
// realistic look — a curve with NO anchor beads — instead of a blank band.
const ff = toWaterfallRows([formFactorMember], realTraces)[0]!;
const ni = toWaterfallRows([unindexedMember], realTraces)[0]!;
export const MIXED_STATES: WaterfallRow[] = [
  toWaterfallRows([realMembers[0]!], realTraces)[0]!,
  { ...ff, trace: realTraces[66]! },
  { ...ni, trace: realTraces[67]! },
];
