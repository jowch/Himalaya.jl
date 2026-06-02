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
export const MIXED_STATES: WaterfallRow[] = toWaterfallRows(
  [realMembers[0]!, formFactorMember, unindexedMember],
  realTraces,
);
