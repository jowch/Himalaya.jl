import type { Trace } from "../../api";
import t37 from "./traces/37.json";
import t65 from "./traces/65.json";
import t66 from "./traces/66.json";
import t67 from "./traces/67.json";
import t93 from "./traces/93.json";

/** Measured traces (q, I, sigma) keyed by exposure id. Feeds the measured-trace
 *  TraceCurve consumers (e.g. the scoping sparkline). */
export const realTraces: Record<number, Trace> = {
  37: t37 as Trace, 65: t65 as Trace, 66: t66 as Trace, 67: t67 as Trace, 93: t93 as Trace,
};
