/**
 * `isSeriesDraftDirty` — does an in-progress builder draft carry unsaved
 * changes vs the saved series? (BU-NAVAWAY-DRAFT.)
 *
 * Compares the normalized SAVE projection of both sides: `buildSeriesSaveBody`
 * strips local row ids, sorts by position, and canonicalizes empty strings, so
 * a pristine fork (Edit just opened, nothing changed) compares EQUAL and is not
 * dirty. Both bodies come from the same builder, so their key order is stable
 * and a JSON-string compare is a sound structural equality here.
 *
 * Drives the builder's `beforeunload` guard: the warning fires only on a real
 * unsaved change (the controls-don't-lie law), never on a no-op draft.
 */
import type { Series } from "../../api";
import type { SeriesDraft } from "./seriesDraft";
import { fromSeries } from "./seriesDraftFactories";
import { buildSeriesSaveBody } from "./buildSeriesSaveBody";

export function isSeriesDraftDirty(draft: SeriesDraft, series: Series): boolean {
  return (
    JSON.stringify(buildSeriesSaveBody(draft)) !==
    JSON.stringify(buildSeriesSaveBody(fromSeries(series)))
  );
}
