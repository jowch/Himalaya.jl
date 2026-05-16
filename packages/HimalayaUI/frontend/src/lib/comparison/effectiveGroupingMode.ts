/**
 * Canonical helper for resolving the active grouping mode
 * (Compare UX C-4 Step 0).
 *
 * Precedence (highest → lowest):
 *   1. Draft's `viewGroupingMode` — owns the value while an edit is active
 *      (incl. viewer escape hatch: toggling without an existing draft
 *      creates one, so the draft always exists when the toggle fires).
 *   2. Server record's `view_grouping_mode` — the author's last-saved pick.
 *   3. Hard default: `"bySample"`.
 *
 * Consumers never spell out the fallback chain themselves; this helper is
 * the single source of truth for that logic.
 */
import type { ActiveDraft } from "./draft";
import type { Comparison } from "../../api";
import type { GroupingMode } from "./coloring";

export type { GroupingMode };

export function effectiveGroupingMode(
  draft: ActiveDraft | null,
  comparison: Comparison | undefined,
): GroupingMode {
  // Draft owns the value during edit (incl. viewer escape hatch where
  // toggling without a draft creates one).
  if (draft?.viewGroupingMode !== undefined) return draft.viewGroupingMode;
  // Otherwise the server record is the source of truth.
  return (comparison?.view_grouping_mode as GroupingMode | null) ?? "bySample";
}
