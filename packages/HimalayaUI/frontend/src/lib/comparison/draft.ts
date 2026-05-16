/**
 * Compare draft type definitions and sessionStorage persistence helpers
 * (Plan §Phase 4, Task 4.3).
 *
 * Optionality is encoded as `T | undefined` (not `T?`) per CLAUDE.md
 * frontend gotchas — `exactOptionalPropertyTypes: true` blocks
 * `set({ field: undefined })` on `field?: T` declarations, which we need for
 * round-trip Zustand updates.
 *
 * baseHash invariant: captured at edit-mode entry and FROZEN until submit.
 * SSE foreign events do NOT refresh it; submission rides this baseline as
 * `expected_content_hash` so a real conflict surfaces instead of being
 * silently swallowed.
 *
 * Brand-new comparison: `id` and `baseHash` are both undefined; submit uses
 * POST /api/comparisons (omit `expected_content_hash`). After the first
 * successful submit the draft is normally discarded; subsequent edits
 * load with `baseHash` set to the server's `content_hash`.
 *
 * **Module structure note:** the snapshot-deriving factories live in
 * `draftFactories.ts` to avoid the import cycle
 *   state.ts → draft.ts → snapshot.ts → queries.ts → state.ts.
 * `state.ts` calls `loadDraftFromSession()` at module-init time, so this
 * file MUST stay free of any transitive `queries.ts` import.
 */
import type { MemberSnapshot } from "../../api";
import type { GroupingMode } from "./coloring";

export type DraftMemberNormalization = "none" | "max" | "area" | "qwindow";

export interface DraftMember {
  /** Existing server-side id (UPDATE on submit) or undefined (INSERT). */
  id: number | undefined;
  /** Null = orphaned member (exposure deleted; FK ON DELETE SET NULL). */
  exposure_id: number | null;
  display_order: number;
  band_height: number;
  y_offset: number;
  normalization: DraftMemberNormalization;
  color_override: string | undefined;
  label_override: string | undefined;
  q_window_min: number | undefined;
  q_window_max: number | undefined;
  peak_display: { hidden: number[]; labeled: number[] } | undefined;
  /**
   * Computed at submit time via `computeMemberSnapshot`; undefined during
   * editing. Stored only after `loadDraftFromComparison` recovery so the
   * UI has something to render before the next mutation.
   */
  snapshot: MemberSnapshot | undefined;
}

export interface ActiveDraft {
  /** Server-assigned comparison id, or undefined for an unsubmitted draft. */
  id: number | undefined;
  /** content_hash captured at edit-mode entry. Frozen until submit. */
  baseHash: string | undefined;
  title: string;
  description: string;
  members: DraftMember[];
  /**
   * Fork lineage (Plan §Phase 11). Set when this draft was started as a fork
   * of an existing comparison; both ride together to `POST /api/comparisons`
   * so the backend can record the parent + the parent's content_hash at
   * fork time. Undefined for non-fork create flows.
   */
  forkedFromId: number | undefined;
  forkedAtHash: string | undefined;
  /**
   * Author's view choices (spec §6.4; Compare UX C-4).
   * `undefined` = inherit from server record / default.
   * These shadow the server's `view_grouping_mode`, `view_show_peak_ticks`,
   * `view_show_peak_labels` while a draft is active; forwarded to the
   * save payload so the server persists the author's latest preference.
   */
  viewGroupingMode: GroupingMode | undefined;
  viewShowPeakTicks: boolean | undefined;
  viewShowPeakLabels: boolean | undefined;
}

export type ActiveDraftSlot = ActiveDraft | null;

export const COMPARE_DRAFT_KEY = "himalaya-ui:compare-draft";
export const COMPARE_DRAFT_VERSION = 1;

interface PersistedEnvelope {
  version: number;
  draft: ActiveDraft;
}

export function persistDraftToSession(draft: ActiveDraft | null): void {
  try {
    if (draft === null) {
      sessionStorage.removeItem(COMPARE_DRAFT_KEY);
      return;
    }
    const env: PersistedEnvelope = { version: COMPARE_DRAFT_VERSION, draft };
    sessionStorage.setItem(COMPARE_DRAFT_KEY, JSON.stringify(env));
  } catch {
    // sessionStorage can throw under quota/access restrictions; the draft
    // is best-effort durable, so swallow and continue.
  }
}

export function loadDraftFromSession(): ActiveDraft | null {
  try {
    const raw = sessionStorage.getItem(COMPARE_DRAFT_KEY);
    if (!raw) return null;
    const parsed = JSON.parse(raw) as Partial<PersistedEnvelope> | null;
    if (!parsed || typeof parsed !== "object") return null;
    if (parsed.version !== COMPARE_DRAFT_VERSION) return null;
    if (!parsed.draft || typeof parsed.draft !== "object") return null;
    return parsed.draft as ActiveDraft;
  } catch {
    return null;
  }
}

export function emptyDraft(): ActiveDraft {
  return {
    id: undefined,
    baseHash: undefined,
    title: "",
    description: "",
    members: [],
    forkedFromId: undefined,
    forkedAtHash: undefined,
    viewGroupingMode: undefined,
    viewShowPeakTicks: undefined,
    viewShowPeakLabels: undefined,
  };
}
