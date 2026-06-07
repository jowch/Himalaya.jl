import type { ReactNode } from "react";
import type { SeriesMember, SeriesSummary } from "../../api";
import type { PhaseSegment } from "../ui/PhaseStrip";
import { dominantPhase, resolveState } from "../waterfall/waterfallModel";

/** Per-member phase-strip segments — SAME phase/state derivation as toWaterfallRows
 *  (imported, not re-derived) so a card's strip and waterfall never disagree.
 *  Coexistence = the confirmed_phases tail beyond the dominant phase. */
export function membersToSegments(members: SeriesMember[]): PhaseSegment[] {
  return members.map((m) => {
    const phase = dominantPhase(m);
    const state = resolveState(m);
    const seg: PhaseSegment = { phase };
    // Coexistence only applies to an indexed (non-null) dominant phase — a
    // form_factor/null member must never carry a stray coexistence array.
    if (phase !== null) {
      const cp = m.snapshot?.confirmed_phases ?? [];
      if (cp.length > 1) seg.coexistWith = cp.slice(1).map((p) => p.phase);
    }
    if (state === "form_factor" || state === "null") seg.state = state;
    return seg;
  });
}

export interface CardChrome {
  figLabel: string; title: string; sampleCount: number; variable: string;
  provenance: ReactNode; editedLabel: string; author: string;
  notice?: { tone: "draft" }; draft: boolean;
}

/** Everything a card shows that is derivable from the LIST summary (no detail fetch).
 *  `position` = 1-based index in the filtered/sorted list (the "Fig. N"). */
export function toCardChrome(s: SeriesSummary, position: number, now: Date): CardChrome {
  const draft = s.content_hash === "";
  const title = s.title.trim() === "" ? "Untitled series" : s.title;
  const chrome: CardChrome = {
    figLabel: draft ? "Recipe" : `Fig. ${position}`,
    title,
    sampleCount: s.member_count,
    variable: s.ordering_variable ?? "",
    provenance: s.spans_experiments ? "↔ cross-experiment · q normalized" : null,
    editedLabel: formatEdited(s.updated_at ?? s.last_event_at, now),
    author: s.author_username ?? "",
    draft,
  };
  if (draft) chrome.notice = { tone: "draft" };
  return chrome;
}

/** Relative-time label (ported verbatim from the legacy SeriesFolioCard). */
export function formatEdited(iso: string | null, now: Date): string {
  if (iso === null) return "recently";
  const then = new Date(iso.replace(" ", "T") + (iso.endsWith("Z") ? "" : "Z"));
  if (Number.isNaN(then.getTime())) return "recently";
  const days = Math.floor((now.getTime() - then.getTime()) / 86_400_000);
  if (days <= 0) return "just now";
  if (days === 1) return "yesterday";
  if (days < 7) return `${days} days ago`;
  if (days < 14) return "1 week ago";
  if (days < 21) return "2 weeks ago";
  return `${Math.floor(days / 7)} weeks ago`;
}
