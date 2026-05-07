/**
 * WarmAddMenu (Plan §Phase 5, Task 5.3).
 *
 * "Add to comparison" affordance shown on the Inspect page so the user can
 * push the active exposure into a draft comparison without leaving Inspect.
 *
 * v1 menu surface (per the spec, simplified):
 *   - "+ New comparison"
 *       Discard any existing draft, start a fresh one pre-populated with
 *       the active exposure, navigate to /experiments/:eid/compare/new.
 *   - "Recent draft: <title>"
 *       Visible only when activeDraft is non-null AND the active exposure
 *       isn't already a member. Adds the exposure to the existing draft
 *       in-tab — no SSE, no server round-trip. Other tabs see no change
 *       until the user submits and the SSE event fires (per spec).
 *   - "Already in current draft" hint
 *       Replaces the recent-draft button when the exposure is already a
 *       member, so the user knows they can't double-add.
 *
 * Deferred: "Pick a comparison..." (a separate comparison-scoped picker)
 * — covered by Phase 5 spec but flagged as v1-optional in the plan. The
 * "+ New" + "Recent draft" pair already delivers the warm-add value; the
 * comparison picker can ship in a follow-up alongside the comparison
 * sidebar improvements in Phase 11.
 *
 * Cross-tab boundary: drafts live in sessionStorage (per-tab) and are NOT
 * broadcast across tabs. A second tab editing the same comparison does
 * not see warm-add additions until the user submits and the SSE event
 * fires. (Cross-tab draft sync would need BroadcastChannel or server-side
 * drafts; deferred to v2 per spec §Warm-add affordance.)
 */
import { useEffect, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { useAppState } from "../state";
import { comparePath } from "../lib/comparison/routes";

interface Props {
  /** The exposure currently being viewed in Inspect; undefined disables the trigger. */
  exposureId: number | undefined;
  /** Active experiment context for navigation to /experiments/:eid/compare/new. */
  experimentId: number | undefined;
}

export function WarmAddMenu({
  exposureId,
  experimentId,
}: Props): JSX.Element {
  const navigate = useNavigate();
  const qc = useQueryClient();
  const draft = useAppState((s) => s.activeDraft);
  const startNewDraft = useAppState((s) => s.startNewDraft);
  const addMember = useAppState((s) => s.addMember);
  const setDraftTitle = useAppState((s) => s.setDraftTitle);

  const [open, setOpen] = useState(false);
  const containerRef = useRef<HTMLDivElement>(null);

  // Click-outside dismissal so the menu doesn't linger after navigation
  // intent shifts. Esc also closes (handled inline below).
  useEffect(() => {
    if (!open) return;
    const handler = (e: MouseEvent): void => {
      if (!containerRef.current) return;
      if (!containerRef.current.contains(e.target as Node)) {
        setOpen(false);
      }
    };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [open]);

  const triggerDisabled = exposureId === undefined;
  const alreadyAdded =
    exposureId !== undefined &&
    !!draft?.members.some((m) => m.exposure_id === exposureId);

  const onPickNew = (): void => {
    if (exposureId === undefined) return;
    // Start a fresh draft, pre-populated with the active exposure. We
    // explicitly call startNewDraft → addMember so existing draft (if
    // any) is discarded. The user can name it on the edit page.
    startNewDraft();
    addMember(exposureId, qc);
    // Default title — spec doesn't pin a format; use a date-ish handle so
    // the sidebar list isn't full of bare "Untitled" rows.
    const now = new Date().toISOString().slice(0, 10);
    setDraftTitle(`Comparison ${now}`);
    setOpen(false);
    // No experiment context → use the global create route
    // (/compare/all/new) so the user has somewhere to land. Comparisons
    // have no FK to experiments per spec, so the global vs experiment
    // scope only affects which sidebar listing they end up in.
    navigate(
      comparePath({
        scope: experimentId !== undefined ? "experiment" : "all",
        eid: experimentId,
        isNew: true,
      }),
    );
  };

  const onPickRecent = (): void => {
    if (exposureId === undefined) return;
    if (alreadyAdded) return;
    addMember(exposureId, qc);
    setOpen(false);
  };

  return (
    <div
      ref={containerRef}
      className="relative inline-block"
      onKeyDown={(e) => { if (e.key === "Escape") setOpen(false); }}
    >
      <button
        type="button"
        data-testid="warm-add-trigger"
        disabled={triggerDisabled}
        onClick={() => setOpen((v) => !v)}
        className="px-2 py-1 rounded border border-border text-sm text-fg-muted
                   hover:text-fg disabled:opacity-50"
      >
        Add to comparison ▾
      </button>

      {open && (
        <div
          data-testid="warm-add-menu"
          role="menu"
          className="absolute z-40 right-0 mt-1 min-w-[240px]
                     bg-bg-elevated border border-border rounded-md shadow-xl
                     py-1 anim-pal-in"
        >
          {draft !== null && alreadyAdded && (
            <div
              data-testid="warm-add-already-added"
              className="px-3 py-1.5 text-xs text-fg-dim italic"
            >
              Already in current draft
            </div>
          )}
          {draft !== null && !alreadyAdded && (
            <button
              type="button"
              role="menuitem"
              data-testid="warm-add-recent"
              onClick={onPickRecent}
              className="w-full text-left px-3 py-1.5 text-sm text-fg
                         hover:bg-bg-hover"
            >
              Add to current draft: <span className="font-medium">
                {draft.title || "(untitled)"}
              </span>
            </button>
          )}
          <button
            type="button"
            role="menuitem"
            data-testid="warm-add-new"
            onClick={onPickNew}
            className="w-full text-left px-3 py-1.5 text-sm text-fg
                       hover:bg-bg-hover"
          >
            + New comparison
          </button>
        </div>
      )}
    </div>
  );
}

