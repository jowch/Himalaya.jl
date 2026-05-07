/**
 * ComparePageEdit — edit-mode shell with save flow (Plan §Phase 4, Task 4.4).
 *
 * Reads URL params via react-router:
 *   /experiments/:eid/compare/new        — empty draft (create flow)
 *   /experiments/:eid/compare/:id/edit   — edit existing comparison
 *
 * Owns:
 *   - Title + description inputs (controlled, bound to Zustand draft).
 *   - Save / Cancel / Discard buttons.
 *   - Snapshot recompute at submit time per Plan §Task 4.3.
 *
 * The plot, member panel, drag handles, picker chip, and conflict modal
 * land in Phases 5–12; this file only mounts placeholders for them.
 *
 * Save flow (Plan §Task 4.4):
 *   - Disabled when `draft.members.length === 0`.
 *   - Computes a fresh `MemberSnapshot` per member at submit (cache-derived).
 *   - Calls `useSaveComparison`. The mutator routes to POST /api/comparisons
 *     (create) or POST /api/comparisons/:id/submit (update) based on the
 *     `id` field; `expected_content_hash` rides from `draft.baseHash`.
 *   - On success: clear draft, navigate to /experiments/:eid/compare/:newId.
 *   - On 409 ConflictError: leave the draft intact. The conflict modal in
 *     Phase 12 picks up the typed throw via `save.error`.
 *
 * The hydration effect aligns the in-memory draft with the URL on every
 * mount: `/new` ⇒ start a fresh empty draft if one isn't already active;
 * `/:id/edit` ⇒ load from the comparison fetch (with cache-derived snapshot
 * recovery) if the draft isn't already this id.
 */
import { useCallback, useEffect, useRef, useState } from "react";
import { useParams, useNavigate } from "react-router-dom";
import { useQueryClient } from "@tanstack/react-query";
import { ComparisonSidebar } from "../components/ComparisonSidebar";
import { ComparisonPicker } from "../components/ComparisonPicker";
import { useAppState } from "../state";
import { useSaveComparison, useComparison } from "../queries";
import { computeMemberSnapshot } from "../lib/comparison/snapshot";
import type { Comparison, ComparisonMemberInput, SaveComparisonBody } from "../api";

export function ComparePageEdit(): JSX.Element {
  const params = useParams<{ eid?: string; id?: string }>();
  const navigate = useNavigate();
  const qc = useQueryClient();

  const eid = params.eid !== undefined ? Number(params.eid) : undefined;
  const id  = params.id  !== undefined ? Number(params.id)  : undefined;

  const draft = useAppState((s) => s.activeDraft);
  const startNewDraft       = useAppState((s) => s.startNewDraft);
  const loadDraftFromComp   = useAppState((s) => s.loadDraftFromComparison);
  const setDraftTitle       = useAppState((s) => s.setDraftTitle);
  const setDraftDescription = useAppState((s) => s.setDraftDescription);
  const discardDraft        = useAppState((s) => s.discardDraft);

  // Hydrate draft from URL on mount / when URL id changes.
  const comparisonQ = useComparison(id);
  useEffect(() => {
    if (id === undefined) {
      // /new: start an empty draft if there isn't one (or one that's tied
      // to a different comparison id). Don't clobber a fresh draft mid-edit.
      if (draft === null || draft.id !== undefined) {
        startNewDraft();
      }
    } else if (comparisonQ.data && (draft === null || draft.id !== id)) {
      // /:id/edit: load from the server fetch with cache-derived snapshots.
      loadDraftFromComp(comparisonQ.data, qc);
    }
    // Intentionally exclude `draft` from deps — we only re-evaluate when
    // the URL id or the fetched comparison changes. Reading `draft` inside
    // the effect for the comparison-id check is fine; including it would
    // re-run the effect after the load and re-load again.
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [id, comparisonQ.data, startNewDraft, loadDraftFromComp, qc]);

  const save = useSaveComparison();
  const pendingSubmitRef = useRef(false);

  const goToReview = useCallback(
    (newId: number) => {
      if (eid !== undefined) navigate(`/experiments/${eid}/compare/${newId}`);
      else navigate(`/compare/all`);
    },
    [navigate, eid],
  );

  const goToList = useCallback(() => {
    if (eid !== undefined) navigate(`/experiments/${eid}/compare`);
    else navigate(`/compare/all`);
  }, [navigate, eid]);

  const handleCancel = useCallback(() => {
    if (id !== undefined) {
      // Cancel from edit-existing → return to that comparison's review page.
      if (eid !== undefined) navigate(`/experiments/${eid}/compare/${id}`);
      else navigate(`/compare/all`);
    } else {
      // Cancel from create → return to list.
      goToList();
    }
  }, [navigate, id, eid, goToList]);

  const handleDiscard = useCallback(() => {
    discardDraft();
    goToList();
  }, [discardDraft, goToList]);

  const handleSave = useCallback(() => {
    if (draft === null) return;
    if (draft.members.length === 0) return;
    // Compute a fresh snapshot per member at submit time (Plan §Task 4.3).
    const members: ComparisonMemberInput[] = draft.members.map((m) => {
      const snapshot = m.exposure_id !== null
        ? computeMemberSnapshot(m.exposure_id, qc)
        : (m.snapshot ?? {
            effective_peaks: [],
            confirmed_index: null,
            analysis_inputs_hash: "",
          });
      const out: ComparisonMemberInput = {
        exposure_id: m.exposure_id,
        display_order: m.display_order,
        band_height: m.band_height,
        y_offset: m.y_offset,
        normalization: m.normalization,
        snapshot,
      };
      if (m.id !== undefined) out.id = m.id;
      if (m.color_override !== undefined) out.color_override = m.color_override;
      if (m.label_override !== undefined) out.label_override = m.label_override;
      if (m.q_window_min !== undefined) out.q_window_min = m.q_window_min;
      if (m.q_window_max !== undefined) out.q_window_max = m.q_window_max;
      if (m.peak_display !== undefined) out.peak_display = m.peak_display;
      return out;
    });
    // useSaveComparison flat-spreads the input into the SaveComparisonBody;
    // see saveComparison mutator's `request: (p) => api.saveComparison(...)`.
    const payload: SaveComparisonBody & { id?: number } = {
      title: draft.title,
      members,
    };
    if (draft.id !== undefined) payload.id = draft.id;
    if (draft.description !== "") payload.description = draft.description;
    if (draft.baseHash !== undefined) payload.expected_content_hash = draft.baseHash;
    // Mark this submit as "in-flight" so the post-success effect knows to
    // navigate. Without the ref, an already-saved success state on remount
    // would re-fire navigation.
    pendingSubmitRef.current = true;
    save.mutate(payload);
  }, [draft, qc, save]);

  // Post-success navigation. Reading `save.data` (the response) lets us
  // navigate to /experiments/:eid/compare/:newId for create flow. We gate
  // on `pendingSubmitRef` so only a save we just initiated triggers nav.
  useEffect(() => {
    if (!save.isSuccess) return;
    if (!pendingSubmitRef.current) return;
    pendingSubmitRef.current = false;
    const response = save.data as Comparison | undefined;
    discardDraft();
    if (response && typeof response.id === "number") {
      goToReview(response.id);
    } else {
      goToList();
    }
  }, [save.isSuccess, save.data, discardDraft, goToReview, goToList]);

  const scope: "all" | "experiment" = eid !== undefined ? "experiment" : "all";

  // Phase 5 Task 5.2 — local state for the picker open/close. Picker open
  // state is purely client-side, no need for Zustand.
  const [pickerOpen, setPickerOpen] = useState(false);
  // Phase 11 wires Edit/Fork visibility against current_user vs. created_by;
  // for now we surface the testid so downstream tests can target it once
  // the gating exists. The button is hidden from the rendered tree until
  // then to keep the edit shell clean.

  return (
    <div
      data-testid="compare-page-edit"
      {...(id !== undefined ? { "data-comparison-id": String(id) } : {})}
      className="flex-1 min-h-0 flex gap-3 px-4 pb-4 pt-2"
    >
      <aside className="card overflow-hidden w-[300px] shrink-0 flex flex-col">
        <ComparisonSidebar
          experimentId={eid}
          scope={scope}
          activeComparisonId={id}
        />
      </aside>
      <section className="card overflow-hidden flex-1 min-h-0 flex flex-col p-4 gap-3">
        <div className="flex items-center gap-2">
          <input
            data-testid="compare-edit-title"
            type="text"
            placeholder="Comparison title"
            value={draft?.title ?? ""}
            onChange={(e) => setDraftTitle(e.target.value)}
            className="flex-1 bg-transparent border border-border rounded px-2 py-1 text-base"
          />
          <button
            type="button"
            data-testid="comparison-save"
            onClick={handleSave}
            disabled={(draft?.members.length ?? 0) === 0 || save.isPending}
            className="px-3 py-1 rounded bg-accent text-bg disabled:opacity-50 text-sm font-medium"
          >
            {save.isPending ? "Saving…" : "Save"}
          </button>
          <button
            type="button"
            data-testid="comparison-cancel"
            onClick={handleCancel}
            className="px-3 py-1 rounded border border-border text-sm"
          >
            Cancel
          </button>
          <button
            type="button"
            data-testid="comparison-discard"
            onClick={handleDiscard}
            className="px-3 py-1 rounded border border-border text-fg-muted text-sm"
          >
            Discard draft
          </button>
        </div>
        <textarea
          data-testid="compare-edit-description"
          placeholder="Description (optional)"
          value={draft?.description ?? ""}
          onChange={(e) => setDraftDescription(e.target.value)}
          className="bg-transparent border border-border rounded px-2 py-1 text-sm resize-none h-16"
        />
        <div
          data-testid="compare-edit-plot-placeholder"
          className="flex-1 min-h-0 flex flex-col items-center justify-center
                     border border-border/40 rounded text-fg-muted text-sm gap-3"
        >
          <div>
            Plot + member panel land in Phases 6–9 ({(draft?.members.length ?? 0)} member{(draft?.members.length ?? 0) === 1 ? "" : "s"})
          </div>
          <button
            type="button"
            data-testid="compare-edit-add-traces"
            onClick={() => setPickerOpen(true)}
            className="px-3 py-1 rounded border border-border text-fg text-sm
                       hover:bg-bg-elevated"
          >
            + Add traces
          </button>
        </div>
      </section>

      <ComparisonPicker
        isOpen={pickerOpen}
        onClose={() => setPickerOpen(false)}
        experimentId={eid}
      />
    </div>
  );
}
