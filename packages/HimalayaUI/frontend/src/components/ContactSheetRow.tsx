import { useCallback, useEffect, useRef, useState } from "react";
import { useNavigate } from "react-router-dom";
import {
  useExposures,
  useSetExposureStatus,
  useAddCorpusSampleTag,
  useRemoveCorpusSampleTag,
} from "../queries";
import type { CorpusSample, Exposure } from "../api";
import { sampleDisplayName } from "../lib/sample/displayName";
import { isSampleScreened } from "../lib/sample/screened";
import { DetectorImage } from "./DetectorImage";
import { RejectXMark } from "./RejectXMark";
import { SampleStatusChip } from "./SampleStatusChip";
import { CullBar } from "./CullBar";

/**
 * A sample's index/phase status column. #160 ships only "not-indexed";
 * a later issue wires the real phase call into this typed seam.
 */
export type SampleStatus = "not-indexed";

/**
 * Shared CSS grid template for the contact sheet — the column header in
 * SamplesPage and every ContactSheetRow use this so the columns align.
 *
 * Column widths track sample-table.html's `.COLS` (L-5): a fixed specimen
 * column, a flexible exposure strip, then narrow kept / tags / status
 * columns. `items-stretch` so the screened mark can top-align in a tall row.
 */
export const CONTACT_SHEET_COLS =
  "grid grid-cols-[15.25rem_minmax(22.5rem,1fr)_4.875rem_10.5rem_9.375rem] gap-4 items-stretch";

interface Props {
  sample: CorpusSample;
  /**
   * Pre-fetched exposures from the parent's bulk hook (R2-M14). Optional —
   * when omitted, the row falls back to its own `useExposures(sample.id)`
   * subscription (the legacy fan-out path) so unit tests and any non-bulk
   * caller keep working. The bulk path's `undefined` value distinguishes
   * "not yet fetched" (the loading skeleton) from "fetched, empty" ([]).
   */
  exposures?: Exposure[] | undefined;
  exposuresLoading?: boolean;
}

/**
 * One exposure thumbnail in the contact-sheet idiom (R1): a square 62px dark
 * window (L-4), a zero-padded frame badge in the corner (M-7), and the
 * grease-pencil ✕ over rejected frames (M-10).
 *
 * R2-M11: the three permanent overlay buttons (checkbox / rep ⊙ / reject ✕)
 * are gone — the existing `ring-print-accent` selection ring, the floating
 * CullBar, and the `R` / `X` keystrokes already cover the function. The rep
 * dot remains as a corner state badge on `is-rep` (mockup `.thumb.is-rep .rep`).
 *
 * `frameNo` is the 1-based position of the exposure in its sample's strip.
 */
function ExposureThumb({
  exposure,
  frameNo,
  selectedForBatch,
  onSelect,
  onOpenLoupe,
}: {
  exposure: Exposure;
  frameNo: number;
  selectedForBatch: boolean;
  /** Body click. `extend` is the shift-click range modifier. */
  onSelect: (id: number, extend: boolean) => void;
  onOpenLoupe: () => void;
}): JSX.Element {
  const isRejected = exposure.status === "rejected";
  const isRepresentative = exposure.selected;
  return (
    <button
      type="button"
      data-testid={`exposure-thumb-${exposure.id}`}
      data-rejected={isRejected ? "true" : undefined}
      data-representative={isRepresentative ? "true" : undefined}
      data-batch-selected={selectedForBatch ? "true" : undefined}
      // R3-S01 (#256, P1 a11y): the thumb root is a real <button> so the
      // contact sheet is keyboard-operable. R2-M11 stripped the per-thumb
      // overlay buttons — the thumb's only keyboard-reachable focus targets —
      // so the window body itself becomes the focusable control.
      aria-label={`Frame ${frameNo}${isRejected ? " (dropped)" : ""}`}
      // L-5 legend "click — select a frame" / "⇧ click — extend the range":
      // the whole thumb body is the click target (shiftKey extends the range).
      // double-click opens the loupe (L-8). With R2-M11 the per-thumb buttons
      // are retired, so no inner controls need stopPropagation.
      onClick={(e) => onSelect(exposure.id, e.shiftKey)}
      onDoubleClick={onOpenLoupe}
      // Keyboard Enter/Space → open the loupe (the navigation affordance an AT
      // user needs to screen a sample; the keyboard analogue of double-click,
      // since there is no single-key analogue of "double-click"). preventDefault
      // suppresses the <button>'s synthesized activation click so onSelect does
      // NOT also fire — one key, one action.
      onKeyDown={(e) => {
        if (e.key === "Enter" || e.key === " ") {
          e.preventDefault();
          onOpenLoupe();
        }
      }}
      className={[
        // square dark window (62px) — L-4. appearance-none + p-0 neutralise UA
        // button chrome; the box is fully specified by the size/border classes.
        "relative h-[62px] w-[62px] shrink-0 cursor-pointer appearance-none p-0",
        "overflow-hidden rounded-[3px] border border-frame-edge bg-frame-edge",
        "focus-visible:ring-2 focus-visible:ring-print-accent focus-visible:ring-offset-1",
        selectedForBatch ? "ring-2 ring-print-accent" : "ring-0",
      ].join(" ")}
    >
      <DetectorImage
        exposureId={exposure.id}
        imagePath={exposure.image_path}
        imageVersion={exposure.image_version}
        size="thumb"
        className={`h-full w-full transition-opacity ${
          isRejected ? "opacity-30" : ""
        }`}
      />
      {/* M-7: zero-padded frame-number badge */}
      <span
        data-testid={`frame-no-${exposure.id}`}
        className={`pointer-events-none absolute bottom-px left-[3px] font-mono
                    text-[8.5px] text-paper/80 ${isRejected ? "opacity-45" : ""}`}
      >
        {String(frameNo).padStart(2, "0")}
      </span>
      {/* representative pick — the frame that goes on to indexing (mockup
          `.thumb.is-rep .rep`). Rendered as a pure state badge; the rep pick
          itself happens in the loupe or by the implicit default-selection. */}
      {isRepresentative && (
        <span
          data-testid={`exposure-rep-dot-${exposure.id}`}
          className="pointer-events-none absolute right-[3px] top-[3px] h-[9px] w-[9px]
                     rounded-full border-[1.5px] border-plate bg-print-accent"
          title="representative exposure"
        />
      )}
      {/* M-10: the grease-pencil reject mark */}
      {isRejected && <RejectXMark />}
    </button>
  );
}

/**
 * ContactSheetRow — one sample row of the corpus contact sheet (#160).
 *
 * R2-M14 (#207): when the parent passes `exposures` from the bulk
 * `useCorpusExposures()` hook, no per-row query observer is mounted — that
 * was the 139× JSON fan-out behind the `ERR_INSUFFICIENT_RESOURCES` on the
 * live surface. The row still owns its own observer when called without the
 * prop (the historical path, kept for backward-compatible callers + tests).
 *
 * Culling is wired here (#162): per-thumb selection via click / shift-click,
 * batch-reject through the floating CullBar + `X` keystroke. The rep pick
 * lives in the loupe (R2-M11 stripped the thumb-overlay rep ⊙ button).
 *
 * Corpus tag add/delete is wired here (#159 / #207): the `+ tag` chip opens
 * a tiny inline form on the row; existing chips carry a remove ✕. Mutations
 * route through `useAddCorpusSampleTag` / `useRemoveCorpusSampleTag`, so the
 * loupe and the sheet share one cache row.
 *
 * The footer-legend affordances are all live here: click / shift-click
 * select & extend a contiguous range over this sample's frames, double-click
 * opens the loupe, and while a selection exists `X` batch-rejects it and
 * `Esc` clears it. The keydown listener is row-scoped and only mounts when
 * this row owns the live selection, so the per-sample model holds (a
 * keyboard action never touches another row's frames).
 */
export function ContactSheetRow({
  sample,
  exposures: bulkExposures,
  exposuresLoading: bulkLoading,
}: Props): JSX.Element {
  const navigate = useNavigate();

  // R2-M14: bulk-fed rows skip the per-row observer entirely. The `enabled`
  // gate flips false when bulk data is in, so the legacy `useExposures` hook
  // doesn't mount a duplicate subscription. The cache row itself is shared
  // (same queryKey), so a bulk-fed row still picks up mutator-side patches.
  const useBulk = bulkExposures !== undefined;
  const exposuresQuery = useExposures(useBulk ? undefined : sample.id);
  const exposures: Exposure[] = useBulk
    ? bulkExposures!
    : (exposuresQuery.data ?? []);
  const exposuresLoading = useBulk
    ? (bulkLoading ?? false)
    : exposuresQuery.isLoading;

  const setStatus = useSetExposureStatus(sample.id);
  const addTag = useAddCorpusSampleTag(sample.id);
  const removeTag = useRemoveCorpusSampleTag(sample.id);

  // Multi-select state — local to this row (selection never crosses samples,
  // matching the per-sample query fan-out). Stale ids are harmless: the batch
  // handler filters to currently-present, non-rejected exposures.
  const [selectedIds, setSelectedIds] = useState<ReadonlySet<number>>(
    () => new Set(),
  );
  // Range-select anchor: the last single-clicked frame id. A shift-click
  // selects the contiguous span between the anchor and the clicked frame.
  const anchorIdRef = useRef<number | null>(null);
  // Always-current frame order for range math (avoids stale closures when the
  // exposures query refetches). Frame order == the rendered exposures array.
  const orderRef = useRef<number[]>([]);
  orderRef.current = exposures.map((e) => e.id);

  const clearSelection = useCallback(() => setSelectedIds(new Set()), []);

  // Select a frame. `extend` (shift-click) selects the contiguous range from
  // the anchor to this frame; a plain click toggles the single frame and
  // re-anchors. Range/anchor stay within this row's frames by construction.
  const handleSelect = useCallback((id: number, extend: boolean) => {
    setSelectedIds((prev) => {
      const order = orderRef.current;
      if (extend && anchorIdRef.current !== null) {
        const a = order.indexOf(anchorIdRef.current);
        const b = order.indexOf(id);
        if (a !== -1 && b !== -1) {
          const [lo, hi] = a <= b ? [a, b] : [b, a];
          const next = new Set(prev);
          for (let i = lo; i <= hi; i++) next.add(order[i]);
          return next;
        }
      }
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      anchorIdRef.current = id;
      return next;
    });
  }, []);

  const openLoupe = useCallback(() => {
    navigate(`/samples/loupe/${sample.id}`);
  }, [navigate, sample.id]);

  // Batch reject — N independent ops, one per currently-kept selected
  // exposure. Each setStatus.mutate() mints its own client_op_id and applies
  // its own optimistic patch to the shared exposures cache, so the patches
  // compose without a batch mutator.
  const handleBatchReject = useCallback(() => {
    for (const exp of exposures) {
      if (selectedIds.has(exp.id) && exp.status !== "rejected") {
        setStatus.mutate({ exposureId: exp.id, status: "rejected" });
      }
    }
    clearSelection();
  }, [exposures, selectedIds, setStatus, clearSelection]);

  // Keyboard affordances (the CullBar / footer-legend `X` and `Esc` keycaps).
  // Only mounted while this row owns a live selection, so it never competes
  // with another row and never fires when nothing is selected. Suppressed
  // while typing in a field, mirroring useGlobalShortcuts.
  const hasSelection = selectedIds.size > 0;
  useEffect(() => {
    if (!hasSelection) return;
    const onKeyDown = (e: KeyboardEvent): void => {
      const t = e.target as HTMLElement | null;
      const editing =
        t &&
        (t.tagName === "INPUT" ||
          t.tagName === "TEXTAREA" ||
          t.isContentEditable);
      if (editing || e.metaKey || e.ctrlKey || e.altKey) return;
      if (e.key === "Escape") {
        e.preventDefault();
        clearSelection();
      } else if (e.key === "x" || e.key === "X") {
        e.preventDefault();
        handleBatchReject();
      }
    };
    window.addEventListener("keydown", onKeyDown);
    return () => window.removeEventListener("keydown", onKeyDown);
  }, [hasSelection, clearSelection, handleBatchReject]);

  const total = exposures.length;
  const kept = exposures.filter((e) => e.status !== "rejected").length;
  const dropped = total - kept;

  // M-2: screened state. Derived from the exposures today; flips to #162's
  // backend flag automatically once it lands (see lib/sample/screened.ts).
  const screened = isSampleScreened(sample, exposures);

  // Route through the shared helper — `||` semantics (empty-string-safe)
  // and one source of truth for "what string do we render for a sample".
  const name = sampleDisplayName(sample);

  // Tag add/delete UI state — the row keeps its own inline form when the user
  // taps `+`. The submit routes through useAddCorpusSampleTag (which is the
  // same mutator the loupe uses, so an add here is visible there immediately).
  const [tagFormOpen, setTagFormOpen] = useState(false);
  const [tagKeyDraft, setTagKeyDraft] = useState("");
  const [tagValDraft, setTagValDraft] = useState("");
  const resetTagForm = useCallback(() => {
    setTagKeyDraft("");
    setTagValDraft("");
    setTagFormOpen(false);
  }, []);
  const submitTag = useCallback(() => {
    const v = tagValDraft.trim();
    if (v === "") return;
    addTag.mutate({ key: tagKeyDraft.trim(), value: v });
    resetTagForm();
  }, [addTag, tagKeyDraft, tagValDraft, resetTagForm]);

  return (
    <div
      data-testid={`sample-row-${sample.id}`}
      data-unscreened={screened ? undefined : "true"}
      className={`group ${CONTACT_SHEET_COLS} border-b border-hair px-4
                  ${screened ? "" : "bg-paper-sunk/60"}`}
    >
      {/* Sample — screened mark (M-2) + identity. The mark top-aligns in the
          tall row, echoing the grease-pencil tick on a screened frame. */}
      <div
        data-testid="sample-cell"
        className="flex min-h-[92px] items-center gap-[11px] py-[13px]"
      >
        <span
          data-testid="screened-mark"
          data-screened={screened ? "true" : undefined}
          title={screened ? "screened" : "not yet screened"}
          className={[
            "flex h-[13px] w-[13px] shrink-0 items-center justify-center rounded-full border-[1.5px]",
            screened
              ? "border-ink bg-ink"
              : "border-hair-strong bg-transparent",
          ].join(" ")}
        >
          {screened && (
            <svg viewBox="0 0 13 13" className="h-full w-full" aria-hidden="true">
              <path
                d="M3.4 6.8l2.1 2.1 4.2-4.6"
                fill="none"
                stroke="var(--color-paper)"
                strokeWidth="1.7"
                strokeLinecap="round"
                strokeLinejoin="round"
              />
            </svg>
          )}
        </span>
        <span className="min-w-0">
          <span className="block truncate text-[13.5px] font-semibold text-ink">
            {name}
          </span>
          <span className="mt-0.5 block font-mono text-[10.5px] text-ink-faint">
            #{sample.id}
          </span>
        </span>
      </div>

      {/* Exposures — square thumbnail strip. The floating cull bar (M-3)
          replaces the inline action bar; it is rendered at the row root so it
          escapes the overflow-x-auto strip's clipping. */}
      <div
        data-testid="exposures-cell"
        className="flex min-h-[92px] items-center py-[13px]"
      >
        <div className="flex flex-row flex-nowrap gap-[7px] overflow-x-auto">
          {exposuresLoading ? (
            <span className="self-center text-xs text-ink-faint">
              Loading frames…
            </span>
          ) : (
            exposures.map((e, i) => (
              <ExposureThumb
                key={e.id}
                exposure={e}
                frameNo={i + 1}
                selectedForBatch={selectedIds.has(e.id)}
                onSelect={handleSelect}
                onOpenLoupe={openLoupe}
              />
            ))
          )}
        </div>
      </div>

      {/* Kept — kept / total, plus an "N dropped" sub-label. */}
      <div
        data-testid="kept-cell"
        className="flex min-h-[92px] flex-col justify-center font-mono text-sm"
      >
        {exposuresLoading ? (
          <span className="text-ink-faint">—</span>
        ) : (
          <>
            <span className="text-ink">
              <span className="text-base">{kept}</span>
              <span className="text-ink-faint"> / {total}</span>
            </span>
            {dropped > 0 && (
              <span className="font-sans text-[10px] font-semibold text-print-accent">
                {dropped} dropped
              </span>
            )}
          </>
        )}
      </div>

      {/* Tags — chips with inline remove + `+ tag` invite that opens a tiny
          add form. The `+ tag` invite shows always when the row is empty;
          when chips exist, the `+` is hover-revealed (mockup `.tags .tag-add
          :not(.invite)` opacity rule). #159 corpus tag mutators are wired. */}
      <div
        data-testid="tags-cell"
        className="flex min-h-[92px] flex-wrap items-center gap-1"
      >
        {sample.tags.map((t) => (
          <span
            key={t.id}
            title={t.key || undefined}
            data-testid={`sample-tag-${t.id}`}
            className="inline-flex items-center gap-1 rounded-full border border-hair
                       bg-plate px-2 py-0.5 text-[10.5px] font-semibold text-ink-soft"
          >
            {t.value}
            <button
              type="button"
              aria-label={`Remove ${t.key || t.value} tag`}
              onClick={() => removeTag.mutate(t.id)}
              className="text-ink-faint hover:text-print-accent"
            >
              ×
            </button>
          </span>
        ))}
        {tagFormOpen ? (
          <span
            data-testid="tag-form"
            className="inline-flex items-center gap-1 rounded-full border border-hair-strong
                       bg-plate px-1.5 py-0.5"
          >
            <input
              aria-label="tag key"
              placeholder="key"
              value={tagKeyDraft}
              onChange={(e) => setTagKeyDraft(e.target.value)}
              className="w-12 bg-transparent text-[10.5px] text-ink outline-none
                         placeholder:text-ink-faint"
            />
            <span className="text-ink-faint">:</span>
            <input
              aria-label="tag value"
              placeholder="value"
              value={tagValDraft}
              onChange={(e) => setTagValDraft(e.target.value)}
              onKeyDown={(e) => {
                if (e.key === "Enter") submitTag();
                else if (e.key === "Escape") resetTagForm();
              }}
              autoFocus
              className="w-16 bg-transparent text-[10.5px] text-ink outline-none
                         placeholder:text-ink-faint"
            />
            <button
              type="button"
              onClick={submitTag}
              className="text-[10.5px] font-semibold text-print-accent
                         hover:underline"
            >
              Add
            </button>
            <button
              type="button"
              aria-label="cancel adding tag"
              onClick={resetTagForm}
              className="text-ink-faint hover:text-print-accent"
            >
              ×
            </button>
          </span>
        ) : (
          <button
            type="button"
            data-testid="tag-add"
            onClick={() => setTagFormOpen(true)}
            title="Add a tag"
            className={[
              "rounded-full border border-dashed border-hair-strong px-2 py-0.5",
              "text-[10.5px] font-semibold text-ink-faint",
              "hover:border-print-accent hover:text-print-accent",
              // Hover-revealed when chips exist; always-visible invite when empty.
              sample.tags.length === 0
                ? ""
                : "opacity-0 group-hover:opacity-100 focus:opacity-100",
            ].join(" ")}
          >
            {sample.tags.length === 0 ? "+ tag" : "+"}
          </button>
        )}
      </div>

      {/* Status — phase chip when a phase is present (M-6), else the hollow-dot
          "Not indexed" affordance. `sample.phase` is a forward-looking seam:
          no corpus-level indexing rollup is wired yet, so it is absent today
          and the cell reads "Not indexed". */}
      <div
        data-testid="status-cell"
        className="flex min-h-[92px] items-center"
      >
        <SampleStatusChip phase={sample.phase} />
      </div>

      {/* M-3: floating cull bar — rendered from the row root (outside the
          grid cell + overflow strip) so it floats at the bottom centre of the
          viewport. Selection stays per-sample by design. */}
      {selectedIds.size > 0 && (
        <CullBar
          count={selectedIds.size}
          onReject={handleBatchReject}
          onClear={clearSelection}
        />
      )}
    </div>
  );
}
