import { useCallback, useState } from "react";
import {
  useExposures,
  useSetExposureStatus,
  useSelectExposure,
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
}

/**
 * One exposure thumbnail with culling affordances (#162) in the contact-sheet
 * idiom (R1): a square 62px dark window (L-4), a zero-padded frame badge in the
 * corner (M-7), and the grease-pencil ✕ over rejected frames (M-10).
 *
 * `frameNo` is the 1-based position of the exposure in its sample's strip.
 */
function ExposureThumb({
  exposure,
  frameNo,
  selectedForBatch,
  onToggleReject,
  onToggleSelected,
  onPickRepresentative,
}: {
  exposure: Exposure;
  frameNo: number;
  selectedForBatch: boolean;
  onToggleReject: (exp: Exposure) => void;
  onToggleSelected: (id: number) => void;
  onPickRepresentative: (exp: Exposure) => void;
}): JSX.Element {
  const isRejected = exposure.status === "rejected";
  const isRepresentative = exposure.selected;
  return (
    <div
      data-testid={`exposure-thumb-${exposure.id}`}
      data-rejected={isRejected ? "true" : undefined}
      data-representative={isRepresentative ? "true" : undefined}
      data-batch-selected={selectedForBatch ? "true" : undefined}
      className={[
        // square dark window (62px) — L-4
        "relative h-[62px] w-[62px] shrink-0 overflow-hidden rounded-[3px]",
        "border border-frame-edge bg-frame-edge",
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
      <button
        type="button"
        data-testid={`exposure-select-${exposure.id}`}
        aria-pressed={selectedForBatch}
        title="Select for batch action"
        onClick={() => onToggleSelected(exposure.id)}
        className={[
          "absolute left-0 top-0 m-0.5 h-3 w-3 rounded-sm border",
          selectedForBatch
            ? "border-print-accent bg-print-accent"
            : "border-hair-strong bg-paper/80",
        ].join(" ")}
      />
      {/* representative pick — the frame that goes on to indexing */}
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
      {!isRepresentative && (
        <button
          type="button"
          data-testid={`exposure-represent-${exposure.id}`}
          title="Make representative"
          onClick={() => onPickRepresentative(exposure)}
          className="absolute bottom-0 left-0 m-0.5 rounded bg-paper/80 px-1
                     text-[10px] leading-none text-ink-faint"
        >
          ⊙
        </button>
      )}
      <button
        type="button"
        data-testid={`exposure-reject-${exposure.id}`}
        title={isRejected ? "Un-reject exposure" : "Reject exposure"}
        onClick={() => onToggleReject(exposure)}
        className="absolute bottom-0 right-0 m-0.5 rounded bg-paper/80 px-1
                   text-[10px] leading-none text-print-accent"
      >
        {isRejected ? "↺" : "✕"}
      </button>
    </div>
  );
}

/**
 * ContactSheetRow — one sample row of the corpus contact sheet (#160).
 *
 * Owns its own useExposures query (per-sample fan-out) so the table fills
 * in row-by-row. The same queryKeys.exposures(sampleId) cache entry is
 * reused by culling (#162) and the loupe (#161).
 *
 * Culling is wired here (#162): per-thumb reject toggle, multi-select batch
 * reject, and representative pick, all through the existing exposure queue
 * hooks. The tag-add button remains inert — sample-tag mutation is #159.
 */
export function ContactSheetRow({ sample }: Props): JSX.Element {
  const exposuresQuery = useExposures(sample.id);
  const exposures = exposuresQuery.data ?? [];

  const setStatus = useSetExposureStatus(sample.id);
  const setRepresentative = useSelectExposure(sample.id);

  // Representative pick. selectExposureMutator's onMutate writes
  // `selected: e.id === exposureId` across the list, so the pick is
  // mutually exclusive (one representative per sample) for free.
  const handlePickRepresentative = useCallback(
    (exp: Exposure) => {
      setRepresentative.mutate(exp.id);
    },
    [setRepresentative],
  );

  // Single-exposure reject toggle. Un-reject sets status to null (matches
  // LoupePage's `status === "rejected" ? null : "rejected"` convention).
  const handleToggleReject = useCallback(
    (exp: Exposure) => {
      setStatus.mutate({
        exposureId: exp.id,
        status: exp.status === "rejected" ? null : "rejected",
      });
    },
    [setStatus],
  );

  // Multi-select state — local to this row (selection never crosses samples,
  // matching the per-sample query fan-out). Stale ids are harmless: the batch
  // handler filters to currently-present, non-rejected exposures.
  const [selectedIds, setSelectedIds] = useState<ReadonlySet<number>>(
    () => new Set(),
  );
  const toggleSelected = useCallback((id: number) => {
    setSelectedIds((prev) => {
      const next = new Set(prev);
      if (next.has(id)) next.delete(id);
      else next.add(id);
      return next;
    });
  }, []);
  const clearSelection = useCallback(() => setSelectedIds(new Set()), []);

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

  const total = exposures.length;
  const kept = exposures.filter((e) => e.status !== "rejected").length;
  const dropped = total - kept;

  // M-2: screened state. Derived from the exposures today; flips to #162's
  // backend flag automatically once it lands (see lib/sample/screened.ts).
  const screened = isSampleScreened(sample, exposures);

  // Route through the shared helper — `||` semantics (empty-string-safe)
  // and one source of truth for "what string do we render for a sample".
  const name = sampleDisplayName(sample);

  return (
    <div
      data-testid={`sample-row-${sample.id}`}
      data-unscreened={screened ? undefined : "true"}
      className={`${CONTACT_SHEET_COLS} border-b border-hair px-4
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
          {exposuresQuery.isLoading ? (
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
                onToggleReject={handleToggleReject}
                onToggleSelected={toggleSelected}
                onPickRepresentative={handlePickRepresentative}
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
        {exposuresQuery.isLoading ? (
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

      {/* Tags — read-only chips + inert add button (mutation is #159). */}
      <div
        data-testid="tags-cell"
        className="flex min-h-[92px] flex-wrap items-center gap-1"
      >
        {sample.tags.map((t) => (
          <span
            key={t.id}
            title={t.key || undefined}
            className="rounded bg-paper-sunk px-1.5 py-0.5 text-xs text-ink-soft"
          >
            {t.value}
          </span>
        ))}
        <button
          type="button"
          data-testid="tag-add"
          disabled
          title="Add a tag (coming soon)"
          className="rounded border border-dashed border-hair-strong px-1.5
                     py-0.5 text-xs text-ink-faint"
        >
          + tag
        </button>
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
