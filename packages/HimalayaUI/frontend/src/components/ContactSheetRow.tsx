import { useCallback, useState } from "react";
import {
  useExposures,
  useSetExposureStatus,
  useSelectExposure,
} from "../queries";
import type { CorpusSample, Exposure } from "../api";
import { sampleDisplayName } from "../lib/sample/displayName";
import { DetectorImage } from "./DetectorImage";

/**
 * A sample's index/phase status column. #160 ships only "not-indexed";
 * a later issue wires the real phase call into this typed seam.
 */
export type SampleStatus = "not-indexed";

/**
 * Shared CSS grid template for the contact sheet — the column header in
 * SamplesPage and every ContactSheetRow use this so the columns align.
 */
export const CONTACT_SHEET_COLS =
  "grid grid-cols-[16rem_1fr_7rem_14rem_8rem] gap-4 items-center";

interface Props {
  sample: CorpusSample;
}

/** One exposure thumbnail with culling affordances (#162). */
function ExposureThumb({
  exposure,
  selectedForBatch,
  onToggleReject,
  onToggleSelected,
  onPickRepresentative,
}: {
  exposure: Exposure;
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
        "relative w-12 shrink-0 aspect-[3/4] overflow-hidden rounded",
        selectedForBatch ? "ring-2 ring-accent" : "ring-1 ring-hair",
        isRejected ? "opacity-40 grayscale" : "",
      ].join(" ")}
    >
      <DetectorImage
        exposureId={exposure.id}
        imagePath={exposure.image_path}
        imageVersion={exposure.image_version}
        size="thumb"
        className="h-full w-full"
      />
      <button
        type="button"
        data-testid={`exposure-select-${exposure.id}`}
        aria-pressed={selectedForBatch}
        title="Select for batch action"
        onClick={() => onToggleSelected(exposure.id)}
        className={[
          "absolute left-0 top-0 m-0.5 h-3 w-3 rounded-sm border",
          selectedForBatch
            ? "border-accent bg-accent"
            : "border-hair-strong bg-paper/80",
        ].join(" ")}
      />
      {isRepresentative && (
        <span
          className="absolute left-0.5 top-0.5 text-[10px] text-print-accent"
          title="representative exposure"
        >
          ⊙
        </span>
      )}
      {isRejected && (
        <span className="absolute inset-0 flex items-center justify-center text-print-accent">
          ✕
        </span>
      )}
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

  // Route through the shared helper — `||` semantics (empty-string-safe)
  // and one source of truth for "what string do we render for a sample".
  const name = sampleDisplayName(sample);

  return (
    <div
      data-testid={`sample-row-${sample.id}`}
      className={`${CONTACT_SHEET_COLS} border-b border-hair px-4 py-3`}
    >
      {/* Sample — identity only (no screened mark; that is #162). */}
      <div data-testid="sample-cell" className="flex flex-col">
        <span className="font-semibold text-ink">{name}</span>
        <span className="text-xs text-ink-faint">#{sample.id}</span>
      </div>

      {/* Exposures — thumbnail strip + (when a selection exists) an action
          bar. The strip is a fixed-height horizontal scroller; the action bar
          is a flex-col SIBLING below it, never a clipped child of the
          overflow-x-auto strip. */}
      <div data-testid="exposures-cell" className="flex flex-col gap-1.5">
        <div className="flex h-16 flex-row gap-2 overflow-x-auto">
          {exposuresQuery.isLoading ? (
            <span className="self-center text-xs text-ink-faint">
              Loading frames…
            </span>
          ) : (
            exposures.map((e) => (
              <ExposureThumb
                key={e.id}
                exposure={e}
                selectedForBatch={selectedIds.has(e.id)}
                onToggleReject={handleToggleReject}
                onToggleSelected={toggleSelected}
                onPickRepresentative={handlePickRepresentative}
              />
            ))
          )}
        </div>
        {selectedIds.size > 0 && (
          <div
            data-testid="contact-sheet-actionbar"
            className="flex items-center gap-2"
          >
            <button
              type="button"
              data-testid="batch-reject"
              onClick={handleBatchReject}
              className="rounded border border-print-accent px-1.5 py-0.5
                         text-xs text-print-accent"
            >
              Reject {selectedIds.size} selected
            </button>
            <button
              type="button"
              data-testid="batch-clear"
              onClick={clearSelection}
              className="text-xs text-ink-faint hover:underline"
            >
              Clear
            </button>
          </div>
        )}
      </div>

      {/* Kept — kept / total, plus an "N dropped" sub-label. */}
      <div data-testid="kept-cell" className="flex flex-col text-sm">
        {exposuresQuery.isLoading ? (
          <span className="text-ink-faint">—</span>
        ) : (
          <>
            <span className="text-ink">
              {kept}
              <span className="text-ink-faint"> / {total}</span>
            </span>
            {dropped > 0 && (
              <span className="text-xs text-print-accent">
                {dropped} dropped
              </span>
            )}
          </>
        )}
      </div>

      {/* Tags — read-only chips + inert add button (mutation is #159). */}
      <div data-testid="tags-cell" className="flex flex-wrap items-center gap-1">
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

      {/* Status — fixed placeholder behind the SampleStatus seam.
          TODO: wire the real phase call when an issue is scoped for it. */}
      <div data-testid="status-cell">
        <span className="text-xs text-ink-faint">Not indexed</span>
      </div>
    </div>
  );
}
