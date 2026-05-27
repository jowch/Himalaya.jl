import { useCallback } from "react";
import { useExposures, useSetExposureStatus } from "../queries";
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
  onToggleReject,
}: {
  exposure: Exposure;
  onToggleReject: (exp: Exposure) => void;
}): JSX.Element {
  const isRejected = exposure.status === "rejected";
  const isRepresentative = exposure.selected;
  return (
    <div
      data-testid={`exposure-thumb-${exposure.id}`}
      data-rejected={isRejected ? "true" : undefined}
      data-representative={isRepresentative ? "true" : undefined}
      className={[
        "relative w-12 shrink-0 aspect-[3/4] overflow-hidden rounded",
        "ring-1 ring-hair",
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
 * Inert affordances: the thumbnails carry no onClick and the tag-add
 * button is disabled — selection (#162) and tag mutation (#159) wire in
 * separately.
 */
export function ContactSheetRow({ sample }: Props): JSX.Element {
  const exposuresQuery = useExposures(sample.id);
  const exposures = exposuresQuery.data ?? [];

  const setStatus = useSetExposureStatus(sample.id);

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

      {/* Exposures — thumbnail strip. */}
      <div
        data-testid="exposures-cell"
        className="flex h-16 flex-row gap-2 overflow-x-auto"
      >
        {exposuresQuery.isLoading ? (
          <span className="self-center text-xs text-ink-faint">
            Loading frames…
          </span>
        ) : (
          exposures.map((e) => (
            <ExposureThumb
              key={e.id}
              exposure={e}
              onToggleReject={handleToggleReject}
            />
          ))
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
