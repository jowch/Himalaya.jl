/**
 * Ingest progress stages — the segments of the scan progress bar.
 *
 * `scan_and_group!` (ingest.jl) reports `(processed, total, stage)`. Each stage
 * is a SEGMENT with its OWN denominator, so the bar advances left-to-right
 * through segments instead of resetting a single shared scale between phases.
 *
 * Segments are EQUAL WIDTH on purpose. Weighting them by expected duration would
 * mean inventing per-phase time estimates that drift with corpus shape (a rescan
 * with few new files spends almost no time analyzing; a cold first scan spends
 * most of it there) — that is how progress bars end up stalling at 40% or
 * jumping. Equal width claims only "these are the stages, here is the one you're
 * in", which is true regardless of corpus.
 *
 * NOTE this is deliberately NOT keyed off `IngestProgress.status`, which selects
 * the SURFACE (`"scanning"` → the grouping takeover, `"analyzing"` → the inline
 * bar on the corpus page) and is driven by the route's `phase` field. Conflating
 * the two would flip the whole surface mid-scan.
 */

/** Wire values emitted by `on_progress`, in execution order. */
export const INGEST_STAGES = ["discovery", "analyzing", "thumbnails"] as const;

export type IngestStage = (typeof INGEST_STAGES)[number];

const STAGE_LABELS: Record<IngestStage, string> = {
  discovery: "Reading files",
  analyzing: "Analyzing",
  thumbnails: "Thumbnails",
};

export function isIngestStage(v: string | undefined): v is IngestStage {
  return v !== undefined && (INGEST_STAGES as readonly string[]).includes(v);
}

export function stageLabel(stage: IngestStage): string {
  return STAGE_LABELS[stage];
}

export interface ProgressSegment {
  key: IngestStage;
  label: string;
  /** Fill ratio, 0..1. */
  fraction: number;
  /** True for the stage currently reporting. */
  active: boolean;
  /** Raw counts for the active segment's caption; undefined for inactive ones. */
  processed?: number;
  total?: number;
}

/**
 * Derive every segment's fill from the single stage the backend is reporting.
 *
 * Stages run in a fixed order, so anything before the reported one is finished
 * and anything after it has not started. Returns `null` when the stage is
 * unrecognised or absent (an `ingest_started` frame, or a backend that predates
 * stage reporting) — callers fall back to the plain single-track bar rather than
 * rendering a misleading empty segment strip.
 */
export function deriveSegments(
  stage: string | undefined,
  processed: number,
  total: number,
): ProgressSegment[] | null {
  if (!isIngestStage(stage)) return null;
  const activeIdx = INGEST_STAGES.indexOf(stage);

  return INGEST_STAGES.map((key, i) => {
    if (i < activeIdx) {
      return { key, label: STAGE_LABELS[key], fraction: 1, active: false };
    }
    if (i > activeIdx) {
      return { key, label: STAGE_LABELS[key], fraction: 0, active: false };
    }
    // A zero total means "nothing to do in this stage" — a clean rescan reports
    // analyzing as 0-of-0. That is COMPLETE, not empty; rendering it at 0% would
    // leave the bar looking stalled for the rest of the scan.
    const fraction = total > 0 ? Math.min(1, Math.max(0, processed / total)) : 1;
    return {
      key,
      label: STAGE_LABELS[key],
      fraction,
      active: true,
      processed,
      total,
    };
  });
}
