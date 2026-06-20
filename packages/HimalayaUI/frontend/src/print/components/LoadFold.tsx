import type { JSX } from "react";
import type { Load, LoadSample } from "../../api";
import { SampleFold } from "./SampleFold";

export interface LoadFoldProps {
  load: Load;
  open: boolean;
  /** Samples to render (already filtered by the page). */
  visibleSamples: LoadSample[];
  openSamples: Set<number>;
  selected: Set<number>;
  onToggleLoad: (loadId: number) => void;
  onToggleSampleOpen: (sampleId: number) => void;
  onToggleSelect: (sampleId: number) => void;
  onRename: (sampleId: number) => void;
  onSplit: (sampleId: number) => void;
  onMerge: (loserId: number, survivorId: number) => void;
  onDismissFlag: (sampleId: number) => void;
  onMoveExposure: (sampleId: number, exposureId: number) => void;
  thumbSrcFor: (exposureId: number) => string | null;
  className?: string;
}

export function LoadFold(p: LoadFoldProps): JSX.Element {
  const { load: l, open } = p;
  const flaggedCount = l.samples.filter((s) => s.flag).length;
  const clean = flaggedCount === 0;
  const totalExposures = l.samples.reduce((a, s) => a + s.exposures.length, 0);
  const time = l.start_time && l.end_time ? `${l.start_time}–${l.end_time}` : "";
  const hidden = l.samples.length - p.visibleSamples.length;

  return (
    <div
      data-testid="load-fold"
      className={[
        "mb-3 overflow-hidden rounded-md border bg-plate",
        clean ? "border-hair-strong" : "border-warning",
        p.className ?? "",
      ].filter(Boolean).join(" ")}
    >
      <button
        type="button"
        aria-expanded={open}
        aria-label={`Load ${l.load_index}`}
        className="flex w-full items-center gap-3 px-4 py-3.5 text-left"
        onClick={() => p.onToggleLoad(l.load_id)}
      >
        <span className="min-w-0 flex-1">
          <span className="text-headline text-ink">Load {l.load_index}</span>
          {hidden > 0 ? (
            <span className="ml-2 text-xs text-ink-faint">
              {"·"} {p.visibleSamples.length} of {l.samples.length} shown
            </span>
          ) : null}
          <span className="block font-mono text-xs text-ink-faint">
            {l.samples.length} samples{" "}{"·"} {totalExposures} exposures
            {time ? ` ${"·"} ${time}` : ""}
          </span>
        </span>
        <span className={`ml-auto text-xs font-semibold${clean ? " text-success" : " text-warning"}`}>
          {clean ? "✓ grouped cleanly" : `⚠ ${flaggedCount} to check`}
        </span>
      </button>

      {open ? (
        <div className="border-t border-hair bg-paper-sunk p-1.5">
          {p.visibleSamples.map((s) => (
            <div key={s.sample_id} className="m-1.5">
              <SampleFold
                sample={s}
                open={p.openSamples.has(s.sample_id) || !!s.flag}
                selected={p.selected.has(s.sample_id)}
                onToggleOpen={p.onToggleSampleOpen}
                onToggleSelect={p.onToggleSelect}
                onRename={p.onRename}
                onSplit={p.onSplit}
                onMerge={p.onMerge}
                onDismissFlag={p.onDismissFlag}
                onMoveExposure={p.onMoveExposure}
                thumbSrcFor={p.thumbSrcFor}
              />
            </div>
          ))}
        </div>
      ) : null}
    </div>
  );
}
