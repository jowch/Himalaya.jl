import { useEffect, useMemo, useRef, useState } from "react";
import { usePeaks, useSpeculativeSnap, useCreateSpeculative } from "../queries";
import { useFocusTrap } from "../hooks/useFocusTrap";
import { KNOWN_PHASES, phaseColor } from "../phases";
import { Button } from "./ui";

export interface SpeculativeBuilderProps {
  exposureId: number;
  onClose: () => void;
}

const DEFAULT_PHASE = "Lamellar";

/** Phases ordered by minpeaks ascending — low-bar phases first since the
 * builder exists primarily to let users hypothesize from sparse patterns. */
const PHASES_BY_MINPEAKS = [
  "Lamellar", "Hexagonal", "Square",
  "Pn3m", "Im3m", "Ia3d", "Fm3m", "Fd3m",
] as const;

export function SpeculativeBuilder({ exposureId, onClose }: SpeculativeBuilderProps): JSX.Element {
  const peaksQ = usePeaks(exposureId);
  // Memoize so the auto-pick effect doesn't re-fire on every render — without
  // this, the array identity churns and `[peaks, anchorPeakId]` deps thrash.
  const peaks  = useMemo(
    () => (peaksQ.data ?? []).filter((p) => !p.excluded).slice().sort((a, b) => a.q - b.q),
    [peaksQ.data],
  );

  const [phase, setPhase]                 = useState<string>(DEFAULT_PHASE);
  const [anchorPeakId, setAnchorPeakId]   = useState<number | null>(null);
  const [anchorRatio, setAnchorRatio]     = useState<number>(1);
  // ratio_position → checked (anchor's ratio is always implicitly true)
  const [included, setIncluded]           = useState<Record<number, boolean>>({});

  const dialogRef = useRef<HTMLDivElement>(null);
  useFocusTrap(dialogRef, true);

  // Reset selection when phase or anchor changes
  useEffect(() => { setIncluded({}); }, [phase, anchorPeakId, anchorRatio]);

  // Auto-pick the strongest peak as anchor on mount (highest sharpness).
  useEffect(() => {
    if (anchorPeakId === null && peaks.length > 0) {
      const strongest = peaks.reduce((acc, p) =>
        (p.sharpness ?? 0) > (acc.sharpness ?? 0) ? p : acc, peaks[0]!);
      setAnchorPeakId(strongest.id);
    }
  }, [peaks, anchorPeakId]);

  const snapQ = useSpeculativeSnap(exposureId, phase, anchorPeakId ?? undefined, anchorRatio);
  const createMut = useCreateSpeculative(exposureId);

  function toggleIncluded(rpos: number): void {
    setIncluded((prev) => ({ ...prev, [rpos]: !prev[rpos] }));
  }

  function handleSave(): void {
    if (anchorPeakId === null) return;
    const additional = (snapQ.data ?? [])
      .filter((s) => !s.is_anchor && s.suggested_peak_id !== null && included[s.ratio_position])
      .map((s) => ({ ratio_position: s.ratio_position, peak_id: s.suggested_peak_id! }));

    // Speculative indices land in Candidates (not the active set) by default;
    // promoting one is an explicit click on the Speculative section.
    createMut.mutate(
      { phase, anchor_peak_id: anchorPeakId, anchor_ratio: anchorRatio, additional, active: false },
      { onSuccess: onClose },
    );
  }

  // Esc to close
  useEffect(() => {
    function onKey(e: KeyboardEvent): void { if (e.key === "Escape") onClose(); }
    window.addEventListener("keydown", onKey);
    return () => { window.removeEventListener("keydown", onKey); };
  }, [onClose]);

  const color = phaseColor(phase);

  return (
    <div
      role="dialog"
      aria-modal="true"
      aria-label="Build speculative index"
      className="fixed inset-0 z-50 flex items-center justify-center bg-black/40 p-4"
      onClick={onClose}
    >
      <div
        ref={dialogRef}
        data-testid="speculative-builder"
        className="bg-bg border border-border-soft rounded-xl shadow-xl w-full max-w-md flex flex-col gap-4 p-5 max-h-[90vh] overflow-y-auto"
        onClick={(e) => { e.stopPropagation(); }}
      >
        <div className="flex items-center justify-between">
          <h2 className="text-lg font-semibold text-fg">Speculative index</h2>
          <button
            onClick={onClose}
            aria-label="Close"
            className="text-fg-dim hover:text-fg text-xl px-2 leading-none"
          >×</button>
        </div>

        <p className="text-xs text-fg-dim">
          Hand-build an index from a phase guess and an anchor peak. Best for sparse
          patterns where auto-analysis falls below the {""}
          <span className="font-mono">minpeaks</span> bar.
        </p>

        {/* Phase picker */}
        <div className="flex flex-col gap-1.5">
          <label htmlFor="spec-phase" className="text-xs text-fg-dim font-semibold">Phase</label>
          <select
            id="spec-phase"
            data-testid="spec-phase-select"
            className="bg-bg-hover border border-border-soft rounded-md p-1.5 text-fg"
            value={phase}
            onChange={(e) => { setPhase(e.target.value); }}
          >
            {PHASES_BY_MINPEAKS.filter((p) => KNOWN_PHASES.includes(p)).map((p) => (
              <option key={p} value={p}>{p}</option>
            ))}
          </select>
        </div>

        {/* Anchor peak picker */}
        <div className="flex flex-col gap-1.5">
          <label htmlFor="spec-anchor" className="text-xs text-fg-dim font-semibold">
            Anchor peak (q-value)
          </label>
          <select
            id="spec-anchor"
            data-testid="spec-anchor-select"
            className="bg-bg-hover border border-border-soft rounded-md p-1.5 text-fg font-mono"
            value={anchorPeakId ?? ""}
            onChange={(e) => { setAnchorPeakId(Number(e.target.value)); }}
          >
            <option value="" disabled>Choose a peak…</option>
            {peaks.map((p) => (
              <option key={p.id} value={p.id}>
                q = {p.q.toFixed(4)} ({p.source})
              </option>
            ))}
          </select>
        </div>

        {/* Anchor ratio position */}
        <div className="flex flex-col gap-1.5">
          <label htmlFor="spec-ratio" className="text-xs text-fg-dim font-semibold">
            Anchor ratio position (1 = first peak of the series)
          </label>
          <input
            id="spec-ratio"
            type="number"
            min={1}
            max={20}
            step={1}
            value={anchorRatio}
            onChange={(e) => { setAnchorRatio(Math.max(1, Number(e.target.value) || 1)); }}
            className="bg-bg-hover border border-border-soft rounded-md p-1.5 text-fg font-mono w-24"
          />
        </div>

        {/* Snap preview */}
        <div className="flex flex-col gap-1.5">
          <span className="text-xs text-fg-dim font-semibold">Predicted ratio positions</span>
          {anchorPeakId === null ? (
            <p className="text-xs text-fg-dim italic">Pick an anchor peak to see predictions.</p>
          ) : snapQ.isLoading ? (
            <p className="text-xs text-fg-dim italic">Computing…</p>
          ) : snapQ.error ? (
            <p className="text-xs text-error">Snap failed: {String(snapQ.error)}</p>
          ) : (
            <ul className="flex flex-col gap-1 max-h-56 overflow-y-auto rounded-md border border-border-soft p-2">
              {(snapQ.data ?? []).map((s) => {
                const checked = s.is_anchor || (s.suggested_peak_id !== null && !!included[s.ratio_position]);
                const disabled = s.is_anchor || s.suggested_peak_id === null;
                const dq = s.suggested_residual !== null && s.predicted_q > 0
                  ? (s.suggested_residual / s.predicted_q * 100).toFixed(2)
                  : null;
                return (
                  <li
                    key={s.ratio_position}
                    data-testid={`spec-snap-row-${s.ratio_position}`}
                    className="flex items-center gap-2 text-xs"
                  >
                    <input
                      type="checkbox"
                      checked={checked}
                      disabled={disabled}
                      onChange={() => { toggleIncluded(s.ratio_position); }}
                      aria-label={`ratio position ${s.ratio_position}`}
                    />
                    <span className="font-mono text-fg-muted w-10">r{s.ratio_position}</span>
                    <span className="font-mono text-fg-dim w-24">
                      ≈ q {s.predicted_q.toFixed(4)}
                    </span>
                    {s.is_anchor ? (
                      <span className="px-1.5 rounded text-xs"
                            style={{ background: `color-mix(in oklab, ${color} 20%, transparent)`, color }}>
                        anchor
                      </span>
                    ) : s.suggested_peak_id !== null ? (
                      <span className="font-mono text-fg-muted truncate">
                        snap q {s.suggested_q?.toFixed(4)} (Δ {dq}%)
                      </span>
                    ) : (
                      <span className="text-fg-dim italic">no peak</span>
                    )}
                  </li>
                );
              })}
            </ul>
          )}
        </div>

        {createMut.error && (
          <p className="text-xs text-error" role="alert">
            {String((createMut.error as Error).message ?? createMut.error)}
          </p>
        )}

        <div className="flex justify-end gap-2 pt-2">
          <Button variant="ghost" onClick={onClose} type="button">Cancel</Button>
          <Button
            variant="primary"
            onClick={handleSave}
            type="button"
            disabled={anchorPeakId === null || createMut.isPending}
            data-testid="spec-save-button"
          >
            {createMut.isPending ? "Saving…" : "Save"}
          </Button>
        </div>
      </div>
    </div>
  );
}
