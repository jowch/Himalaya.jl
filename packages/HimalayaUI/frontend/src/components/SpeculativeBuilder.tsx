import { useEffect, useMemo, useRef, useState } from "react";
import type { SpeculativeSnap } from "../api";
import { usePeaks, useSpeculativeSnap, useCreateSpeculative } from "../queries";
import { useExposureHasPendingPeakOps } from "../lib/queue/hooks";
import { KNOWN_PHASES, phaseColor } from "../phases";
import { Button, ModalShell } from "./ui";

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

  // M2.4 gating: while a peak op is pending the snap query is disabled, but we
  // keep showing the last good snap with an "updating to latest…" subtext so
  // the dialog doesn't blank out mid-edit.
  const blocked = useExposureHasPendingPeakOps(exposureId);
  const snapQ = useSpeculativeSnap(exposureId, phase, anchorPeakId ?? undefined, anchorRatio);
  const lastGoodSnapRef = useRef<SpeculativeSnap[] | null>(null);
  if (snapQ.data !== undefined) {
    lastGoodSnapRef.current = snapQ.data;
  }
  const displaySnap: SpeculativeSnap[] | null =
    snapQ.data ?? lastGoodSnapRef.current;

  const createMut = useCreateSpeculative(exposureId);

  // M2.4: useQueueMutation does not expose per-call onSuccess; observe
  // `isSuccess` (set after the mutationFn resolves) and close the dialog
  // exactly once per successful save.
  const closedOnSuccessRef = useRef(false);
  useEffect(() => {
    if (createMut.isSuccess && !closedOnSuccessRef.current) {
      closedOnSuccessRef.current = true;
      onClose();
    }
  }, [createMut.isSuccess, onClose]);

  function toggleIncluded(rpos: number): void {
    setIncluded((prev) => ({ ...prev, [rpos]: !prev[rpos] }));
  }

  function handleSave(): void {
    if (anchorPeakId === null) return;
    const additional = (displaySnap ?? [])
      .filter((s) => !s.is_anchor && s.suggested_peak_id !== null && included[s.ratio_position])
      .map((s) => ({ ratio_position: s.ratio_position, peak_id: s.suggested_peak_id! }));

    // Speculative indices land in Candidates (not the active set) by default;
    // promoting one is an explicit click on the Speculative section.
    createMut.mutate({
      phase,
      anchor_peak_id: anchorPeakId,
      anchor_ratio: anchorRatio,
      additional,
      active: false,
    });
  }

  const color = phaseColor(phase);
  const errorMessage: string | null = createMut.error
    ? (createMut.error instanceof Error
        ? createMut.error.message
        : String(createMut.error))
    : null;

  return (
    <ModalShell
      open
      onClose={onClose}
      size="sm"
      aria-label="Build speculative index"
      testId="speculative-builder"
      className="max-h-[90vh] overflow-y-auto gap-4 p-5"
    >
      <div className="flex items-center justify-between">
        <h2 className="text-lg font-semibold text-ink">Speculative index</h2>
        <button
          onClick={onClose}
          aria-label="Close"
          className="text-ink-faint hover:text-ink text-xl px-2 leading-none"
        >×</button>
      </div>

        <p className="text-xs text-ink-faint">
          Hand-build an index from a phase guess and an anchor peak. Best for sparse
          patterns where auto-analysis falls below the {""}
          <span className="font-mono">minpeaks</span> bar.
        </p>

        {/* Phase picker */}
        <div className="flex flex-col gap-1.5">
          <label htmlFor="spec-phase" className="text-xs text-ink-faint font-semibold">Phase</label>
          <select
            id="spec-phase"
            data-testid="spec-phase-select"
            className="bg-paper-sunk border border-hair rounded-md p-1.5 text-ink"
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
          <label htmlFor="spec-anchor" className="text-xs text-ink-faint font-semibold">
            Anchor peak (q-value)
          </label>
          <select
            id="spec-anchor"
            data-testid="spec-anchor-select"
            className="bg-paper-sunk border border-hair rounded-md p-1.5 text-ink font-mono"
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
          <label htmlFor="spec-ratio" className="text-xs text-ink-faint font-semibold">
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
            className="bg-paper-sunk border border-hair rounded-md p-1.5 text-ink font-mono w-24"
          />
        </div>

        {/* Snap preview */}
        <div className="flex flex-col gap-1.5">
          <span className="text-xs text-ink-faint font-semibold">Predicted ratio positions</span>
          {anchorPeakId === null ? (
            <p className="text-xs text-ink-faint italic">Pick an anchor peak to see predictions.</p>
          ) : displaySnap === null && snapQ.isLoading ? (
            <p className="text-xs text-ink-faint italic">Computing…</p>
          ) : snapQ.error ? (
            <p className="text-xs text-error">
              Snap failed: {snapQ.error instanceof Error ? snapQ.error.message : String(snapQ.error)}
            </p>
          ) : (
            <>
              {blocked && (
                <p
                  data-testid="spec-snap-updating"
                  className="text-xs text-ink-faint italic"
                >
                  Updating to latest…
                </p>
              )}
              <ul className="flex flex-col gap-1 max-h-56 overflow-y-auto rounded-md border border-hair p-2">
                {(displaySnap ?? []).map((s) => {
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
                      <span className="font-mono text-ink-soft w-10">r{s.ratio_position}</span>
                      <span className="font-mono text-ink-faint w-24">
                        ≈ q {s.predicted_q.toFixed(4)}
                      </span>
                      {s.is_anchor ? (
                        <span className="px-1.5 rounded text-xs"
                              style={{ background: `color-mix(in oklab, ${color} 20%, transparent)`, color }}>
                          anchor
                        </span>
                      ) : s.suggested_peak_id !== null ? (
                        <span className="font-mono text-ink-soft truncate">
                          snap q {s.suggested_q?.toFixed(4)} (Δ {dq}%)
                        </span>
                      ) : (
                        <span className="text-ink-faint italic">no peak</span>
                      )}
                    </li>
                  );
                })}
              </ul>
            </>
          )}
        </div>

        {errorMessage !== null && (
          <p className="text-xs text-error" role="alert">
            {errorMessage}
          </p>
        )}

        <div className="flex justify-end gap-2 pt-2">
          <Button variant="ghost" onClick={onClose} type="button">Cancel</Button>
          <Button
            variant="solid"
            onClick={handleSave}
            type="button"
            disabled={anchorPeakId === null || createMut.isPending}
            data-testid="spec-save-button"
          >
            {createMut.isPending ? "Saving…" : "Save"}
          </Button>
        </div>
    </ModalShell>
  );
}
