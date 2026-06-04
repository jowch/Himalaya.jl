import { ModalShell, Button, IconButton, PhaseChip } from "../ui";
import type { SpeculativeSnap } from "../../api";

export interface SpeculativeDialogProps {
  open: boolean;
  onClose: () => void;
  // phase picker
  phases: readonly string[];
  phase: string;
  onPhaseChange: (phase: string) => void;
  // anchor peak picker
  peaks: { id: number; q: number; source: string }[];
  anchorPeakId: number | undefined; // undefined → "Choose a peak…" placeholder selected
  onAnchorChange: (peakId: number) => void;
  // anchor ratio
  anchorRatio: number;
  onAnchorRatioChange: (ratio: number) => void;
  // snap preview rows
  snap: SpeculativeSnap[];
  included: Record<number, boolean>; // ratio_position -> checked
  onToggleIncluded: (ratioPosition: number) => void;
  snapLoading?: boolean; // → "Computing…"
  snapError?: string | null; // → "Snap failed: …"
  blocked?: boolean; // → "Updating to latest…" subtext
  // save
  saving?: boolean; // Save → "Saving…" + disabled
  error?: string | null; // create error → role="alert"
  onCreate: () => void;
  className?: string; // PLACEMENT ONLY
}

export function SpeculativeDialog(props: SpeculativeDialogProps): JSX.Element | null {
  const {
    open,
    onClose,
    phases,
    phase,
    onPhaseChange,
    peaks,
    anchorPeakId,
    onAnchorChange,
    anchorRatio,
    onAnchorRatioChange,
    snap,
    included,
    onToggleIncluded,
    snapLoading = false,
    snapError = null,
    blocked = false,
    saving = false,
    error = null,
    onCreate,
    className,
  } = props;

  return (
    <ModalShell
      open={open}
      onClose={onClose}
      size="sm"
      aria-label="Build speculative index"
      testId="speculative-builder"
      className={`max-h-[90vh] overflow-y-auto gap-4 p-5${className ? ` ${className}` : ""}`}
    >
      <div className="flex items-center justify-between">
        <h2 className="text-lg font-semibold text-ink">Speculative index</h2>
        <IconButton label="Close" dismiss tone="ghost" onClick={onClose} />
      </div>

      <p className="text-xs text-ink-faint">
        Hand-build an index from a phase guess and an anchor peak. Best for sparse
        patterns where auto-analysis falls below the{" "}
        <span className="font-mono">minpeaks</span> bar.
      </p>

      {/* Phase picker */}
      <div className="flex flex-col gap-1.5">
        <label htmlFor="spec-phase" className="text-xs text-ink-faint font-semibold">
          Phase
        </label>
        <select
          id="spec-phase"
          data-testid="spec-phase-select"
          className="bg-paper-sunk border border-hair rounded-md p-1.5 text-ink"
          value={phase}
          onChange={(e) => {
            onPhaseChange(e.target.value);
          }}
        >
          {phases.map((p) => (
            <option key={p} value={p}>
              {p}
            </option>
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
          onChange={(e) => {
            onAnchorChange(Number(e.target.value));
          }}
        >
          <option value="" disabled>
            Choose a peak…
          </option>
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
          onChange={(e) => {
            onAnchorRatioChange(Math.max(1, Number(e.target.value) || 1));
          }}
          className="bg-paper-sunk border border-hair rounded-md p-1.5 text-ink font-mono w-24"
        />
      </div>

      {/* Snap preview */}
      <div className="flex flex-col gap-1.5">
        <span className="text-xs text-ink-faint font-semibold">Predicted ratio positions</span>
        {anchorPeakId === undefined ? (
          <p className="text-xs text-ink-faint italic">Pick an anchor peak to see predictions.</p>
        ) : snapLoading && snap.length === 0 ? (
          <p className="text-xs text-ink-faint italic">Computing…</p>
        ) : snapError ? (
          <p className="text-xs text-error">Snap failed: {snapError}</p>
        ) : (
          <>
            {blocked && (
              <p data-testid="spec-snap-updating" className="text-xs text-ink-faint italic">
                Updating to latest…
              </p>
            )}
            <ul className="flex flex-col gap-1 max-h-56 overflow-y-auto rounded-md border border-hair p-2">
              {snap.map((s) => {
                const checked =
                  s.is_anchor || (s.suggested_peak_id !== null && !!included[s.ratio_position]);
                const disabled = s.is_anchor || s.suggested_peak_id === null;
                const dq =
                  s.suggested_residual !== null && s.predicted_q > 0
                    ? ((s.suggested_residual / s.predicted_q) * 100).toFixed(2)
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
                      onChange={() => {
                        onToggleIncluded(s.ratio_position);
                      }}
                      aria-label={`ratio position ${s.ratio_position}`}
                    />
                    <span className="font-mono text-ink-soft w-10">r{s.ratio_position}</span>
                    <span className="font-mono text-ink-faint w-24">
                      ≈ q {s.predicted_q.toFixed(4)}
                    </span>
                    {s.is_anchor ? (
                      <span className="inline-flex items-center gap-1">
                        <PhaseChip phase={phase} size="sm" aria-hidden="true" />
                        <span className="text-ink-soft">anchor</span>
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

      {error !== null && (
        <p className="text-xs text-error" role="alert">
          {error}
        </p>
      )}

      <div className="flex justify-end gap-2 pt-2">
        <Button variant="ghost" onClick={onClose} type="button">
          Cancel
        </Button>
        <Button
          variant="solid"
          onClick={onCreate}
          type="button"
          disabled={anchorPeakId === undefined || saving}
          data-testid="spec-save-button"
        >
          {saving ? "Saving…" : "Save"}
        </Button>
      </div>
    </ModalShell>
  );
}
