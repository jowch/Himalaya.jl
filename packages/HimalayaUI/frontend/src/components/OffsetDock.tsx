interface OffsetDockProps {
  /** Only renders when true — shown only in full-bleed (rail collapsed). */
  show: boolean;
  value: number;
  onChange: (next: number) => void;
}

/**
 * OffsetDock — floating offset control kept reachable when the rail is
 * collapsed into full-bleed (R8 / B-G). Mockup `series-builder.html` `.dock`:
 * a plate-surfaced pill pinned bottom-right with an "Offset" label, the
 * slider, and the live value. Mirrors the rail OffsetSlider's range/step so
 * the page can keep both in sync.
 */
export function OffsetDock({ show, value, onChange }: OffsetDockProps): JSX.Element | null {
  if (!show) return null;
  return (
    <div
      data-testid="offset-dock"
      className="fixed bottom-6 right-6 z-10 flex items-center gap-3.5 rounded-xl border border-hair-strong bg-plate px-4 py-3 shadow-[0_8px_26px_-10px_rgba(60,52,40,.34)]"
    >
      <span className="text-[10px] font-bold uppercase tracking-wide text-ink-faint">Offset</span>
      <input
        type="range"
        data-testid="offset-dock-slider"
        aria-label="Trace offset"
        min="0.4"
        max="1.4"
        step="0.05"
        value={value}
        onChange={(e) => onChange(parseFloat(e.target.value))}
        className="w-32 cursor-pointer appearance-none rounded-full bg-hair-strong accent-print-accent"
      />
      <span
        data-testid="offset-dock-value"
        className="min-w-[40px] text-xs font-bold text-ink"
      >
        {value.toFixed(2)}&times;
      </span>
    </div>
  );
}
