/**
 * The floating cull bar (M-3) — sample-table.html's `.cullbar`. A fixed,
 * bottom-centre ink pill that appears while a frame selection exists, offering
 * batch reject / clear. Restyled from the inline action bar; selection stays
 * per-sample by design (it never crosses rows — see ContactSheetRow), but the
 * bar reads as the mockup's floating control.
 *
 * Rendered conditionally by the owning row, so at most one bar shows at a time
 * (clicking a frame in another row clears the prior row's selection on its own
 * local state; only the row with a live selection mounts the bar).
 */
export function CullBar({
  count,
  onReject,
  onClear,
}: {
  count: number;
  onReject: () => void;
  onClear: () => void;
}): JSX.Element {
  return (
    <div
      data-testid="cull-bar"
      className="fixed bottom-7 left-1/2 z-30 flex -translate-x-1/2 items-center
                 gap-1 rounded-[10px] bg-ink py-[7px] pl-4 pr-2 text-paper
                 shadow-[0_18px_40px_-14px_rgba(40,33,24,.5)]"
    >
      <span className="mr-2.5 text-[12.5px] font-semibold">
        <b data-testid="cull-count" className="font-mono">
          {count}
        </b>{" "}
        frames selected
      </span>
      <button
        type="button"
        data-testid="batch-reject"
        onClick={onReject}
        className="rounded-[7px] bg-print-accent px-3 py-[7px] text-xs font-semibold
                   text-paper hover:opacity-90"
      >
        Drop
        <span className="ml-1 font-mono text-[9.5px] opacity-60">X</span>
      </button>
      <button
        type="button"
        data-testid="batch-clear"
        onClick={onClear}
        className="rounded-[7px] px-3 py-[7px] text-xs font-semibold text-paper/70
                   hover:text-paper"
      >
        Clear
        <span className="ml-1 font-mono text-[9.5px] opacity-60">Esc</span>
      </button>
    </div>
  );
}
