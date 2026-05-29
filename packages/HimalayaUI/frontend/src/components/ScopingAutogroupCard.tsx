interface Props {
  memberCount: number;
  keyLabel: string;
  flagCount: number;
}

/**
 * The machine's grouping summary (series-scoping.html `.autogroup`): one
 * confident breath naming what Himalaya did — grouped N samples, read the
 * order from <key>, M parsed cleanly, K need a look.
 */
export function ScopingAutogroupCard({ memberCount, keyLabel, flagCount }: Props): JSX.Element {
  const clean = memberCount - flagCount;
  const flagPhrase =
    flagCount === 0
      ? "all parsed cleanly"
      : flagCount === 1
        ? `${clean} parsed cleanly, one needs a look`
        : `${clean} parsed cleanly, ${flagCount} need a look`;
  return (
    <div
      data-testid="scoping-autogroup"
      className="mt-4 flex items-start gap-2.5 rounded-md border border-hair bg-paper-sunk px-3 py-3"
    >
      <svg viewBox="0 0 16 16" className="mt-px h-[15px] w-[15px] shrink-0" aria-hidden>
        <path
          d="M8 1.4l1.7 4.2 4.5.3-3.5 2.9 1.2 4.4L8 10.9 4.1 13.2l1.2-4.4L1.8 5.9l4.5-.3z"
          fill="var(--color-print-accent)"
        />
      </svg>
      <p className="text-xs leading-relaxed text-ink-soft">
        You selected <b className="font-semibold text-ink">{memberCount} samples</b>. Himalaya
        grouped them from their names and read the order from{" "}
        <b className="font-semibold text-ink">{keyLabel}</b>: {flagPhrase}.
      </p>
    </div>
  );
}
