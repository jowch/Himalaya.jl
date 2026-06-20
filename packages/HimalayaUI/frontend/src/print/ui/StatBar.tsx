import type { ReactNode } from "react";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface StatBarStat {
  /** Stable react key + test handle. */
  key: string;
  /** Uppercase micro-caption under the value. */
  caption: string;
  /** The value (already formatted by the caller; mono unless pending). */
  value: ReactNode;
  /** When true, the value renders as a faint italic placeholder
   *  ("discovered-so-far" / span-pending; spec §8.3). */
  pending?: boolean;
}

export interface StatBarProps {
  /** Required: names the band for assistive tech. */
  "aria-label": string;
  stats: ReadonlyArray<StatBarStat>;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * StatBar — the quiet stat ledger (DESIGN: hairline-divided cells, no heavy
 * rule). Each cell is a mono value over an uppercase micro-caption; a thin left
 * hairline divides cells (first has none). A `pending` cell italicizes a faint
 * placeholder for the live-ingest discovered-so-far state. Closed look; the
 * consumer's `className` is placement-only.
 */
export function StatBar({
  "aria-label": ariaLabel,
  stats,
  className = "",
}: StatBarProps): JSX.Element {
  return (
    <div
      data-testid="statbar"
      role="group"
      aria-label={ariaLabel}
      className={cx("flex", className)}
    >
      {stats.map((s, i) => (
        <div
          key={s.key}
          data-testid="statbar-cell"
          data-pending={s.pending ? "true" : undefined}
          className={cx(
            "flex flex-col gap-1 px-7 py-0.5",
            i === 0 ? "pl-0" : "border-l border-hair",
          )}
        >
          <span
            className={cx(
              "leading-none",
              s.pending
                ? "text-sm italic text-ink-faint"
                : "font-mono text-xl text-ink",
            )}
          >
            {s.value}
          </span>
          <span className="text-xs font-bold uppercase tracking-wide text-ink-faint">
            {s.caption}
          </span>
        </div>
      ))}
    </div>
  );
}
