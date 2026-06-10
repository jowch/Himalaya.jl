import { FlagButton, GripHandle } from "../ui";
import { Sparkline } from "../plot/Sparkline";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface ScopeSampleRowProps {
  name: string;
  sampleId: string;
  /** Sparkline data. */
  trace: { q: number[]; I: number[] };
  /** Dominant phase → spark hue; null/undefined → unindexed. */
  phase?: string | null;
  /** The parsed ordering value, e.g. "1 : 0.25". */
  value: string;
  /** Skipped read → FlagButton flagged look + faint accent row wash. */
  flagged?: boolean;
  /** Forwarded to FlagButton onClick. */
  onToggleFlag?: () => void;
  /** PLACEMENT-ONLY. */
  className?: string;
}

/**
 * One row of the series-scoping worksheet (mockup `.srow`): a drag grip, a small
 * trace sparkline, the sample's name + id, and a trailing parsed-value control.
 *
 * When `flagged` (the user skipped this read from the batch write), the whole
 * row takes a faint accent wash and the value control marks the exclusion. The row
 * carries `group` so the grip brightens on row hover; it also draws its own
 * bottom hairline (the parent strips the last). Drag-reorder is page-deferred —
 * the grip here is the visual affordance only.
 */
export function ScopeSampleRow({
  name,
  sampleId,
  trace,
  phase,
  value,
  flagged,
  onToggleFlag,
  className,
}: ScopeSampleRowProps): JSX.Element {
  return (
    <div
      data-testid="scope-sample-row"
      data-flagged={flagged ? "true" : "false"}
      className={cx(
        "group flex items-center gap-3 px-2 py-2.5 border-b border-hair",
        flagged && "bg-accent/5",
        className,
      )}
    >
      <GripHandle />
      <Sparkline trace={trace} {...(phase != null ? { phase } : {})} />
      <div className="flex-1 min-w-0">
        <div className="text-body font-semibold text-ink">{name}</div>
        <div className="text-caption font-mono text-ink-soft">{sampleId}</div>
      </div>
      <FlagButton
        value={value}
        {...(flagged ? { flagged: true } : {})}
        {...(onToggleFlag ? { onClick: onToggleFlag } : {})}
      />
    </div>
  );
}
