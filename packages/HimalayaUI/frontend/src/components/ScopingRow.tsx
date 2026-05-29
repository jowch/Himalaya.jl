import type { Trace } from "../api";
import type { OrderingRow } from "../lib/scoping/proposeOrdering";
import { ScopingSparkline } from "./ScopingSparkline";
import { ScopingValueCell } from "./ScopingValueCell";

interface Props {
  row: OrderingRow;
  trace: Trace | undefined;
  phase: string | null;
  onChangeValue: (value: string) => void;
  onToggleFlag: () => void;
}

/**
 * One member row of the scoping worksheet (series-scoping.html `.srow`):
 * the trace sparkline, the sample name + id, and the confirm-not-fill-out
 * value cell. Flagged rows get an amber row tint. (M2 removed a decorative
 * drag grip: scoping rows are ordered by parsed tag value, not by hand, so a
 * grip reordered nothing and contradicted the value-sort model.)
 */
export function ScopingRow({ row, trace, phase, onChangeValue, onToggleFlag }: Props): JSX.Element {
  return (
    <div
      data-testid={`scoping-row-${row.sampleId}`}
      data-flagged={row.flagged ? "true" : undefined}
      className={`flex items-center gap-3 border-b border-hair px-2 py-2 last:border-b-0 ${
        row.flagged ? "bg-print-accent/5" : ""
      }`}
    >
      <ScopingSparkline trace={trace} phase={phase} />
      <span className="min-w-0 flex-1">
        <span className="block truncate text-base font-semibold text-ink">{row.sampleName}</span>
        <span className="block font-mono text-xs text-ink-faint">smp_{row.sampleId}</span>
      </span>
      <ScopingValueCell
        sampleId={row.sampleId}
        sampleName={row.sampleName}
        value={row.value}
        flagged={row.flagged}
        onChange={onChangeValue}
        onToggleFlag={onToggleFlag}
      />
    </div>
  );
}
