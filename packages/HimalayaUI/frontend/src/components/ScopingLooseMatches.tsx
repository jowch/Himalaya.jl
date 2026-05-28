import { useState } from "react";
import type { Trace } from "../api";
import type { OrderingRow } from "../lib/scoping/proposeOrdering";
import { ScopingSparkline } from "./ScopingSparkline";

interface Props {
  rows: OrderingRow[];
  traces: Map<number, Trace>;
  phases: Map<number, string | null>;
  onAdd: (sampleId: number) => void; // fold a loose match into the series
}

// Keep the worksheet legible: the corpus picker returns every sample (there is
// no contact-sheet selection plumbing yet), so loose matches can run long.
// Show a tidy head and tuck the rest behind a "show N more" — the mockup's
// candidate list is short and confident.
const COLLAPSED_COUNT = 4;

/**
 * "Himalaya also found" (series-scoping.html `.candidates`): samples the machine
 * noticed but would not assume into the series — here, corpus samples that lack
 * a value for the ordering variable. Each is addable in one click; coherence is
 * the human's call.
 */
export function ScopingLooseMatches({ rows, traces, phases, onAdd }: Props): JSX.Element {
  const [expanded, setExpanded] = useState(false);
  const visible = expanded ? rows : rows.slice(0, COLLAPSED_COUNT);
  const hiddenCount = rows.length - visible.length;
  return (
    <div className="mt-5 border-t border-hair pt-4">
      <div className="mb-2 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">
        Himalaya also found
      </div>
      {rows.length === 0 ? (
        <div data-testid="scoping-loose-empty" className="text-[11.5px] italic text-ink-faint">
          Nothing else in the corpus matches this grouping.
        </div>
      ) : (
        <div className="flex flex-col gap-2">
          {visible.map((r) => (
            <div
              key={r.sampleId}
              data-testid={`scoping-loose-${r.sampleId}`}
              className="flex items-center gap-3 rounded-md border border-dashed border-hair-strong px-2 py-2"
            >
              <span className="inline-flex shrink-0 opacity-70">
                <ScopingSparkline trace={traces.get(r.sampleId)} phase={phases.get(r.sampleId) ?? null} />
              </span>
              <div className="min-w-0 flex-1">
                <div className="truncate text-[12.5px] font-semibold text-ink-soft">{r.sampleName}</div>
                <div className="text-[11px] text-ink-faint">
                  No <span className="font-semibold text-print-accent">value</span> for the ordering
                  variable — add it to include.
                </div>
              </div>
              <button
                type="button"
                data-testid={`scoping-loose-add-${r.sampleId}`}
                onClick={() => onAdd(r.sampleId)}
                className="shrink-0 rounded-md border border-hair-strong bg-plate px-2.5 py-1.5 text-[11.5px] font-semibold text-ink hover:bg-paper-sunk"
              >
                + Add to series
              </button>
            </div>
          ))}
          {hiddenCount > 0 ? (
            <button
              type="button"
              data-testid="scoping-loose-more"
              onClick={() => setExpanded(true)}
              className="self-start text-[11.5px] font-semibold text-print-accent hover:underline"
            >
              Show {hiddenCount} more
            </button>
          ) : null}
        </div>
      )}
    </div>
  );
}
