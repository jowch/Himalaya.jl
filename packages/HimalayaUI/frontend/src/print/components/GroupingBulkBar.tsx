import type { JSX } from "react";
import { Button } from "../ui/Button";

export interface GroupingBulkBarProps {
  count: number;
  /** singular noun; pluralized with a naive +"s". */
  noun: string;
  primaryLabel: string;
  /** Primary action is disabled below this many selected (merge needs >=2). */
  primaryEnabled: boolean;
  onPrimary: () => void;
  onClear: () => void;
  className?: string;
}

/** Dark floating bulk-action bar for the grouping selection. Mirrors the
 *  CullBar's dark `bg-ink` idiom but is parameterized for the merge gesture.
 *  Presentational; the page owns the selection + handlers. */
export function GroupingBulkBar(p: GroupingBulkBarProps): JSX.Element {
  return (
    <div
      data-testid="bulk-bar"
      // bottom-20 + z-40 clears the grouping Confirm footer (fixed bottom-0
      // z-30) so the bar isn't covered by it.
      className={`fixed inset-x-0 bottom-20 z-40 mx-auto flex w-fit items-center gap-4 rounded-md bg-ink px-5 py-3 text-paper shadow-lg${p.className ? ` ${p.className}` : ""}`}
    >
      <span className="font-mono text-sm">
        {p.count} {p.noun}{p.count === 1 ? "" : "s"} selected
      </span>
      <Button variant="accent" disabled={!p.primaryEnabled} onClick={p.onPrimary}>{p.primaryLabel}</Button>
      <Button variant="ghostInverse" onClick={p.onClear}>Clear</Button>
    </div>
  );
}
