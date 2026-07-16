import type { ReactNode } from "react";
import { cx } from "../../lib/cx";

export interface MetaEntry {
  key: string;
  value: ReactNode;
}

interface MetaListProps {
  entries: MetaEntry[];
  /** Placement-only className (spacing, max-width, etc). Appearance is fixed. */
  className?: string;
}


/** Mono key/value list — measured values / ids / timestamps paired with their
 *  labels. Semantic <dl>/<dt>/<dd> carries the term/definition pairing into the
 *  a11y tree (the accessibility win). Monospace is correct here (rule E): these
 *  are measured values, not prose.
 *
 *  Impeccable: A/B/C/D N/A — non-interactive (no states, no touch target) and no
 *  color-coded meaning (key/value distinguished by position + faint/ink tone, a
 *  second channel). H — flat, no plate. */
export function MetaList({
  entries,
  className = "",
}: MetaListProps): JSX.Element {
  return (
    <dl
      data-testid="meta-list"
      className={cx("flex flex-col gap-1.5 font-mono text-sm", className)}
    >
      {entries.map((e) => (
        <div key={e.key} className="flex justify-between gap-4">
          <dt className="text-ink-soft">{e.key}</dt>
          {/* m-0 neutralizes the UA default margin-inline-start on <dd> so
              justify-between aligns the value flush to the right edge. */}
          <dd className="text-ink m-0">{e.value}</dd>
        </div>
      ))}
    </dl>
  );
}
