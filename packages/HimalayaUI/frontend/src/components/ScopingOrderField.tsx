interface Props {
  keyLabel: string;
}

/**
 * The "Ordered by" field (series-scoping.html `.order-field`): the one real
 * decision on this surface. Read-only in this milestone — re-selecting the
 * ordering variable is a follow-up; the field still presents as a decision.
 */
export function ScopingOrderField({ keyLabel }: Props): JSX.Element {
  return (
    <>
      <div className="mb-1.5 mt-5 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">
        Ordered by
      </div>
      <div
        data-testid="scoping-order-field"
        className="flex items-center justify-between rounded-md border border-hair-strong bg-plate px-3.5 py-3"
      >
        <span className="text-[15px] font-semibold text-ink">{keyLabel}</span>
        <span className="text-xs text-ink-faint" aria-hidden>▾</span>
      </div>
      <p className="mt-1.5 text-[11px] text-ink-faint">Read from the sample names.</p>
    </>
  );
}
