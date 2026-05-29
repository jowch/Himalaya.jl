interface Props {
  keyLabel: string;
}

/**
 * The "Ordered by" field (series-scoping.html `.order-field`): a read-out of
 * the ordering variable Himalaya inferred from the sample names. M2 dropped the
 * decorative `▾` chevron that made it present as an openable dropdown — there
 * is no picker behind it. Re-selecting the ordering variable is a real future
 * feature (it ripples through the proposal pipeline); until then this is an
 * honest static read-out, not a control.
 */
export function ScopingOrderField({ keyLabel }: Props): JSX.Element {
  return (
    <>
      <div className="mb-1.5 mt-5 text-[10.5px] font-bold uppercase tracking-wider text-ink-faint">
        Ordered by
      </div>
      <div
        data-testid="scoping-order-field"
        className="flex items-center rounded-md border border-hair-strong bg-plate px-3.5 py-3"
      >
        <span className="text-[15px] font-semibold text-ink">{keyLabel}</span>
      </div>
      <p className="mt-1.5 text-[11px] text-ink-faint">Read from the sample names.</p>
    </>
  );
}
