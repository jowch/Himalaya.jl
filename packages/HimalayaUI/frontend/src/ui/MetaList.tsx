export function MetaList({ rows }: { rows: readonly [string, string][] }): JSX.Element {
  return (
    <dl className="grid grid-cols-[auto_minmax(0,1fr)] gap-x-3 gap-y-[3px] text-xs">
      {rows.map(([k, v]) => (
        <div key={k} className="contents">
          <dt className="text-ink-faint">{k}</dt>
          <dd className="m-0 min-w-0 truncate text-right text-base" title={v}>{v}</dd>
        </div>
      ))}
    </dl>
  );
}
