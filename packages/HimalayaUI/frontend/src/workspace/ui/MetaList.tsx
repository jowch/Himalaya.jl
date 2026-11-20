export function MetaList({ rows }: { rows: readonly [string, string][] }): JSX.Element {
  return (
    <dl className="grid grid-cols-[auto_1fr] gap-x-3 gap-y-1 text-body">
      {rows.map(([k, v]) => (
        <div key={k} className="contents">
          <dt className="text-ink-soft">{k}</dt>
          <dd className="truncate text-right" title={v}>{v}</dd>
        </div>
      ))}
    </dl>
  );
}
