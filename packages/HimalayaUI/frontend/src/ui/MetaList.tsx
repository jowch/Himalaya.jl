import { Term } from "./Term";
import { Tooltip } from "./Tooltip";

interface Props {
  rows: readonly [string, string][];
  explain?: Readonly<Record<string, string>>;
}

export function MetaList({ rows, explain = {} }: Props): JSX.Element {
  return (
    <dl className="grid grid-cols-[auto_minmax(0,1fr)] gap-x-3 gap-y-[3px] text-xs">
      {rows.map(([k, v]) => (
        <div key={k} className="contents">
          <dt className="text-ink-faint">{explain[k] ? <Term explain={explain[k]}>{k}</Term> : k}</dt>
          <Tooltip content={v} overflowOnly>
            <dd className="m-0 min-w-0 truncate text-right text-base">{v}</dd>
          </Tooltip>
        </div>
      ))}
    </dl>
  );
}
