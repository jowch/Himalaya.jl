import { Link } from "react-router";
import { Kbd } from "./Kbd";
import { Tooltip } from "./Tooltip";

interface Props {
  experimentName: string;
  reviewed: number;
  total: number;
  onHelp: () => void;
}

export function WorkspaceHeader({ experimentName, reviewed, total, onHelp }: Props): JSX.Element {
  const pct = (reviewed / Math.max(total, 1)) * 100;
  return (
    <>
      <span className="flex min-w-0 items-baseline gap-3.5">
        <Tooltip content="All experiments">
          <Link to="/experiments" className="text-lg font-semibold text-ink no-underline">
            Himalaya
          </Link>
        </Tooltip>
        <h1 className="min-w-0 truncate text-base font-normal text-ink-soft">{experimentName}</h1>
      </span>
      <span className="ml-auto" aria-hidden="true" />
      <span
        className={"flex items-center gap-2 whitespace-nowrap text-base text-ink-soft" + (total > 0 ? "" : " invisible")}
        data-testid={total > 0 ? "workspace-progress" : undefined}
        data-reserved={total > 0 ? undefined : ""}
        aria-hidden={total > 0 ? undefined : true}
      >
        <span className="block h-[5px] w-[90px] overflow-hidden rounded-sm bg-paper-sunk" aria-hidden="true">
          <i className="block h-full bg-success" style={{ width: `${pct}%` }} />
        </span>
        {reviewed} of {total} reviewed
      </span>
      <Tooltip content="Keyboard shortcuts">
        <button
          type="button"
          onClick={onHelp}
          aria-label="Keyboard shortcuts"
          className="flex items-center rounded-sm border border-hair bg-paper px-2 py-1 hover:border-hair-strong"
        >
          <Kbd>?</Kbd>
        </button>
      </Tooltip>
    </>
  );
}
