import { Link } from "react-router";
import { MOD_KEY } from "../lib/modKey";
import { Kbd } from "./Kbd";
import { Tooltip } from "./Tooltip";
import { Wordmark } from "./Wordmark";

interface Props {
  experimentName: string;
  onActions: () => void;
  onHelp: () => void;
}

export function WorkspaceHeader({ experimentName, onActions, onHelp }: Props): JSX.Element {
  return (
    <>
      <span className="flex min-w-0 items-baseline gap-4">
        <Link
          to="/experiments"
          aria-label="Himalaya SAXS, go to experiments"
          className="shrink-0 rounded no-underline focus-visible:outline focus-visible:outline-2 focus-visible:outline-offset-2 focus-visible:outline-accent"
        >
          <Wordmark tail="SAXS">Himalaya</Wordmark>
        </Link>
        <nav aria-label="Breadcrumb" className="flex min-w-0 items-baseline gap-1.5 text-base">
          <Link to="/experiments" className="shrink-0 text-ink-soft no-underline hover:text-ink">
            Experiments
          </Link>
          <span aria-hidden="true" className="text-ink-faint">/</span>
          <Tooltip content={experimentName} overflowOnly>
            <h1 aria-current="page" className="m-0 min-w-0 truncate text-base font-semibold text-ink">{experimentName}</h1>
          </Tooltip>
        </nav>
      </span>
      <span className="ml-auto" aria-hidden="true" />
      <button
        type="button"
        onClick={onActions}
        aria-label={`Actions, ${MOD_KEY}K`}
        className="flex items-center gap-2 rounded-sm border border-hair bg-plate py-1 pl-2.5 pr-1 text-base text-ink-soft hover:border-hair-strong hover:text-ink"
      >
        Actions
        <span className="flex gap-0.5">
          <Kbd>{MOD_KEY}</Kbd>
          <Kbd>K</Kbd>
        </span>
      </button>
      <Tooltip content="Keyboard shortcuts">
        <button
          type="button"
          onClick={onHelp}
          aria-label="Keyboard shortcuts"
          className="flex items-center rounded-sm border border-hair bg-plate px-1 py-1 hover:border-hair-strong"
        >
          <Kbd>?</Kbd>
        </button>
      </Tooltip>
    </>
  );
}
