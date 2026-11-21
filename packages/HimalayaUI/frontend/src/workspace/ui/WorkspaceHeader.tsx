import { Link } from "react-router-dom";
import { Kbd } from "./Kbd";

interface Props {
  experimentId: number;
  experimentName: string;
  reviewed: number;
  total: number;
  onHelp: () => void;
}

export function WorkspaceHeader({ experimentId, experimentName, reviewed, total, onHelp }: Props): JSX.Element {
  const pct = total > 0 ? (reviewed / total) * 100 : 0;
  return (
    <>
      <Link to={`/experiments/${experimentId}/corpus`} title="Back to the corpus" className="text-lg font-semibold leading-none text-ink no-underline">
        Himalaya
      </Link>
      <h1 className="min-w-0 truncate text-base font-normal text-ink-soft">{experimentName}</h1>
      <span className="ml-auto flex items-center gap-2 whitespace-nowrap text-base text-ink-soft" data-testid="workspace-progress">
        <span className="block h-[5px] w-[90px] overflow-hidden rounded-sm bg-paper-sunk" aria-hidden="true">
          <i className="block h-full bg-success" style={{ width: `${pct}%` }} />
        </span>
        {reviewed} of {total} reviewed
      </span>
      <button
        type="button"
        onClick={onHelp}
        aria-label="Keyboard shortcuts"
        title="Keyboard shortcuts"
        className="flex items-center rounded-sm border border-hair bg-paper px-2 py-1 hover:border-hair-strong"
      >
        <Kbd>?</Kbd>
      </button>
    </>
  );
}
