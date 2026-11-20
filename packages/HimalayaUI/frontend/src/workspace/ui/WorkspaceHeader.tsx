import { Link } from "react-router-dom";

interface Props {
  experimentId: number;
  experimentName: string;
  reviewed: number;
  total: number;
}

export function WorkspaceHeader({ experimentId, experimentName, reviewed, total }: Props): JSX.Element {
  return (
    <>
      <Link to={`/experiments/${experimentId}/corpus`} className="text-label text-ink-soft no-underline">
        ← Corpus
      </Link>
      <h1 className="text-title truncate">{experimentName}</h1>
      <span className="text-caption" data-testid="workspace-progress">
        {reviewed} of {total} reviewed
      </span>
    </>
  );
}
