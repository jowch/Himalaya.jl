import { Link } from "react-router-dom";
import { InlineEditableText } from "./InlineEditableText";
import { relativeTime } from "../lib/comparison/relativeTime";

interface Props {
  title: string;
  description: string | null;
  memberCount: number;
  authorUsername: string | null;
  isCurrentUserAuthor: boolean;
  lastEventAt: string | null;
  forkedFromTitle: string | null;
  forkedFromHref: string | null;
  onTitleChange: (s: string) => void;
  onDescChange: (s: string) => void;
  /**
   * When true, title and description render as plain text (no click-to-edit
   * affordance). Used in review mode so an idle click doesn't silently
   * discard the user's typing. Default `false`.
   */
  readOnly?: boolean;
  now?: number; // injected for tests
}

export function CompareTitleStrip(p: Props): JSX.Element {
  const byline = p.isCurrentUserAuthor
    ? "by you"
    : p.authorUsername !== null
      ? `by ${p.authorUsername}`
      : "by —";
  const rel = relativeTime(p.lastEventAt, p.now);
  return (
    <div data-testid="compare-title-strip" className="flex flex-col gap-1">
      <h2 className="text-lg font-semibold text-fg">
        <InlineEditableText
          value={p.title}
          onCommit={p.onTitleChange}
          placeholder="Untitled comparison"
          testId="compare-title"
          readOnly={p.readOnly ?? false}
        />
      </h2>
      <div className="text-sm text-fg-dim flex items-center gap-1 flex-wrap">
        <span>{byline}</span>
        <span aria-hidden="true">·</span>
        <span>{rel === null ? "—" : `edited ${rel}`}</span>
        <span aria-hidden="true">·</span>
        <span>{p.memberCount} traces</span>
        {p.forkedFromHref !== null && p.forkedFromTitle !== null && (
          <>
            <span aria-hidden="true">·</span>
            <span>
              forked from{" "}
              <Link to={p.forkedFromHref} className="text-accent hover:underline">
                {p.forkedFromTitle}
              </Link>
            </span>
          </>
        )}
      </div>
      {p.description !== null && (
        <div className="text-sm text-fg-dim">
          <InlineEditableText
            value={p.description ?? ""}
            onCommit={p.onDescChange}
            placeholder="Add a description…"
            testId="compare-description"
            readOnly={p.readOnly ?? false}
          />
        </div>
      )}
    </div>
  );
}
