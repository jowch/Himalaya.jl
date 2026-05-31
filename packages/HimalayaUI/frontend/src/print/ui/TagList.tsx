import { TagPill } from "./TagPill";

interface TagListProps {
  tags: string[];
  onAdd?: () => void;
  onRemove?: (tag: string) => void;
  editable?: boolean;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** A wrapping list of {@link TagPill}s plus an optional dashed "+ tag" add invite
 *  (mockup `.tags` + `.tag-add`). Composes TagPill (DRY).
 *
 *  B — the rationed accent appears ONLY on the add invite's hover/focus
 *  (→accent), never on the resting list or pills. C — the add button ships
 *  default / hover / focus-visible; the named-group reveal is keyboard-reachable
 *  because `focus-visible:opacity-100` un-hides it on focus even before row hover.
 *  When tags exist the invite is hidden until row hover (`group-hover/tags`); when
 *  the list is EMPTY it is always visible (the empty-state invite). G — only a
 *  tokenized OPACITY/color transition, never a layout/size animation. */
export function TagList({
  tags,
  onAdd,
  onRemove,
  editable = false,
  className = "",
}: TagListProps): JSX.Element {
  return (
    <div
      data-testid="tag-list"
      className={cx("flex items-center flex-wrap gap-1.5 group/tags", className)}
    >
      {tags.map((t) => (
        <TagPill
          key={t}
          {...(editable && onRemove ? { onRemove: () => onRemove(t) } : {})}
        >
          {t}
        </TagPill>
      ))}
      {onAdd && (
        <button
          type="button"
          data-testid="tag-add"
          aria-label="Add tag"
          onClick={onAdd}
          className={cx(
            "text-xs font-semibold rounded-full border border-dashed border-hair-strong px-2 py-px text-ink-faint transition-colors",
            "hover:text-accent hover:border-accent focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent",
            tags.length > 0 &&
              "opacity-0 group-hover/tags:opacity-100 focus-visible:opacity-100",
          )}
        >
          + tag
        </button>
      )}
    </div>
  );
}
