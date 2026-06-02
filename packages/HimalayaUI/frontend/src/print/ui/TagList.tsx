import { useState } from "react";
import { TagPill } from "./TagPill";
import { Chip } from "./Chip";
import type { ChipSize } from "./Chip";
import { TagEditor } from "./TagEditor";
import { Popover } from "./Popover";
import type { Tag } from "./tag";

interface TagListProps {
  tags: Tag[];
  onAdd?: (tag: Tag) => void;
  onRemove?: (tag: Tag) => void;
  editable?: boolean;
  /** Size axis, forwarded to every pill AND the add invite. Defaults to `"sm"`. */
  size?: ChipSize;
  /** When set and `tags.length > maxVisible`, only the first `maxVisible` pills
   *  render, followed by a single muted "+N more" trigger that opens a popover
   *  revealing the hidden tags as real {@link TagPill}s. Unset → all tags render
   *  (the unbounded, wrapping case). */
  maxVisible?: number;
  className?: string;
}

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

/** A wrapping list of key-value {@link TagPill}s plus an optional add affordance.
 *  Composes TagPill, the base {@link Chip} `add` variant, and {@link TagEditor}.
 *
 *  B — the rationed accent appears ONLY on the add invite's hover/focus (Chip's
 *  `add` variant), never on the resting list or pills. G — when tags exist the
 *  invite is hidden until row hover via a tokenized OPACITY transition
 *  (`group-hover/tags` + `focus-within`), never a layout/size animation; when the
 *  list is EMPTY it is always visible (the empty-state invite). */
export function TagList({
  tags,
  onAdd,
  onRemove,
  editable = false,
  size = "sm",
  maxVisible,
  className = "",
}: TagListProps): JSX.Element {
  const [adding, setAdding] = useState(false);

  const capped = maxVisible !== undefined && tags.length > maxVisible;
  const visible = capped ? tags.slice(0, maxVisible) : tags;
  const hidden = capped ? tags.slice(maxVisible) : [];

  return (
    <div
      data-testid="tag-list"
      className={cx("flex items-center flex-wrap gap-1.5 group/tags", className)}
    >
      {visible.map((t, i) => (
        <TagPill
          key={i}
          tag={t}
          size={size}
          {...(editable && onRemove ? { onRemove: () => onRemove(t) } : {})}
        />
      ))}
      {capped && (
        <Popover
          label="More tags"
          trigger={
            <button
              type="button"
              data-testid="tag-overflow"
              className="text-ink-faint text-xs hover:text-ink cursor-pointer focus-visible:outline focus-visible:outline-1 focus-visible:outline-accent rounded-sm"
            >
              +{hidden.length} more
            </button>
          }
        >
          <div className="flex flex-wrap gap-1.5">
            {hidden.map((t, i) => (
              <TagPill
                key={i}
                tag={t}
                size={size}
                {...(editable && onRemove ? { onRemove: () => onRemove(t) } : {})}
              />
            ))}
          </div>
        </Popover>
      )}
      {onAdd &&
        (adding ? (
          <TagEditor
            onCommit={(tag) => {
              onAdd(tag);
              setAdding(false);
            }}
            onCancel={() => setAdding(false)}
          />
        ) : (
          <Chip
            variant="add"
            size={size}
            onClick={() => setAdding(true)}
            className={cx(
              tags.length > 0 &&
                "opacity-0 group-hover/tags:opacity-100 focus-within:opacity-100",
            )}
          >
            + tag
          </Chip>
        ))}
    </div>
  );
}
