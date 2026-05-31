import { useState } from "react";
import { TagPill } from "./TagPill";
import { Chip } from "./Chip";
import { TagEditor } from "./TagEditor";
import type { Tag } from "./tag";

interface TagListProps {
  tags: Tag[];
  onAdd?: (tag: Tag) => void;
  onRemove?: (tag: Tag) => void;
  editable?: boolean;
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
  className = "",
}: TagListProps): JSX.Element {
  const [adding, setAdding] = useState(false);

  return (
    <div
      data-testid="tag-list"
      className={cx("flex items-center flex-wrap gap-1.5 group/tags", className)}
    >
      {tags.map((t, i) => (
        <TagPill
          key={i}
          tag={t}
          {...(editable && onRemove ? { onRemove: () => onRemove(t) } : {})}
        />
      ))}
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
