import { Chip } from "./Chip";
import type { ChipSize } from "./Chip";
import type { Tag } from "./tag";

interface TagPillProps {
  tag: Tag;
  onRemove?: () => void;
  /** Accessible label for the remove button (default `"Remove"`). Pass a
   *  per-tag string such as `"Remove dose 10"` to make duplicate pills
   *  individually named for assistive technology — see LO-TAGDUP. */
  removeLabel?: string;
  /** Size axis, forwarded to the base Chip. Defaults to `"sm"` — tags are dense. */
  size?: ChipSize;
  className?: string;
}

/** A single key-value tag pill, composed on the base {@link Chip}. A key+value
 *  tag reads as KEY (faint sans label) + VALUE (mono ink measurement); a
 *  key-only tag is just the bare key.
 *
 *  E — type carries the semantics: the value is a measurement, so it is set in
 *  `font-mono text-ink`; the key is a label, so it is faint sans (`text-ink-soft`).
 *  B — the remove × is Chip's NEUTRAL affordance (no accent), reached only when
 *  `onRemove` is supplied (variant flips static → removable). */
export function TagPill({
  tag,
  onRemove,
  removeLabel,
  size = "sm",
  className = "",
}: TagPillProps): JSX.Element {
  const content =
    tag.value !== undefined ? (
      <>
        <span data-role="tag-key" className="text-ink-soft">
          {tag.key}
        </span>
        <span data-role="tag-value" className="font-mono text-ink ml-1">
          {tag.value}
        </span>
      </>
    ) : (
      <span data-role="tag-key">{tag.key}</span>
    );

  return (
    <Chip
      variant={onRemove ? "removable" : "static"}
      size={size}
      testId="tag-pill"
      {...(onRemove ? { onRemove } : {})}
      {...(removeLabel !== undefined ? { removeLabel } : {})}
      className={className}
    >
      {content}
    </Chip>
  );
}
