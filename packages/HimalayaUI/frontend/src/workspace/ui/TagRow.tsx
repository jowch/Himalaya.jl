import { useEffect, useLayoutEffect, useRef, useState } from "react";
import type { SampleLabel } from "../../api";
import { displayTag, standardTag } from "../../lib/labels";
import { TagPicker } from "./TagPicker";

interface Props {
  tags: readonly SampleLabel[];
  vocabulary: readonly string[];
  editable: boolean;
  adding: boolean;
  onOpenAdd: () => void;
  onCloseAdd: () => void;
  onAdd: (label: string) => void;
  onRemove: (label: string) => void;
}

const GAP = 6;
const MORE_WIDTH = 44;

interface ChipProps {
  tag: SampleLabel;
  editable: boolean;
  measure?: boolean;
  onAdd?: (label: string) => void;
  onRemove?: (label: string) => void;
}

function Chip({ tag, editable, measure = false, onAdd, onRemove }: ChipProps): JSX.Element {
  const suggested = tag.state === "suggested";
  const standard = standardTag(tag.label) !== undefined;
  const name = displayTag(tag.label);
  const tone = standard ? "text-tag-std-ink" : "text-tag-ink";
  const fill = suggested
    ? "border-dashed border-current bg-transparent"
    : standard
      ? "border-transparent bg-tag-std-wash"
      : "border-transparent bg-tag-wash";
  return (
    <span
      data-tag={measure ? undefined : tag.label}
      data-tag-state={measure ? undefined : tag.state}
      data-tag-standard={standard && !measure ? "true" : undefined}
      className={`inline-flex shrink-0 items-center gap-0.5 whitespace-nowrap rounded-full border pb-0 pl-2.5 pr-1 pt-1 text-caption ${tone} ${fill}`}
    >
      {suggested && editable && onAdd ? (
        <button type="button" onClick={() => onAdd(tag.label)} title="Suggested. Click to confirm" aria-label={`Confirm tag ${name}`}>
          {name}
        </button>
      ) : (
        <span>
          {name}
          {suggested && <span className="sr-only">, suggested</span>}
        </span>
      )}
      {editable && (
        <button
          type="button"
          aria-label={`Remove tag ${name}`}
          onClick={() => onRemove?.(tag.label)}
          className="rounded-full px-1 opacity-70 hover:opacity-100"
        >
          ×
        </button>
      )}
    </span>
  );
}

function fitCount(widths: readonly number[], room: number): number {
  let used = 0;
  for (let i = 0; i < widths.length; i++) {
    const next = used + (i > 0 ? GAP : 0) + widths[i]!;
    const reserve = i < widths.length - 1 ? GAP + MORE_WIDTH : 0;
    if (next + reserve > room) return i;
    used = next;
  }
  return widths.length;
}

export function TagRow({ tags, vocabulary, editable, adding, onOpenAdd, onCloseAdd, onAdd, onRemove }: Props): JSX.Element {
  const last = tags.filter((t) => t.state === "confirmed").at(-1);
  const box = useRef<HTMLDivElement | null>(null);
  const ruler = useRef<HTMLDivElement | null>(null);
  const moreRoot = useRef<HTMLDivElement | null>(null);
  const moreButton = useRef<HTMLButtonElement | null>(null);
  const [shown, setShown] = useState(tags.length);
  const [moreOpen, setMoreOpen] = useState(false);

  useLayoutEffect(() => {
    const measure = (): void => {
      const room = box.current?.clientWidth ?? 0;
      const chips = ruler.current ? Array.from(ruler.current.children) : [];
      const widths = chips.map((c) => (c as HTMLElement).offsetWidth);
      setShown(room > 0 && widths.some((w) => w > 0) ? fitCount(widths, room) : tags.length);
    };
    measure();
    if (typeof ResizeObserver === "undefined" || !box.current) return;
    const ro = new ResizeObserver(measure);
    ro.observe(box.current);
    if (ruler.current) ro.observe(ruler.current);
    return () => ro.disconnect();
  }, [tags, editable]);

  const hidden = tags.slice(shown);
  useEffect(() => {
    if (hidden.length === 0) setMoreOpen(false);
  }, [hidden.length]);
  useEffect(() => {
    if (moreOpen) moreRoot.current?.querySelector<HTMLElement>("button")?.focus();
  }, [moreOpen]);
  useEffect(() => {
    if (!moreOpen) return;
    const onDown = (ev: MouseEvent): void => {
      if (moreRoot.current && ev.target instanceof Node && !moreRoot.current.contains(ev.target)) setMoreOpen(false);
    };
    document.addEventListener("mousedown", onDown, true);
    return () => document.removeEventListener("mousedown", onDown, true);
  }, [moreOpen]);

  return (
    <div className="relative flex h-7 flex-nowrap items-center gap-1.5 px-3.5" data-testid="workspace-tags">
      <div ref={box} className="flex min-w-0 flex-1 flex-nowrap items-center gap-1.5 overflow-hidden" data-testid="workspace-tag-chips">
        {tags.slice(0, shown).map((t) => (
          <Chip key={t.id} tag={t} editable={editable} onAdd={onAdd} onRemove={onRemove} />
        ))}
        {hidden.length > 0 && (
          <button
            ref={moreButton}
            type="button"
            aria-expanded={moreOpen}
            aria-label={`${hidden.length} more ${hidden.length === 1 ? "tag" : "tags"}`}
            onClick={() => setMoreOpen((o) => !o)}
            data-testid="workspace-tags-more"
            className="inline-flex shrink-0 items-center whitespace-nowrap rounded-full border border-hair-strong px-2.5 pb-0 pt-1 text-caption text-ink-soft hover:text-ink"
          >
            +{hidden.length}
          </button>
        )}
      </div>
      <div aria-hidden="true" className="pointer-events-none invisible absolute inset-0 overflow-hidden">
        <div ref={ruler} className="absolute left-0 top-0 flex w-max gap-1.5" data-testid="workspace-tag-ruler">
          {tags.map((t) => (
            <Chip key={t.id} tag={t} editable={editable} measure />
          ))}
        </div>
      </div>
      {editable && (
        <button
          type="button"
          onClick={onOpenAdd}
          className="inline-flex shrink-0 items-center whitespace-nowrap rounded-full border border-dashed border-hair-strong px-2.5 pb-0 pt-1 text-caption text-ink-soft hover:text-ink"
        >
          + Tag
        </button>
      )}
      {moreOpen && hidden.length > 0 && (
        <div
          ref={moreRoot}
          role="dialog"
          aria-label="More tags"
          data-keys-trap
          data-testid="workspace-tags-overflow"
          className="absolute right-3.5 top-full z-20 mt-1 flex max-w-80 flex-wrap gap-1.5 rounded-sm border border-hair bg-plate p-2 shadow-lg"
          onKeyDown={(e) => {
            if (e.key === "Escape") {
              e.preventDefault();
              e.stopPropagation();
              setMoreOpen(false);
              moreButton.current?.focus();
            }
          }}
        >
          {hidden.map((t) => (
            <Chip key={t.id} tag={t} editable={editable} onAdd={onAdd} onRemove={onRemove} />
          ))}
        </div>
      )}
      {adding && editable && (
        <TagPicker
          vocabulary={vocabulary}
          onSample={tags.filter((t) => t.state === "confirmed").map((t) => t.label)}
          onPick={onAdd}
          onRemoveLast={() => {
            if (last) onRemove(last.label);
          }}
          onClose={onCloseAdd}
        />
      )}
    </div>
  );
}
