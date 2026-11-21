import { useEffect, useLayoutEffect, useRef, useState } from "react";
import type { SampleLabel, TagColor, TagColorEntry } from "../../api";
import { defaultTagColor, displayTag, sameLabel, standardTag, tagColor, tagColorIsDefault } from "../../lib/labels";
import { TagColorMenu } from "./TagColorMenu";
import { TagPicker } from "./TagPicker";
import { TAG_SWATCH } from "./tagSwatches";

interface Props {
  tags: readonly SampleLabel[];
  vocabulary: readonly string[];
  editable: boolean;
  adding: boolean;
  onOpenAdd: () => void;
  onCloseAdd: () => void;
  onAdd: (label: string) => void;
  onRemove: (label: string) => void;
  colors: readonly TagColorEntry[];
  coloring: string | undefined;
  onOpenColor: (label: string) => void;
  onCloseColor: () => void;
  onSetColor: (label: string, color: TagColor | null) => void;
}

const GAP = 6;
const MORE_WIDTH = 44;

interface ChipProps {
  tag: SampleLabel;
  editable: boolean;
  color: TagColor;
  measure?: boolean;
  onAdd?: (label: string) => void;
  onRemove?: (label: string) => void;
  onOpenColor?: (label: string) => void;
}

function Chip({ tag, editable, color, measure = false, onAdd, onRemove, onOpenColor }: ChipProps): JSX.Element {
  const suggested = tag.state === "suggested";
  const standard = standardTag(tag.label) !== undefined;
  const name = displayTag(tag.label);
  const tone = TAG_SWATCH[color].ink;
  const fill = suggested ? "border-dashed border-current bg-transparent" : `border-transparent ${TAG_SWATCH[color].wash}`;
  return (
    <span
      data-tag={measure ? undefined : tag.label}
      data-tag-state={measure ? undefined : tag.state}
      data-tag-standard={standard && !measure ? "true" : undefined}
      data-tag-color={measure ? undefined : color}
      onContextMenu={
        editable && onOpenColor
          ? (e) => {
              e.preventDefault();
              onOpenColor(tag.label);
            }
          : undefined
      }
      className={`inline-flex shrink-0 items-center gap-0.5 whitespace-nowrap rounded-full border pt-0.5 pb-[3px] pl-2.5 pr-1 text-caption ${tone} ${fill}`}
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

export function TagRow({
  tags,
  vocabulary,
  editable,
  adding,
  onOpenAdd,
  onCloseAdd,
  onAdd,
  onRemove,
  colors,
  coloring,
  onOpenColor,
  onCloseColor,
  onSetColor,
}: Props): JSX.Element {
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
  const colored = coloring === undefined ? undefined : tags.find((t) => sameLabel(t.label, coloring));
  const [menuLeft, setMenuLeft] = useState(14);
  const opener = useRef<HTMLElement | null>(null);
  const restoreFocus = useRef(false);
  const addButton = useRef<HTMLButtonElement | null>(null);
  const chipFor = (label: string): HTMLElement | null =>
    box.current?.querySelector<HTMLElement>(`[data-tag="${CSS.escape(label)}"]`) ?? null;
  useLayoutEffect(() => {
    if (!colored) return;
    const at = chipFor(colored.label) ?? moreButton.current;
    setMenuLeft(at ? at.offsetLeft : 14);
  }, [colored, shown]);
  useLayoutEffect(() => {
    if (coloring === undefined) return;
    const active = document.activeElement;
    opener.current = active instanceof HTMLElement && active !== document.body ? active : null;
  }, [coloring]);
  useEffect(() => {
    if (coloring !== undefined) return;
    const back = opener.current;
    opener.current = null;
    if (!restoreFocus.current) return;
    restoreFocus.current = false;
    const safe = back && back.isConnected && !back.closest("[data-tag]") ? back : null;
    (safe ?? addButton.current ?? moreButton.current)?.focus();
  }, [coloring]);
  const closeColor = (restore: boolean): void => {
    restoreFocus.current = restore;
    onCloseColor();
  };
  const openColor = (label: string): void => {
    setMoreOpen(false);
    onOpenColor(label);
  };
  const chipProps = (t: SampleLabel) => ({
    tag: t,
    editable,
    color: tagColor(colors, t.label),
    onAdd,
    onRemove,
    onOpenColor: openColor,
  });
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
          <Chip key={t.id} {...chipProps(t)} />
        ))}
        {hidden.length > 0 && (
          <button
            ref={moreButton}
            type="button"
            aria-expanded={moreOpen}
            aria-label={`${hidden.length} more ${hidden.length === 1 ? "tag" : "tags"}`}
            onClick={() => setMoreOpen((o) => !o)}
            data-testid="workspace-tags-more"
            className="inline-flex shrink-0 items-center whitespace-nowrap rounded-full border border-hair-strong px-2.5 pt-0.5 pb-[3px] text-caption text-ink-soft hover:text-ink"
          >
            +{hidden.length}
          </button>
        )}
      </div>
      <div aria-hidden="true" className="pointer-events-none invisible absolute inset-0 overflow-hidden">
        <div ref={ruler} className="absolute left-0 top-0 flex w-max gap-1.5" data-testid="workspace-tag-ruler">
          {tags.map((t) => (
            <Chip key={t.id} tag={t} editable={editable} color={tagColor(colors, t.label)} measure />
          ))}
        </div>
      </div>
      {editable && (
        <button
          type="button"
          ref={addButton}
          onClick={onOpenAdd}
          className="inline-flex shrink-0 items-center whitespace-nowrap rounded-full border border-dashed border-hair-strong px-2.5 pt-0.5 pb-[3px] text-caption text-ink-soft hover:text-ink"
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
            <Chip key={t.id} {...chipProps(t)} />
          ))}
        </div>
      )}
      {colored && editable && (
        <TagColorMenu
          name={displayTag(colored.label)}
          current={tagColor(colors, colored.label)}
          isDefault={tagColorIsDefault(colors, colored.label)}
          left={menuLeft}
          onPick={(c) => onSetColor(colored.label, c === defaultTagColor(colored.label) ? null : c)}
          onClose={closeColor}
        />
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
