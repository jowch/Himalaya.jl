import { useEffect, useLayoutEffect, useRef, useState } from "react";
import type { SampleLabel, LabelColor, LabelColorEntry } from "../api";
import { defaultLabelColor, displayLabel, sameLabel, standardLabel, labelColor, labelColorIsDefault } from "../lib/labels";
import { Kbd } from "./Kbd";
import { LabelColorMenu } from "./LabelColorMenu";
import { LabelPicker } from "./LabelPicker";
import { LABEL_SWATCH } from "./labelSwatches";

interface Props {
  labels: readonly SampleLabel[];
  vocabulary: readonly string[];
  editable: boolean;
  adding: boolean;
  onOpenAdd: () => void;
  onCloseAdd: () => void;
  onAdd: (label: string) => void;
  onRemove: (label: string) => void;
  colors: readonly LabelColorEntry[];
  coloring: string | undefined;
  onOpenColor: (label: string) => void;
  onCloseColor: () => void;
  onSetColor: (label: string, color: LabelColor | null) => void;
}

const GAP = 6;
const MORE_WIDTH = 44;

interface ChipProps {
  entry: SampleLabel;
  editable: boolean;
  color: LabelColor;
  measure?: boolean;
  onAdd?: (label: string) => void;
  onRemove?: (label: string) => void;
  onOpenColor?: (label: string) => void;
}

function Chip({ entry, editable, color, measure = false, onAdd, onRemove, onOpenColor }: ChipProps): JSX.Element {
  const suggested = entry.state === "suggested";
  const standard = standardLabel(entry.label) !== undefined;
  const name = displayLabel(entry.label);
  const tone = LABEL_SWATCH[color].ink;
  const fill = suggested ? "border-dashed border-current bg-transparent" : `border-transparent ${LABEL_SWATCH[color].wash}`;
  return (
    <span
      data-label={measure ? undefined : entry.label}
      data-label-state={measure ? undefined : entry.state}
      data-label-standard={standard && !measure ? "true" : undefined}
      data-label-color={measure ? undefined : color}
      onContextMenu={
        editable && onOpenColor && !standard
          ? (e) => {
              e.preventDefault();
              onOpenColor(entry.label);
            }
          : undefined
      }
      className={`inline-flex shrink-0 items-center gap-0.5 whitespace-nowrap rounded-full border pt-0.5 pb-[3px] pl-2.5 pr-1 text-caption ${tone} ${fill}`}
    >
      {suggested && editable && onAdd ? (
        <button type="button" onClick={() => onAdd(entry.label)} title="Suggested. Click to confirm" aria-label={`Confirm label ${name}`}>
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
          aria-label={`Remove label ${name}`}
          onClick={() => onRemove?.(entry.label)}
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

export function LabelRow({
  labels,
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
  const last = labels.filter((t) => t.state === "confirmed").at(-1);
  const box = useRef<HTMLDivElement | null>(null);
  const ruler = useRef<HTMLDivElement | null>(null);
  const moreRoot = useRef<HTMLDivElement | null>(null);
  const moreButton = useRef<HTMLButtonElement | null>(null);
  const [shown, setShown] = useState(labels.length);
  const [moreOpen, setMoreOpen] = useState(false);

  useLayoutEffect(() => {
    const measure = (): void => {
      const room = box.current?.clientWidth ?? 0;
      const chips = ruler.current ? Array.from(ruler.current.children) : [];
      const widths = chips.map((c) => (c as HTMLElement).offsetWidth);
      setShown(room > 0 && widths.some((w) => w > 0) ? fitCount(widths, room) : labels.length);
    };
    measure();
    if (typeof ResizeObserver === "undefined" || !box.current) return;
    const ro = new ResizeObserver(measure);
    ro.observe(box.current);
    if (ruler.current) ro.observe(ruler.current);
    return () => ro.disconnect();
  }, [labels, editable]);

  const hidden = labels.slice(shown);
  const colored = coloring === undefined ? undefined : labels.find((t) => sameLabel(t.label, coloring));
  const [menuLeft, setMenuLeft] = useState(14);
  const opener = useRef<HTMLElement | null>(null);
  const restoreFocus = useRef(false);
  const addButton = useRef<HTMLButtonElement | null>(null);
  const chipFor = (label: string): HTMLElement | null =>
    box.current?.querySelector<HTMLElement>(`[data-label="${CSS.escape(label)}"]`) ?? null;
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
    const safe = back && back.isConnected && !back.closest("[data-label]") ? back : null;
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
    entry: t,
    editable,
    color: labelColor(colors, t.label),
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
    <div className="relative flex h-6 w-full min-w-0 flex-nowrap items-center justify-end gap-1.5" data-testid="workspace-labels">
      <div ref={box} className={"flex w-0 flex-1 flex-nowrap items-center justify-end gap-1.5 overflow-hidden " + (labels.length > 0 ? "min-w-11" : "min-w-0")} data-testid="workspace-label-chips">
        {labels.slice(0, shown).map((t) => (
          <Chip key={t.id} {...chipProps(t)} />
        ))}
        {hidden.length > 0 && (
          <button
            ref={moreButton}
            type="button"
            aria-expanded={moreOpen}
            aria-label={`${hidden.length} more ${hidden.length === 1 ? "label" : "labels"}`}
            onClick={() => setMoreOpen((o) => !o)}
            data-testid="workspace-labels-more"
            className="inline-flex shrink-0 items-center whitespace-nowrap rounded-full border border-hair-strong px-2.5 pt-0.5 pb-[3px] text-caption text-ink-soft hover:text-ink"
          >
            +{hidden.length}
          </button>
        )}
      </div>
      <div aria-hidden="true" className="pointer-events-none invisible absolute inset-0 overflow-hidden">
        <div ref={ruler} className="absolute left-0 top-0 flex w-max gap-1.5" data-testid="workspace-label-ruler">
          {labels.map((t) => (
            <Chip key={t.id} entry={t} editable={editable} color={labelColor(colors, t.label)} measure />
          ))}
        </div>
      </div>
      {editable && (
        <button
          type="button"
          ref={addButton}
          onClick={onOpenAdd}
          className="inline-flex shrink-0 items-center gap-1 whitespace-nowrap rounded-full border border-dashed border-hair-strong pl-2.5 pr-1 text-caption text-ink-soft hover:text-ink"
        >
          + <span className="@max-[34rem]:sr-only">Label</span> <Kbd>L</Kbd>
        </button>
      )}
      {moreOpen && hidden.length > 0 && (
        <div
          ref={moreRoot}
          role="dialog"
          aria-label="More labels"
          data-keys-trap
          data-testid="workspace-labels-overflow"
          className="absolute right-0 top-full z-20 mt-1 flex max-w-80 flex-wrap gap-1.5 rounded-sm border border-hair bg-plate p-2 shadow-lg"
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
        <LabelColorMenu
          name={displayLabel(colored.label)}
          current={labelColor(colors, colored.label)}
          isDefault={labelColorIsDefault(colors, colored.label)}
          left={menuLeft}
          onPick={(c) => onSetColor(colored.label, c === defaultLabelColor(colored.label) ? null : c)}
          onClose={closeColor}
        />
      )}
      {adding && editable && (
        <LabelPicker
          vocabulary={vocabulary}
          onSample={labels.filter((t) => t.state === "confirmed").map((t) => t.label)}
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
