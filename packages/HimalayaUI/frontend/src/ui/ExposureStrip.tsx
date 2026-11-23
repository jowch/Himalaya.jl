import type { MouseEvent } from "react";
import type { Exposure } from "../api";
import { imageUrl } from "../workspace/model";
import { detectorLutFilter } from "../detector/detectorLut";
import { Checkbox } from "./Checkbox";
import { Tooltip } from "./Tooltip";

const DETECTOR_COLOUR = { filter: detectorLutFilter("magma") };

export interface ExposureStripGroup {
  key: number;
  title: string | undefined;
  exposures: readonly Exposure[];
}

export interface ExposureSelection {
  selectedIds: ReadonlySet<number>;
  selectableIds: ReadonlySet<number>;
  onToggle: (id: number) => void;
  onRange: (id: number) => void;
  onAll: () => void;
  onNone: () => void;
}

interface Props {
  groups: readonly ExposureStripGroup[];
  currentId: number | undefined;
  selection?: ExposureSelection | undefined;
  onPick: (id: number) => void;
  onToggleDrop: (e: Exposure) => void;
}

const ACTION = "flex-1 rounded-[4px] border border-hair bg-paper py-0.5 text-xs text-ink hover:border-hair-strong";
const TOOL = "rounded-[4px] border border-hair bg-paper px-2 py-0.5 text-xs text-ink hover:border-hair-strong disabled:opacity-40 disabled:hover:border-hair";

function selectionCaption(selection: ExposureSelection): string {
  const n = selection.selectedIds.size;
  if (n === 0) return "No images selected";
  if (n === 1) return "1 image selected";
  return `Average of ${n} images`;
}

export function ExposureStrip({ groups, currentId, selection, onPick, onToggleDrop }: Props): JSX.Element {
  const click = (ev: MouseEvent, id: number): void => {
    const canSelect = selection?.selectableIds.has(id) ?? false;
    if (canSelect && ev.shiftKey) selection?.onRange(id);
    else if (canSelect && (ev.metaKey || ev.ctrlKey)) selection?.onToggle(id);
    else onPick(id);
  };
  const all = selection ? [...selection.selectableIds].every((id) => selection.selectedIds.has(id)) : true;
  return (
    <div className="flex flex-col gap-1.5" data-testid="workspace-exposures">
      {selection && (
        <div className="flex items-center gap-1.5" data-testid="exposure-selection-tools">
          <span className="mr-auto text-xs text-ink-soft" data-testid="exposure-selection-caption">{selectionCaption(selection)}</span>
          <button type="button" className={TOOL} disabled={all} onClick={selection.onAll}>Select all</button>
          <button type="button" className={TOOL} disabled={selection.selectedIds.size === 0} onClick={selection.onNone}>Select none</button>
        </div>
      )}
      {groups.map((g, gi) => (
        <section key={g.key} aria-label={g.title} className="flex flex-col gap-1">
          {g.title !== undefined && (
            <Tooltip content={g.title} overflowOnly>
              <h4 data-testid="frame-group-header" className="m-0 truncate text-xs font-semibold text-ink-soft">{g.title}</h4>
            </Tooltip>
          )}
          <ul className="m-0 grid list-none grid-cols-6 gap-2 p-0">
            {g.exposures.map((e, i) => {
              const cur = e.id === currentId;
              const dropped = e.status === "rejected";
              const selectable = selection?.selectableIds.has(e.id) ?? false;
              const selected = selection?.selectedIds.has(e.id) ?? false;
              const src = imageUrl(e, true);
              const name = gi === 0 || g.title === undefined ? `Image ${i + 1}` : `Image ${i + 1} of ${g.title}`;
              return (
                <li
                  key={e.id}
                  data-selected={selected ? "true" : undefined}
                  className={
                    "group flex min-w-0 flex-col gap-1 rounded-sm p-1 " +
                    (cur ? "bg-row-cursor shadow-[inset_3px_0_0_var(--color-accent)] " : "") +
                    (dropped ? "opacity-50" : "")
                  }
                >
                  <div className="relative">
                    <button
                      type="button"
                      aria-current={cur ? "true" : undefined}
                      aria-label={`${name}${dropped ? ", dropped" : ""}`}
                      onClick={(ev) => click(ev, e.id)}
                      className="grid aspect-square w-full place-items-center overflow-hidden rounded-[3px] bg-frame-edge"
                    >
                      {src ? (
                        <img src={src} alt="" loading="lazy" style={DETECTOR_COLOUR} className="h-full w-full object-cover" />
                      ) : (
                        <span className="text-xs text-frame-tag">No TIFF</span>
                      )}
                    </button>
                    {selectable && selection && (
                      <span className="absolute left-0.5 top-0.5 rounded-full bg-paper">
                        <Checkbox
                          checked={selected}
                          aria-label={`Analyze ${name.toLowerCase()}`}
                          onChange={() => selection.onToggle(e.id)}
                          onClickCapture={(ev) => {
                            if (ev.shiftKey) {
                              ev.stopPropagation();
                              selection.onRange(e.id);
                            }
                          }}
                        />
                      </span>
                    )}
                  </div>
                  {dropped && <span className="text-xs text-ink-soft line-through">dropped</span>}
                  <div className={(cur ? "visible" : "invisible group-hover:visible group-focus-within:visible") + " flex gap-1"}>
                    <button type="button" onClick={() => onToggleDrop(e)} className={ACTION}>
                      {dropped ? "Restore" : "Drop"}
                    </button>
                  </div>
                </li>
              );
            })}
          </ul>
        </section>
      ))}
    </div>
  );
}
