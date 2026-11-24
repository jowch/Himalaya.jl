import type { MouseEvent } from "react";
import type { Exposure } from "../api";
import { imageUrl } from "../workspace/model";
import { detectorLutFilter } from "../detector/detectorLut";
import { Button } from "./Button";
import { Checkbox } from "./Checkbox";
import { Tooltip } from "./Tooltip";
import { releaseMouseFocus, useRovingFocus } from "../hooks/useRovingFocus";

const DETECTOR_COLOUR = { filter: detectorLutFilter("magma") };

export interface ExposureStripGroup {
  key: number;
  title: string | undefined;
  exposures: readonly Exposure[];
  selection?: ExposureSelection | undefined;
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
  onPick: (id: number) => void;
  onToggleDrop: (e: Exposure) => void;
  warningOf?: (e: Exposure) => string | undefined;
}

const ACTION = "shrink-0 rounded-[4px] px-2 border border-hair bg-paper py-0.5 text-xs text-ink hover:border-hair-strong";

function selectionCaption(selection: ExposureSelection): string {
  const n = selection.selectedIds.size;
  if (n === 0) return "No images selected";
  if (n === 1) return "1 image selected";
  return `Average of ${n} images`;
}

function SelectionTools({ selection, of }: { selection: ExposureSelection; of: string | undefined }): JSX.Element {
  const all = [...selection.selectableIds].every((id) => selection.selectedIds.has(id));
  const suffix = of === undefined ? "" : ` of ${of}`;
  return (
    <div className="flex items-center gap-1.5" data-testid="exposure-selection-tools">
      <span className="mr-auto text-xs text-ink-soft" data-testid="exposure-selection-caption">{selectionCaption(selection)}</span>
      <Button type="button" variant="outline" disabled={all} aria-label={`Select all${suffix}`} onClick={selection.onAll}>Select all</Button>
      <Button type="button" variant="outline" disabled={selection.selectedIds.size === 0} aria-label={`Select none${suffix}`} onClick={selection.onNone}>Select none</Button>
    </div>
  );
}

export function ExposureStrip({ groups, currentId, onPick, onToggleDrop, warningOf }: Props): JSX.Element {
  const ref = useRovingFocus<HTMLDivElement>(currentId);
  const hasCurrent = groups.some((g) => g.exposures.some((e) => e.id === currentId));
  const click = (ev: MouseEvent<HTMLElement>, id: number, selection: ExposureSelection | undefined): void => {
    releaseMouseFocus(ev);
    const canSelect = selection?.selectableIds.has(id) ?? false;
    if (canSelect && ev.shiftKey) selection?.onRange(id);
    else if (canSelect && (ev.metaKey || ev.ctrlKey)) selection?.onToggle(id);
    else onPick(id);
  };
  return (
    <div ref={ref} className="flex flex-col gap-1.5" data-testid="workspace-exposures">
      {groups.map((g, gi) => {
        const selection = g.selection;
        return (
        <section key={g.key} aria-label={g.title} className="flex flex-col gap-1">
          {g.title !== undefined && (
            <Tooltip content={g.title} overflowOnly>
              <h4 data-testid="frame-group-header" className="m-0 truncate text-xs font-semibold text-ink-soft">{g.title}</h4>
            </Tooltip>
          )}
          {selection && <SelectionTools selection={selection} of={gi === 0 ? undefined : g.title} />}
          <ul className="-mx-1 my-0 grid list-none grid-cols-[repeat(auto-fill,136px)] gap-2 p-0">
            {g.exposures.map((e, i) => {
              const cur = e.id === currentId;
              const stop = cur || (!hasCurrent && gi === 0 && i === 0) ? 0 : -1;
              const dropped = e.status === "rejected";
              const selectable = selection?.selectableIds.has(e.id) ?? false;
              const selected = selection?.selectedIds.has(e.id) ?? false;
              const src = imageUrl(e, true);
              const warning = warningOf?.(e);
              const of = gi === 0 || g.title === undefined ? "" : ` of ${g.title}`;
              const name = `Image ${i + 1}${of}`;
              return (
                <li
                  key={e.id}
                  data-selected={selected ? "true" : undefined}
                  data-warned={warning !== undefined ? "" : undefined}
                  className={
                    "group flex min-w-0 flex-col gap-1 rounded-sm p-1 " +
                    (cur ? "bg-row-cursor " : "") +
                    (dropped ? "opacity-50" : "")
                  }
                >
                  <div className="relative">
                    <button
                      type="button"
                      aria-current={cur ? "true" : undefined}
                      aria-label={`${name}${dropped ? ", dropped" : ""}${warning !== undefined ? `, ${warning.toLowerCase()}` : ""}`}
                      tabIndex={stop}
                      data-roving="image"
                      onClick={(ev) => click(ev, e.id, selection)}
                      className={
                        "relative grid aspect-square w-full place-items-center overflow-hidden rounded-[3px] bg-frame-edge" +
                        (warning !== undefined ? " outline outline-2 -outline-offset-2 outline-warning" : "")
                      }
                    >
                      {src ? (
                        <img src={src} alt="" loading="lazy" style={DETECTOR_COLOUR} className="h-full w-full object-cover" />
                      ) : (
                        <span className="text-xs text-frame-tag">No TIFF</span>
                      )}
                    </button>
                    {selectable && selection && (
                      <span
                        data-testid="exposure-select-plate"
                        onClick={(ev) => {
                          if (ev.target !== ev.currentTarget) return;
                          if (ev.shiftKey) selection.onRange(e.id);
                          else selection.onToggle(e.id);
                        }}
                        className="absolute left-1 top-1 flex h-6 w-6 cursor-pointer items-center justify-center rounded-sm border border-hair bg-paper"
                      >
                        <Checkbox
                          checked={selected}
                          tabIndex={stop}
                          data-roving="select"
                          aria-label={`Analyze image ${i + 1}${of}`}
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
                  <div data-tile-footer="" className="flex h-6 min-w-0 items-center gap-1">
                    <span className={"min-w-0 flex-1 truncate text-xs " + (warning !== undefined ? "text-warning" : "text-ink-soft")} aria-hidden={dropped ? undefined : "true"}>
                      {dropped ? "Dropped" : warning ?? `Image ${i + 1}`}
                    </span>
                    <button
                      type="button"
                      tabIndex={stop}
                      data-roving="drop"
                      onClick={(ev) => {
                        releaseMouseFocus(ev);
                        onToggleDrop(e);
                      }}
                      className={(cur ? "visible " : "invisible group-hover:visible group-focus-within:visible ") + ACTION}
                    >
                      {dropped ? "Restore" : "Drop"}
                    </button>
                  </div>
                </li>
              );
            })}
          </ul>
        </section>
        );
      })}
    </div>
  );
}
