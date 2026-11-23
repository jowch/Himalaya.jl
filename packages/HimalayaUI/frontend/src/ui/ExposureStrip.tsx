import type { MouseEvent } from "react";
import type { Exposure } from "../api";
import { imageUrl } from "../workspace/model";
import { detectorLutFilter } from "../detector/detectorLut";
import { Button } from "./Button";
import { Checkbox } from "./Checkbox";
import { Tooltip } from "./Tooltip";
import { DetectorRings, type DetectorRingsView } from "./DetectorRings";
import { releaseMouseFocus, useRovingFocus } from "../hooks/useRovingFocus";
import { useTileColumns } from "./useTileColumns";
import { FULL_IMAGE_TILE } from "./tileColumns";

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
  warningOf?: (e: Exposure) => string | undefined;
  ringsOf?: ((e: Exposure) => DetectorRingsView | null) | undefined;
}

const GRID_COLUMNS: Record<number, string> = {
  2: "grid-cols-2",
  3: "grid-cols-3",
  4: "grid-cols-4",
  6: "grid-cols-6",
};

const ACTION = "flex-1 rounded-[4px] border border-hair bg-paper py-0.5 text-xs text-ink hover:border-hair-strong";

function selectionCaption(selection: ExposureSelection): string {
  const n = selection.selectedIds.size;
  if (n === 0) return "No images selected";
  if (n === 1) return "1 image selected";
  return `Average of ${n} images`;
}

export function ExposureStrip({ groups, currentId, selection, onPick, onToggleDrop, warningOf, ringsOf }: Props): JSX.Element {
  const ref = useRovingFocus<HTMLDivElement>(currentId);
  const fit = useTileColumns(ref);
  const fullImages = fit.tile > FULL_IMAGE_TILE;
  const hasCurrent = groups.some((g) => g.exposures.some((e) => e.id === currentId));
  const click = (ev: MouseEvent<HTMLElement>, id: number): void => {
    releaseMouseFocus(ev);
    const canSelect = selection?.selectableIds.has(id) ?? false;
    if (canSelect && ev.shiftKey) selection?.onRange(id);
    else if (canSelect && (ev.metaKey || ev.ctrlKey)) selection?.onToggle(id);
    else onPick(id);
  };
  const all = selection ? [...selection.selectableIds].every((id) => selection.selectedIds.has(id)) : true;
  return (
    <div ref={ref} className="flex flex-col gap-1.5" data-testid="workspace-exposures">
      {selection && (
        <div className="flex items-center gap-1.5" data-testid="exposure-selection-tools">
          <span className="mr-auto text-xs text-ink-soft" data-testid="exposure-selection-caption">{selectionCaption(selection)}</span>
          <Button type="button" variant="outline" disabled={all} onClick={selection.onAll}>Select all</Button>
          <Button type="button" variant="outline" disabled={selection.selectedIds.size === 0} onClick={selection.onNone}>Select none</Button>
        </div>
      )}
      {groups.map((g, gi) => (
        <section key={g.key} aria-label={g.title} className="flex flex-col gap-1">
          {g.title !== undefined && (
            <Tooltip content={g.title} overflowOnly>
              <h4 data-testid="frame-group-header" className="m-0 truncate text-xs font-semibold text-ink-soft">{g.title}</h4>
            </Tooltip>
          )}
          <ul data-columns={fit.columns} className={"m-0 grid list-none gap-2 p-0 " + (GRID_COLUMNS[fit.columns] ?? "grid-cols-6")}>
            {g.exposures.map((e, i) => {
              const cur = e.id === currentId;
              const stop = cur || (!hasCurrent && gi === 0 && i === 0) ? 0 : -1;
              const dropped = e.status === "rejected";
              const selectable = selection?.selectableIds.has(e.id) ?? false;
              const selected = selection?.selectedIds.has(e.id) ?? false;
              const src = imageUrl(e, !fullImages);
              const warning = warningOf?.(e);
              const ringView = src ? ringsOf?.(e) ?? null : null;
              const name = gi === 0 || g.title === undefined ? `Image ${i + 1}` : `Image ${i + 1} of ${g.title}`;
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
                      onClick={(ev) => click(ev, e.id)}
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
                      {ringView && <DetectorRings view={ringView} fit="cover" />}
                    </button>
                    {selectable && selection && (
                      <span className="absolute left-0.5 top-0.5 rounded-full bg-paper">
                        <Checkbox
                          checked={selected}
                          tabIndex={stop}
                          data-roving="select"
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
                  {warning !== undefined && <span aria-hidden="true" className="truncate text-xs text-warning">{warning}</span>}
                  <div className={(cur ? "visible" : "invisible group-hover:visible group-focus-within:visible") + " flex gap-1"}>
                    <button
                      type="button"
                      tabIndex={stop}
                      data-roving="drop"
                      onClick={(ev) => {
                        releaseMouseFocus(ev);
                        onToggleDrop(e);
                      }}
                      className={ACTION}
                    >
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
