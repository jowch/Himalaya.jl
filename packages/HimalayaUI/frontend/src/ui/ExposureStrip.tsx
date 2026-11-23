import type { Exposure } from "../api";
import { imageUrl } from "../workspace/model";
import { detectorLutFilter } from "../detector/detectorLut";
import { Tooltip } from "./Tooltip";

const DETECTOR_COLOUR = { filter: detectorLutFilter("magma") };

export interface ExposureStripGroup {
  key: number;
  title: string | undefined;
  exposures: readonly Exposure[];
}

interface Props {
  groups: readonly ExposureStripGroup[];
  currentId: number | undefined;
  onlyId?: number | undefined;
  canShowOnly?: boolean;
  onPick: (id: number) => void;
  onToggleDrop: (e: Exposure) => void;
  onToggleOnly?: (id: number) => void;
}

const ACTION = "flex-1 rounded-[4px] border py-0.5 text-xs ";
const ACTION_OFF = "border-hair bg-paper text-ink hover:border-hair-strong";
const ACTION_ON = "border-ink bg-ink text-paper";

export function ExposureStrip({ groups, currentId, onlyId, canShowOnly = false, onPick, onToggleDrop, onToggleOnly }: Props): JSX.Element {
  return (
    <div className="flex flex-col gap-1.5" data-testid="workspace-exposures">
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
              const src = imageUrl(e, true);
              const name = gi === 0 || g.title === undefined ? `Image ${i + 1}` : `Image ${i + 1} of ${g.title}`;
              return (
                <li
                  key={e.id}
                  className={
                    "group flex min-w-0 flex-col gap-1 rounded-sm border p-1 " +
                    (cur ? "border-hair-strong " : "border-transparent ") +
                    (dropped ? "opacity-50" : "")
                  }
                >
                  <button
                    type="button"
                    aria-current={cur ? "true" : undefined}
                    aria-label={`${name}${dropped ? ", dropped" : ""}`}
                    onClick={() => onPick(e.id)}
                    className="relative grid aspect-square w-full place-items-center overflow-hidden rounded-[3px] bg-frame-edge"
                  >
                    {src ? (
                      <img src={src} alt="" loading="lazy" style={DETECTOR_COLOUR} className="h-full w-full object-cover" />
                    ) : (
                      <span className="text-xs text-frame-tag">No TIFF</span>
                    )}
                  </button>
                  {dropped && <span className="text-xs text-ink-soft line-through">dropped</span>}
                  <div className={(cur ? "visible" : "invisible group-hover:visible group-focus-within:visible") + " flex gap-1"}>
                    <button
                      type="button"
                      onClick={() => onToggleDrop(e)}
                      className={ACTION + ACTION_OFF}
                    >
                      {dropped ? "Restore" : "Drop"}
                    </button>
                    {canShowOnly && onToggleOnly && (
                      <button
                        type="button"
                        aria-pressed={e.id === onlyId}
                        onClick={() => onToggleOnly(e.id)}
                        className={ACTION + (e.id === onlyId ? ACTION_ON : ACTION_OFF)}
                      >
                        Show only
                      </button>
                    )}
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
