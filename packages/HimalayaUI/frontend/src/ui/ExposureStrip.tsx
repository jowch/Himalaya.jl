import type { Exposure } from "../api";
import { imageUrl } from "../workspace/model";
import { detectorLutFilter } from "../detector/detectorLut";
import { Tooltip } from "./Tooltip";
import { DetectorRings, type DetectorRingsView } from "./DetectorRings";

const DETECTOR_COLOUR = { filter: detectorLutFilter("magma") };

export interface ExposureStripGroup {
  key: number;
  title: string | undefined;
  exposures: readonly Exposure[];
}

interface Props {
  groups: readonly ExposureStripGroup[];
  currentId: number | undefined;
  onPick: (id: number) => void;
  onToggleDrop: (e: Exposure) => void;
  ringsOf?: ((e: Exposure) => DetectorRingsView | null) | undefined;
}

export function ExposureStrip({ groups, currentId, onPick, onToggleDrop, ringsOf }: Props): JSX.Element {
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
              const ringView = src ? ringsOf?.(e) ?? null : null;
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
                    {ringView && <DetectorRings view={ringView} fit="cover" />}
                  </button>
                  {dropped && <span className="text-xs text-ink-soft line-through">dropped</span>}
                  <div className={(cur ? "visible" : "invisible group-hover:visible group-focus-within:visible") + " flex gap-1"}>
                    <button
                      type="button"
                      onClick={() => onToggleDrop(e)}
                      className="flex-1 rounded-[4px] border border-hair bg-paper py-0.5 text-xs text-ink hover:border-hair-strong"
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
