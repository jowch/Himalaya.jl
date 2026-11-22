import type { Exposure } from "../api";
import { imageUrl } from "../workspace/model";
import { detectorLutFilter } from "../detector/detectorLut";

const DETECTOR_COLOUR = { filter: detectorLutFilter("magma") };

interface Props {
  exposures: readonly Exposure[];
  currentId: number | undefined;
  onPick: (id: number) => void;
  onToggleDrop: (e: Exposure) => void;
}

export function ExposureStrip({ exposures, currentId, onPick, onToggleDrop }: Props): JSX.Element {
  return (
    <ul className="m-0 grid list-none grid-cols-6 gap-2 p-0" data-testid="workspace-exposures">
      {exposures.map((e, i) => {
        const cur = e.id === currentId;
        const dropped = e.status === "rejected";
        const src = imageUrl(e, true);
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
              aria-label={`Frame ${i + 1}${dropped ? ", dropped" : ""}`}
              onClick={() => onPick(e.id)}
              className="relative grid aspect-square w-full place-items-center overflow-hidden rounded-[3px] bg-frame-edge"
            >
              {src ? (
                <img src={src} alt="" loading="lazy" style={DETECTOR_COLOUR} className="h-full w-full object-cover" />
              ) : (
                <span className="text-xs text-frame-tag">No image</span>
              )}
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
  );
}
