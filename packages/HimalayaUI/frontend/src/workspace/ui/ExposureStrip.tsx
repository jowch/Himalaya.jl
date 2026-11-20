import type { Exposure } from "../../api";
import { imageUrl } from "../model";
import { Kbd } from "./Kbd";

interface Props {
  exposures: readonly Exposure[];
  currentId: number | undefined;
  onPick: (id: number) => void;
  onToggleDrop: (e: Exposure) => void;
}

export function ExposureStrip({ exposures, currentId, onPick, onToggleDrop }: Props): JSX.Element {
  return (
    <ul className="grid grid-cols-[repeat(auto-fill,minmax(112px,1fr))] gap-2" data-testid="workspace-exposures">
      {exposures.map((e, i) => {
        const cur = e.id === currentId;
        const dropped = e.status === "rejected";
        const src = imageUrl(e, true);
        return (
          <li key={e.id} className={"flex flex-col gap-1 rounded-sm p-1 " + (cur ? "outline outline-2 outline-accent" : "")}>
            <button
              type="button"
              aria-current={cur ? "true" : undefined}
              aria-label={`Frame ${i + 1}${dropped ? ", dropped" : ""}`}
              onClick={() => onPick(e.id)}
              className={"relative aspect-square w-full overflow-hidden rounded-sm bg-frame-edge " + (dropped ? "opacity-40" : "")}
            >
              {src ? (
                <img src={src} alt="" loading="lazy" className="h-full w-full object-cover" />
              ) : (
                <span className="text-xs text-frame-tag">No image</span>
              )}
            </button>
            {cur && (
              <div className="flex items-center justify-between gap-1 text-xs">
                <span className="text-ink-soft">{dropped ? "Dropped" : ""}</span>
                <button
                  type="button"
                  onClick={() => onToggleDrop(e)}
                  className="flex items-center gap-1 rounded-sm border border-hair-strong bg-plate px-1.5 text-ink hover:bg-paper-sunk"
                >
                  {dropped ? "Restore" : "Drop"} <Kbd>X</Kbd>
                </button>
              </div>
            )}
          </li>
        );
      })}
    </ul>
  );
}
