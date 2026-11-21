import type { Exposure } from "../../api";
import { imageUrl } from "../model";
import { Kbd } from "./Kbd";
import { Overlay } from "./Overlay";

const CLOSE_KEYS = ["Escape", " "] as const;

interface Props {
  exposures: readonly Exposure[];
  currentId: number;
  onStep: (dir: 1 | -1) => void;
  onToggleDrop: (() => void) | undefined;
  onClose: () => void;
}

export function DetectorOverlay({ exposures, currentId, onStep, onToggleDrop, onClose }: Props): JSX.Element {
  const index = exposures.findIndex((e) => e.id === currentId);
  const current = exposures[index];
  const dropped = current?.status === "rejected";
  const onKey = (key: string, repeat: boolean): boolean => {
    if (key === "ArrowRight" || key === "ArrowLeft") {
      onStep(key === "ArrowRight" ? 1 : -1);
      return true;
    }
    if ((key === "x" || key === "X") && onToggleDrop) {
      if (!repeat) onToggleDrop();
      return true;
    }
    return false;
  };
  return (
    <Overlay label="Detector image" onClose={onClose} testId="workspace-detector-big" closeKeys={CLOSE_KEYS} onKey={onKey}>
      <div className="flex flex-col items-center gap-2 p-3">
        <div className="relative aspect-square w-[min(420px,80vw)] overflow-hidden rounded-[4px] bg-frame-edge">
          {exposures.map((e) => {
            const src = imageUrl(e, false);
            if (!src) return null;
            const shown = e.id === currentId;
            return (
              <img
                key={e.id}
                src={src}
                alt=""
                data-exposure-id={e.id}
                aria-hidden={shown ? undefined : true}
                className={"absolute inset-0 h-full w-full object-contain " + (shown ? (dropped ? "opacity-40" : "") : "invisible")}
              />
            );
          })}
          {current && !imageUrl(current, false) && (
            <p className="absolute inset-0 m-0 flex items-center justify-center text-sm text-ink-faint">No detector image</p>
          )}
        </div>
        <p className="m-0 flex w-full items-baseline justify-between gap-3 text-sm text-ink-soft" data-testid="workspace-detector-caption">
          <span>
            Frame {index + 1} of {exposures.length}
            <span className="text-ink-faint"> · {current?.filename ?? `Exposure ${currentId}`}</span>
          </span>
          {dropped && <span className="font-semibold text-ink line-through">dropped</span>}
        </p>
        <p className="m-0 flex flex-wrap items-center gap-x-3 gap-y-1 text-xs text-ink-faint">
          <span><Kbd>←</Kbd> <Kbd>→</Kbd> frames</span>
          {onToggleDrop && <span><Kbd>X</Kbd> {dropped ? "restore" : "drop"}</span>}
          <span><Kbd>Esc</Kbd> close</span>
        </p>
      </div>
    </Overlay>
  );
}
