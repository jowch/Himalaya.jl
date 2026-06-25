import type { JSX } from "react";
import type { LoadExposure } from "../../api";
import { Thumbnail } from "./Thumbnail";
import { IconButton } from "../ui/IconButton";

export interface ExposureLeafProps {
  exposure: LoadExposure;
  /** Detector thumbnail URL for this exposure, or null (Thumbnail draws a
   *  placeholder). The page derives it; Thumbnail takes `src`, NOT an id. */
  thumbSrc: string | null;
  /** Opens the Move... picker for this exposure.
   *  `anchorEl` is the trigger button — the page uses it to position the menu. */
  onMove: (exposureId: number, anchorEl: HTMLElement) => void;
  className?: string;
}

/** One exposure row in the grouping fold: thumb | filename | H | time | Move.
 *  Appearance lives in the composed primitives; this is layout only. */
export function ExposureLeaf({ exposure, thumbSrc, onMove, className }: ExposureLeafProps): JSX.Element {
  const h = exposure.horizontal_position;
  return (
    <div
      data-testid="exposure-leaf"
      className={`flex items-center gap-3 rounded-sm px-2.5 py-1.5${className ? ` ${className}` : ""}`}
    >
      <Thumbnail src={thumbSrc} size="xs" frameNo={exposure.id} />
      <span className="min-w-0 flex-1 truncate font-mono text-xs text-ink">{exposure.filename}</span>
      <span className="w-28 shrink-0 font-mono text-xs text-ink-soft">
        {h === null ? "—" : `H ${h.toFixed(1)} mm`}
      </span>
      <span className="w-20 shrink-0 font-mono text-xs text-ink-faint">{exposure.timestamp ?? "—"}</span>
      <IconButton label="Move to another sample" onClick={(e) => onMove(exposure.id, e.currentTarget)}>&#x22EF;</IconButton>
    </div>
  );
}
