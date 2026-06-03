import { DetectorImage } from "../detector/DetectorImage";
import { RejectOverlay } from "../ui/RejectOverlay";
import { Dot } from "../ui/Dot";

export interface ThumbnailProps {
  /** Detector thumbnail source URL; null → DetectorImage placeholder. */
  src: string | null;
  /** Frame number label, e.g. "65". */
  frameNo?: string | number;
  representative?: boolean;
  rejected?: boolean;
  selected?: boolean;
  /** `"xs"` → 30px (dense Focus exposure strip); `"sm"` → 62px (contact-sheet strip); `"lg"` → 70px (loupe strip). */
  size?: "xs" | "sm" | "lg";
  onClick?: () => void;
  title?: string;
  /** PLACEMENT ONLY. */
  className?: string;
}

const SIZE_PX: Record<"xs" | "sm" | "lg", number> = {
  xs: 30,
  sm: 62,
  lg: 70,
};

function buildDataState(props: {
  selected?: boolean;
  rejected?: boolean;
  representative?: boolean;
}): string {
  const tokens: string[] = [];
  if (props.selected) tokens.push("selected");
  if (props.rejected) tokens.push("rejected");
  if (props.representative) tokens.push("representative");
  if (tokens.length === 0) tokens.push("normal");
  return tokens.join(" ");
}

export function Thumbnail({
  src,
  frameNo,
  representative = false,
  rejected = false,
  selected = false,
  size = "sm",
  onClick,
  title,
  className,
}: ThumbnailProps): JSX.Element {
  const px = SIZE_PX[size];
  const dataState = buildDataState({ selected, rejected, representative });

  return (
    <button
      data-testid="thumbnail"
      data-state={dataState}
      data-size={size}
      onClick={onClick}
      title={title}
      style={{ width: px, height: px }}
      // Selection/hover use an INSET ring, not an outset one: the gallery wraps
      // thumbnails in `overflow-x-auto`, which clips an outset box-shadow ring to
      // the scroll content box (it showed on only one or two sides, and not at
      // all on edge thumbnails). An inset ring is painted inside the element's
      // own box, so it always renders the full perimeter on every thumbnail.
      className={`relative inline-block flex-shrink-0 overflow-hidden rounded-sm bg-frame-edge border border-frame-edge p-0 cursor-pointer${selected ? " inset-ring-[3px] inset-ring-accent" : " hover:inset-ring-2 hover:inset-ring-hair-strong"}${className ? ` ${className}` : ""}`}
    >
      {/* Detector image — dimmed when rejected */}
      <div
        data-dimmed={rejected ? "true" : "false"}
        className={rejected ? "w-full h-full opacity-40" : "w-full h-full"}
      >
        <DetectorImage src={src} size="thumb" className="w-full h-full" />
      </div>

      {/* Frame-number label */}
      {frameNo != null && (
        <span
          data-role="thumb-fno"
          className={`absolute left-1 bottom-0.5 font-mono text-frame-tag text-xs${rejected ? " opacity-50" : " opacity-80"}`}
        >
          {frameNo}
        </span>
      )}

      {/* Representative marker */}
      {representative && (
        <span data-role="thumb-rep" className="absolute top-0.5 right-0.5 flex">
          <Dot tone="accent" size="md" bordered aria-hidden="true" />
        </span>
      )}

      {/* Reject overlay */}
      {rejected && <RejectOverlay />}
    </button>
  );
}
