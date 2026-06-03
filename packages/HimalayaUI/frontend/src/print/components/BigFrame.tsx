import type { ReactNode } from "react";
import { DetectorImage } from "../detector";
import type { DetectorLutVariant } from "../detector/detectorLut";
import { RejectOverlay } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface BigFrameProps {
  /** Detector image URL; null → DetectorImage frame-window placeholder. */
  src: string | null;
  /** Mono caption set over the dark frame (the `.frame-tag`), e.g. "frame 65 · 0.40 s". */
  caption?: ReactNode;
  /** true → dims the image, shows the "Dropped" pill + the grease-pencil ✕. */
  rejected?: boolean;
  /** Forwarded to DetectorImage (colormap). */
  lutVariant?: DetectorLutVariant;
  /** PLACEMENT ONLY. */
  className?: string;
}

export function BigFrame({
  src,
  caption,
  rejected,
  lutVariant,
  className,
}: BigFrameProps): JSX.Element {
  return (
    <div
      data-testid="big-frame"
      data-rejected={rejected ? "true" : undefined}
      className={cx(
        "relative aspect-square rounded overflow-hidden bg-frame-edge border border-frame-edge max-w-[500px] mx-auto",
        className,
      )}
    >
      <DetectorImage
        src={src}
        size="full"
        className={cx("h-full w-full", rejected && "opacity-40")}
        {...(lutVariant !== undefined ? { lutVariant } : {})}
      />
      {rejected && (
        <span
          data-role="dropped-tag"
          className="absolute left-3 top-3 bg-accent text-plate uppercase font-bold tracking-wide text-xs rounded-sm px-2 py-[3px]"
        >
          Dropped
        </span>
      )}
      {rejected && <RejectOverlay className="absolute inset-0" />}
      {caption != null && (
        <span
          data-role="frame-caption"
          className="absolute bottom-2.5 left-3 font-mono text-frame-tag text-sm tracking-wide"
        >
          {caption}
        </span>
      )}
    </div>
  );
}
