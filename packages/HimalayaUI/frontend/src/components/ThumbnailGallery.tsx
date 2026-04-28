import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";

interface Props {
  exposures: Exposure[];
  selectedId: number | undefined;
  onSelect: (id: number) => void;
  className?: string;
}

/*
 * ThumbnailGallery — horizontal filmstrip.
 *
 * Cells fill the container height; width derived from aspect-ratio (3/4).
 * Overflow scrolls horizontally. The gallery is always rendered inside a
 * bounded-height slot (e.g. ~140px in InspectPage), so a wrapping grid mode
 * is no longer needed — the page-level WorkspaceGrid handles layout.
 */
export function ThumbnailGallery({
  exposures,
  selectedId,
  onSelect,
  className,
}: Props): JSX.Element {
  const containerClass = "flex flex-row gap-2 overflow-x-auto h-full";

  return (
    <div className={`${containerClass} ${className ?? ""}`}>
      {exposures.map((e) => {
        const isViewing  = e.id === selectedId;
        const isRejected = e.status === "rejected";
        const isIndexing = e.selected;

        /* border-radius: thumb is rounded-md (8px). chip offset is 5px.
           chip border-radius = 8 - 5 = 3px → rounded-[3px].
           Keep in sync if rounded-md or offset changes. */
        return (
          <div
            key={e.id}
            data-testid={`thumb-cell-${e.id}`}
            data-rejected={isRejected ? "true" : undefined}
            onClick={() => onSelect(e.id)}
            className={[
              "relative flex flex-col items-center gap-1 cursor-pointer group",
              // h-full fills strip; aspect-[3/4] derives width from height
              "aspect-[3/4] h-full",
              isRejected ? "opacity-40 grayscale" : "",
            ].join(" ")}
          >
            <div
              className={[
                "relative w-full overflow-hidden rounded-md transition-all duration-150",
                // flex-1 fills cell height minus label
                "flex-1 min-h-0",
                // Unified ring widths (always ring-2) prevent layout shift on hover.
                // No ring-offset → ring sits flush against the rounded corner.
                isViewing
                  ? "ring-2 ring-accent shadow-sm shadow-accent/30"
                  : "ring-1 ring-border/50 group-hover:ring-2 group-hover:ring-fg-muted/40",
              ].join(" ")}
            >
              <DetectorImage
                exposureId={e.id}
                imagePath={e.image_path}
                imageVersion={e.image_version}
                size="thumb"
                className="w-full h-full"
              />
              {isIndexing && (
                <span
                  className="absolute top-[5px] left-[5px] rounded-[3px]
                             bg-accent/85 text-bg text-xs font-semibold
                             px-1.5 py-0.5 leading-snug backdrop-blur-sm"
                >
                  ⊙ Indexing
                </span>
              )}
            </div>
            <span className="text-caption truncate w-full text-center">
              {e.filename?.replace(/\.dat$/i, "") ?? `#${e.id}`}
            </span>
          </div>
        );
      })}
    </div>
  );
}
