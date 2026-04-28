import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";

interface Props {
  exposures: Exposure[];
  selectedId: number | undefined;
  onSelect: (id: number) => void;
  className?: string;
}

/*
 * ThumbnailGallery — responsive layout via CSS only:
 *
 *   < 1400px  — thin horizontal filmstrip (flex-row).
 *               Cells fill the container height; width derived from aspect-ratio.
 *               No vertical scroll needed.
 *
 *   ≥ 1400px  — auto-fill wrapping grid (minmax 90px columns).
 *               Height derived from column width via aspect-ratio.
 *               Overflows vertically so all exposures are reachable.
 */
export function ThumbnailGallery({
  exposures,
  selectedId,
  onSelect,
  className,
}: Props): JSX.Element {
  const containerClass = [
    // Base (< 1400px): horizontal filmstrip
    "flex flex-row gap-2 overflow-x-auto h-full",
    // Large (≥ 1400px): wrapping auto-fill grid
    "min-[1400px]:grid min-[1400px]:flex-none",
    "min-[1400px]:grid-cols-[repeat(auto-fill,minmax(90px,1fr))]",
    "min-[1400px]:overflow-x-hidden min-[1400px]:overflow-y-auto",
    "min-[1400px]:h-auto min-[1400px]:content-start",
  ].join(" ");

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
            onClick={() => onSelect(e.id)}
            className={[
              "relative flex flex-col items-center gap-1 cursor-pointer group",
              // filmstrip: h-full fills strip; aspect-[3/4] derives width from height
              // grid (≥1400px): h-auto; width from grid column, height from aspect-ratio
              "aspect-[3/4] h-full min-[1400px]:h-auto",
              isRejected ? "opacity-40 grayscale" : "",
            ].join(" ")}
          >
            <div
              className={[
                "relative w-full overflow-hidden rounded-md transition-all duration-150",
                // filmstrip: flex-1 fills cell height minus label
                // grid (≥1400px): aspect-ratio derives height from column width
                "flex-1 min-h-0 min-[1400px]:flex-none min-[1400px]:aspect-[3/4]",
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
                size="thumb"
                className="w-full h-full"
              />
              {isIndexing && (
                <span
                  className="absolute top-[5px] left-[5px] rounded-[3px]
                             bg-accent/85 text-bg text-[9px] font-semibold
                             px-1.5 py-0.5 leading-snug backdrop-blur-sm"
                >
                  ⊙ Indexing
                </span>
              )}
            </div>
            <span className="text-[9px] text-fg-muted truncate w-full text-center">
              {e.filename?.replace(/\.dat$/i, "") ?? `#${e.id}`}
            </span>
          </div>
        );
      })}
    </div>
  );
}
