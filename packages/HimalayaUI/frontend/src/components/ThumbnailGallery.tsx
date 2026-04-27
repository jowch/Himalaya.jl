import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";

interface Props {
  exposures: Exposure[];
  selectedId: number | undefined;
  onSelect: (id: number) => void;
  /** 1 = horizontal strip, 2 = two-column grid, "auto" = auto-fill grid */
  columns: 1 | 2 | "auto";
  className?: string;
}

export function ThumbnailGallery({
  exposures,
  selectedId,
  onSelect,
  columns,
  className,
}: Props): JSX.Element {
  const isFilmstrip = columns === 1;

  const gridClass = isFilmstrip
    ? "grid grid-flow-col grid-rows-1 min-[1400px]:grid-rows-2 gap-2 overflow-x-auto h-full"
    : columns === 2
      ? "grid grid-cols-2 gap-2"
      : "grid grid-cols-[repeat(auto-fill,minmax(72px,1fr))] gap-2";

  return (
    <div className={`${gridClass} ${className ?? ""}`}>
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
              isFilmstrip ? "aspect-[3/4]" : "shrink-0",
              isRejected ? "opacity-40 grayscale" : "",
            ].join(" ")}
          >
            <div
              className={[
                "relative w-full overflow-hidden rounded-md transition-all duration-150",
                isFilmstrip ? "flex-1 min-h-0" : "aspect-[3/4]",
                isViewing
                  ? "ring-2 ring-accent ring-offset-1 ring-offset-bg-elevated"
                  : "group-hover:ring-1 group-hover:ring-border",
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
