import { Thumbnail } from "./Thumbnail";

export interface GalleryExposure {
  id: number;
  src: string | null;
  frameNo?: string | number;
  representative?: boolean;
  rejected?: boolean;
}

export interface ThumbnailGalleryProps {
  exposures: GalleryExposure[];
  /** Currently-selected exposure id (drives child Thumbnail `selected`). */
  selectedId?: number;
  onSelect?: (id: number) => void;
  variant?: "sheet" | "loupe";
  /** PLACEMENT ONLY. */
  className?: string;
}

export function ThumbnailGallery({
  exposures,
  selectedId,
  onSelect,
  variant = "sheet",
  className,
}: ThumbnailGalleryProps): JSX.Element {
  const baseClass =
    variant === "loupe"
      ? "flex justify-center gap-2 mt-3"
      : "flex flex-nowrap overflow-x-auto gap-2";

  return (
    <div
      data-testid="thumbnail-gallery"
      data-variant={variant}
      className={`${baseClass}${className ? ` ${className}` : ""}`}
    >
      {exposures.map((exposure) => (
        <Thumbnail
          key={exposure.id}
          src={exposure.src}
          {...(exposure.frameNo != null ? { frameNo: exposure.frameNo } : {})}
          {...(exposure.representative != null ? { representative: exposure.representative } : {})}
          {...(exposure.rejected != null ? { rejected: exposure.rejected } : {})}
          selected={selectedId === exposure.id}
          size={variant}
          {...(onSelect ? { onClick: () => onSelect(exposure.id) } : {})}
        />
      ))}
    </div>
  );
}
