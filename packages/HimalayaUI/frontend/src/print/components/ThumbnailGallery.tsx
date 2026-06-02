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
  /** Child thumbnail size; forwarded to every Thumbnail. Default `"sm"`. */
  size?: "sm" | "lg";
  /** Horizontal alignment of the strip. `"center"` adds `justify-center`
   *  (e.g. the loupe strip). Default `"start"`. */
  align?: "start" | "center";
  /** PLACEMENT ONLY. */
  className?: string;
}

export function ThumbnailGallery({
  exposures,
  selectedId,
  onSelect,
  size = "sm",
  align = "start",
  className,
}: ThumbnailGalleryProps): JSX.Element {
  const baseClass =
    align === "center"
      ? "flex flex-nowrap overflow-x-auto gap-2 justify-center"
      : "flex flex-nowrap overflow-x-auto gap-2";

  return (
    <div
      data-testid="thumbnail-gallery"
      data-size={size}
      data-align={align}
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
          size={size}
          {...(onSelect ? { onClick: () => onSelect(exposure.id) } : {})}
        />
      ))}
    </div>
  );
}
