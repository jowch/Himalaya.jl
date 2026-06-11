import { Thumbnail } from "./Thumbnail";

export interface GalleryExposure {
  id: number;
  src: string | null;
  frameNo?: string | number;
  representative?: boolean;
  rejected?: boolean;
  /** Screened-in (status "accepted"); unscreened (null status) frames omit
   *  BOTH kept and rejected — the SA-SCREENED tri-state. */
  kept?: boolean;
}

export interface ThumbnailGalleryProps {
  exposures: GalleryExposure[];
  /** Single currently-selected exposure id (the loupe "current frame" model). */
  selectedId?: number;
  /** Multi-select highlight set (the contact-sheet *cull* model, where several
   *  frames in one row can be flagged for dropping at once). A thumb renders
   *  `selected` if its id is in this set OR equals `selectedId` — so a gallery
   *  can use either model, or both. */
  selectedIds?: ReadonlySet<number>;
  onSelect?: (id: number) => void;
  /** Double-click handler — opens the loupe for the clicked exposure. */
  onActivate?: (id: number) => void;
  /** Child thumbnail size; forwarded to every Thumbnail. Default `"sm"`. */
  size?: "xs" | "sm" | "lg";
  /** Horizontal alignment of the strip. `"center"` adds `justify-center`
   *  (e.g. the loupe strip). Default `"start"`. */
  align?: "start" | "center";
  /** PLACEMENT ONLY. */
  className?: string;
}

export function ThumbnailGallery({
  exposures,
  selectedId,
  selectedIds,
  onSelect,
  onActivate,
  size = "sm",
  align = "start",
  className,
}: ThumbnailGalleryProps): JSX.Element {
  const baseClass =
    align === "center"
      ? "flex flex-nowrap overflow-x-auto gap-2 justify-center-safe min-w-0"
      : "flex flex-nowrap overflow-x-auto gap-2 min-w-0";

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
          {...(exposure.kept != null ? { kept: exposure.kept } : {})}
          selected={selectedId === exposure.id || (selectedIds?.has(exposure.id) ?? false)}
          size={size}
          {...(onSelect ? { onClick: () => onSelect(exposure.id) } : {})}
          {...(onActivate ? { onDoubleClick: () => onActivate(exposure.id) } : {})}
        />
      ))}
    </div>
  );
}
