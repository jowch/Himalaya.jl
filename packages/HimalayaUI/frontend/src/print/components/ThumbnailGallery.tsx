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
  /** Multi-select highlight set (a contact-sheet model where several frames are
   *  flagged at once). A thumb renders `selected` if its id is in this set OR
   *  equals `selectedId`. Optional — the Corpus page drives per-exposure verdicts
   *  off the frame cursor, not this set. */
  selectedIds?: ReadonlySet<number>;
  /** The roving keyboard cursor's active frame (distinct from selection). Renders
   *  a double border on that thumb so the cursor is legible even on a selected
   *  one. */
  cursoredId?: number;
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
  cursoredId,
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
      // Scrollable containers (overflow-x-auto) are auto-added to the sequential
      // Tab order by Chrome/blink (so a sighted keyboard user can scroll them) —
      // EVEN without a tabindex. An explicit tabindex=-1 opts this scroller OUT
      // of that auto Tab-stop. It does NOT affect children (the thumbnail buttons
      // stay tabbable in the Samples grid's interaction mode), and keyboard
      // scroll-into-view is preserved. Without this, Tab from the Samples grid's
      // single active cell would stop on the next overflowing gallery scroller
      // and never exit the grid (SA-ROVING BUG D).
      tabIndex={-1}
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
          {...(cursoredId === exposure.id ? { cursored: true } : {})}
          size={size}
          {...(onSelect ? { onClick: () => onSelect(exposure.id) } : {})}
          {...(onActivate ? { onDoubleClick: () => onActivate(exposure.id) } : {})}
        />
      ))}
    </div>
  );
}
