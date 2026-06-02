import { SpecCell } from "./SpecCell";
import { ThumbnailGallery } from "./ThumbnailGallery";
import type { GalleryExposure } from "./ThumbnailGallery";
import { KeptCell } from "./KeptCell";
import { StatusCell } from "./StatusCell";
import { TagList } from "../ui/TagList";
import type { Tag } from "../ui/tag";

export interface SampleTableRowProps {
  name: string;
  sampleId: string;
  screened?: boolean;
  exposures: GalleryExposure[];
  selectedExposureId?: number;
  onSelectExposure?: (id: number) => void;
  kept: number;
  total: number;
  dropped?: number;
  tags: Tag[];
  onAddTag?: (t: Tag) => void;
  onRemoveTag?: (t: Tag) => void;
  phase?: string | null;
  /** PLACEMENT-ONLY. Appended last to the root. */
  className?: string;
}

/** The 5-column grid track shared by every row AND the table header — keeping
 *  them identical is what proves column alignment. Sample / Exposures / Kept /
 *  Tags / Status. */
export const SAMPLE_TABLE_COLS =
  "grid-cols-[244px_minmax(360px,1fr)_78px_168px_150px]";

/** Per-cell grid-child wrapper: vertical centering + the cell gutter + the row's
 *  resting height floor. All placement/layout (allowed by the design guard). */
const CELL = "flex items-center px-4 py-[13px] min-h-[92px]";

/** One full row of the contact-sheet samples table. Composes the cell composites
 *  (SpecCell / KeptCell / StatusCell), the ThumbnailGallery, and TagList across a
 *  fixed 5-column grid. Unscreened rows carry a faint resting tint; screened rows
 *  sit transparent over the table plate. */
export function SampleTableRow({
  name,
  sampleId,
  screened = false,
  exposures,
  selectedExposureId,
  onSelectExposure,
  kept,
  total,
  dropped,
  tags,
  onAddTag,
  onRemoveTag,
  phase,
  className,
}: SampleTableRowProps): JSX.Element {
  const restTint = screened ? "" : " bg-paper-sunk";

  return (
    <div
      data-testid="sample-table-row"
      data-screened={screened ? "true" : "false"}
      className={`border-b border-hair hover:bg-paper-sunk${restTint}${className ? ` ${className}` : ""}`}
    >
      <div className={`grid ${SAMPLE_TABLE_COLS}`}>
        <div className={CELL}>
          <SpecCell name={name} sampleId={sampleId} screened={screened} />
        </div>
        <div className={CELL}>
          <ThumbnailGallery
            variant="sheet"
            exposures={exposures}
            {...(selectedExposureId != null ? { selectedId: selectedExposureId } : {})}
            {...(onSelectExposure ? { onSelect: onSelectExposure } : {})}
          />
        </div>
        <div className={CELL}>
          <KeptCell
            kept={kept}
            total={total}
            {...(dropped != null ? { dropped } : {})}
          />
        </div>
        <div className={CELL}>
          <TagList
            tags={tags}
            editable
            {...(onAddTag ? { onAdd: onAddTag } : {})}
            {...(onRemoveTag ? { onRemove: onRemoveTag } : {})}
          />
        </div>
        <div className={CELL}>
          <StatusCell {...(phase !== undefined ? { phase } : {})} />
        </div>
      </div>
    </div>
  );
}
