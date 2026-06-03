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
  /** Single highlighted exposure (loupe "current frame" model). */
  selectedExposureId?: number;
  /** Multi-select highlight set (contact-sheet *cull* model — several frames in
   *  this row flagged for dropping at once). Forwarded to the gallery's
   *  `selectedIds`; OR'd with `selectedExposureId`. */
  selectedExposureIds?: ReadonlySet<number>;
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
 *  Tags / Status. Applied via inline `style={{ gridTemplateColumns }}` rather
 *  than a `grid-cols-[…]` arbitrary so it is immune to Tailwind dev-JIT misses
 *  on a brand-new arbitrary class mid-session.
 *
 *  Per-column MIN-WIDTHS (not fixed px that clip): atomic columns are sized to
 *  fit their longest realistic value so they NEVER truncate — Kept holds
 *  "10 / 12" + "2 dropped" (96px), Status holds "Not indexed" (150px). The
 *  identity / exposures / tags columns carry a min + grow (fr). The row's
 *  intrinsic min-width is the sum of mins (~1018px); a future SheetTable wraps
 *  the rows in a horizontal-scroll container with a sticky Sample column, so a
 *  narrow viewport SCROLLS rather than clipping (Carbon/Polaris). */
export const SAMPLE_TABLE_COLS =
  "minmax(244px,1.4fr) minmax(360px,2fr) 96px minmax(168px,1fr) 150px";

/** Per-cell grid-child wrapper: vertical centering + the cell gutter + a DEFINED
 *  row height (not a min — keeps every row uniform). `min-w-0` lets the exposures
 *  gallery scroll instead of widening the track. NO `overflow-hidden` here on
 *  purpose: a cell-level clip silently truncates data AND clips escaping UI like
 *  the tags "+N" tooltip. Instead EACH cell owns its overflow — the name clamps
 *  (line-clamp-2), the id truncates, KeptCell self-clips its short data, the
 *  gallery scrolls (overflow-x-auto), and TagList caps with "+N". All cells thus
 *  self-bound to the row height. Placement/layout only (design-guard clean). */
const CELL = "flex items-center px-3 h-[92px] min-w-0";

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
  selectedExposureIds,
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
      <div className="grid" style={{ gridTemplateColumns: SAMPLE_TABLE_COLS }}>
        <div className={CELL}>
          <SpecCell name={name} sampleId={sampleId} screened={screened} />
        </div>
        <div className={CELL}>
          <ThumbnailGallery
            size="sm"
            className="w-full"
            exposures={exposures}
            {...(selectedExposureId != null ? { selectedId: selectedExposureId } : {})}
            {...(selectedExposureIds != null ? { selectedIds: selectedExposureIds } : {})}
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
            maxVisible={2}
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
