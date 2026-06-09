import { SpecCell } from "./SpecCell";
import { ThumbnailGallery } from "./ThumbnailGallery";
import type { GalleryExposure } from "./ThumbnailGallery";
import { KeptCell } from "./KeptCell";
import { StatusCell } from "./StatusCell";
import { TagList } from "../ui/TagList";
import type { Tag } from "../ui/tag";
import { Checkbox } from "../ui";

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
  /** When set, the sample name in the SpecCell becomes a button that opens the
   *  loupe view for this sample. */
  onOpenLoupe?: () => void;
  /** When set, the status cell becomes a keyboard-accessible Focus door — a
   *  button that opens the indexing workspace for this sample. Absent → the
   *  status cell is a plain div (other consumers unaffected). */
  onOpenFocus?: () => void;
  /** When set, double-clicking a thumbnail fires this with the exposure id,
   *  opening the loupe at that frame. Forwarded to ThumbnailGallery.onActivate. */
  onActivateExposure?: (id: number) => void;
  /** When provided, a checkbox column is rendered as the left-most cell. */
  checked?: boolean;
  indeterminate?: boolean;
  onCheck?: () => void;
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
export const SAMPLE_TABLE_COLS_BASE =
  "minmax(244px,1.4fr) minmax(360px,2fr) 96px minmax(168px,1fr) 150px";
/** Returns the grid-template-columns string for a row/header.
 *  `withCheckbox=true` prepends a 36px checkbox track as the left-most column.
 *  BOTH the SheetTable header AND every SampleTableRow body call this with the
 *  same boolean — the shared computed template is the alignment invariant. */
export function sampleTableCols(withCheckbox: boolean): string {
  return withCheckbox ? `36px ${SAMPLE_TABLE_COLS_BASE}` : SAMPLE_TABLE_COLS_BASE;
}
// Backward-compat alias — existing callers (stories) import this as a string.
export const SAMPLE_TABLE_COLS = SAMPLE_TABLE_COLS_BASE;

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
  onOpenLoupe,
  onOpenFocus,
  onActivateExposure,
  checked,
  indeterminate,
  onCheck,
  className,
}: SampleTableRowProps): JSX.Element {
  const restTint = screened ? "" : " bg-paper-sunk";

  return (
    <div
      data-testid="sample-table-row"
      data-screened={screened ? "true" : "false"}
      className={`border-b border-hair hover:bg-paper-sunk${restTint}${className ? ` ${className}` : ""}`}
    >
      <div
        className="grid"
        style={{ gridTemplateColumns: sampleTableCols(onCheck !== undefined) }}
      >
        {onCheck !== undefined && (
          <div className="flex items-center justify-center px-1 h-[92px]">
            <Checkbox
              checked={checked ?? false}
              {...(indeterminate ? { indeterminate: true } : {})}
              onChange={onCheck}
              aria-label="Select sample"
            />
          </div>
        )}
        <div className={CELL}>
          <SpecCell
            name={name}
            sampleId={sampleId}
            screened={screened}
            {...(onOpenLoupe ? { onOpenLoupe } : {})}
          />
        </div>
        <div className={CELL}>
          <ThumbnailGallery
            size="sm"
            className="w-full"
            exposures={exposures}
            {...(selectedExposureId != null ? { selectedId: selectedExposureId } : {})}
            {...(selectedExposureIds != null ? { selectedIds: selectedExposureIds } : {})}
            {...(onSelectExposure ? { onSelect: onSelectExposure } : {})}
            {...(onActivateExposure ? { onActivate: onActivateExposure } : {})}
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
          {onOpenFocus ? (
            <button
              type="button"
              onClick={onOpenFocus}
              aria-label={phase ? `Open indexing for ${name} (${phase})` : `Index ${name}`}
              className="flex items-center rounded-md px-2 -mx-2 transition-colors hover:bg-plate/60"
            >
              <StatusCell {...(phase !== undefined ? { phase } : {})} />
            </button>
          ) : (
            <StatusCell {...(phase !== undefined ? { phase } : {})} />
          )}
        </div>
      </div>
    </div>
  );
}
