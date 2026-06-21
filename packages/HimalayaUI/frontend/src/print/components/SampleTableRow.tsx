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
  /** The sample's exposures are LOADED and empty (SA-ZEROEXP) — there is nothing
   *  to index. Renders a terminal "No exposures" status and suppresses the Focus
   *  door even when `onOpenFocus` is wired. Supplied by the page only once the
   *  exposure fetch has resolved empty (never inferred from `total === 0`, which
   *  is also true mid-load). Default false. */
  noExposures?: boolean;
  tags: Tag[];
  onAddTag?: (t: Tag) => void;
  onRemoveTag?: (t: Tag) => void;
  phase?: string | null;
  /** The representative exposure is declared form_factor → the status cell reads
   *  "Form factor" rather than the phase chip or the "Not indexed" invite. */
  formFactor?: boolean;
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
 *  intrinsic min-width is the sum of mins (~1018px); SheetTable wraps the
 *  header + rows in ONE shared horizontal-scroll container, so a narrow
 *  viewport SCROLLS rather than clipping (Carbon/Polaris). The contract is
 *  split: SheetTable provides the scroller (`data-testid="sheet-scroll"`) and
 *  the sticky column HEADERS; every row provides its own sticky identity
 *  cells — checkbox at left 0, Sample at left CHECKBOX_TRACK_PX (or 0 when no
 *  checkbox column) — with an opaque background + z-10 so scrolled content
 *  slides beneath, and a right hairline so the frozen edge reads. */
export const SAMPLE_TABLE_COLS_BASE =
  "minmax(244px,1.4fr) minmax(360px,2fr) 96px minmax(168px,1fr) 150px";
/** Width of the checkbox grid track, in px. SINGLE SOURCE for both the grid
 *  template's leading track AND the sticky `left` offset of the Sample column
 *  (here and in SheetTable's header) — they cannot diverge. */
export const CHECKBOX_TRACK_PX = 36;
/** Returns the grid-template-columns string for a row/header.
 *  `withCheckbox=true` prepends the checkbox track as the left-most column.
 *  BOTH the SheetTable header AND every SampleTableRow body call this with the
 *  same boolean — the shared computed template is the alignment invariant. */
export function sampleTableCols(withCheckbox: boolean): string {
  return withCheckbox
    ? `${CHECKBOX_TRACK_PX}px ${SAMPLE_TABLE_COLS_BASE}`
    : SAMPLE_TABLE_COLS_BASE;
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
  noExposures = false,
  tags,
  onAddTag,
  onRemoveTag,
  phase,
  formFactor = false,
  onOpenLoupe,
  onOpenFocus,
  onActivateExposure,
  checked,
  indeterminate,
  onCheck,
  className,
}: SampleTableRowProps): JSX.Element {
  const restTint = screened ? "" : " bg-paper-sunk";
  const hasCheckbox = onCheck !== undefined;
  // SA-ZEROEXP: a sample CONFIRMED to have zero exposures has nothing to index,
  // so the status cell is NOT a Focus door — it reads a terminal "No exposures"
  // status instead of a live door into an empty workspace, and the door is
  // suppressed even when the page wires onOpenFocus. `noExposures` is supplied by
  // the page (exposures loaded AND empty); it must NOT be inferred from
  // `total === 0`, which is also true while a sample's exposures are still
  // loading — that would flash the empty status across the whole sheet mid-load.
  const isDoor = onOpenFocus != null && !noExposures;
  // Sticky identity cells (SheetTable owns the scroller; rows own the frozen
  // cells). The opaque background must mirror the row's own surface — bg-plate
  // for screened rows (transparent over the Card plate), bg-paper-sunk for the
  // unscreened resting tint — and follow the row hover via group-hover, so the
  // frozen column is indistinguishable from the row at wide viewports.
  const stickyBg = screened ? "bg-plate" : "bg-paper-sunk";
  const sticky = `sticky z-10 ${stickyBg} group-hover/row:bg-paper-sunk`;
  const sampleLeft = hasCheckbox ? CHECKBOX_TRACK_PX : 0;

  return (
    <div
      data-testid="sample-table-row"
      role="row"
      data-screened={screened ? "true" : "false"}
      className={`group/row border-b border-hair hover:bg-paper-sunk${restTint}${className ? ` ${className}` : ""}`}
    >
      <div
        role="presentation"
        className="grid"
        style={{ gridTemplateColumns: sampleTableCols(hasCheckbox) }}
      >
        {hasCheckbox && (
          <div
            role="cell"
            data-sticky="true"
            className={`${sticky} left-0 flex items-center justify-center px-1 h-[92px]`}
          >
            <Checkbox
              checked={checked ?? false}
              {...(indeterminate ? { indeterminate: true } : {})}
              onChange={onCheck}
              aria-label="Select sample"
            />
          </div>
        )}
        <div
          role="cell"
          data-sticky="true"
          className={`${CELL} ${sticky} border-r border-hair-strong`}
          style={{ left: sampleLeft }}
        >
          <SpecCell
            name={name}
            sampleId={sampleId}
            screened={screened}
            {...(onOpenLoupe ? { onOpenLoupe } : {})}
          />
        </div>
        <div
          role="cell"
          className={CELL}
        >
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
        <div role="cell" className={CELL}>
          <KeptCell
            kept={kept}
            total={total}
            {...(dropped != null ? { dropped } : {})}
          />
        </div>
        <div
          role="cell"
          className={CELL}
        >
          <TagList
            tags={tags}
            maxVisible={2}
            editable
            {...(onAddTag ? { onAdd: onAddTag } : {})}
            {...(onRemoveTag ? { onRemove: onRemoveTag } : {})}
          />
        </div>
        <div
          role="cell"
          className={CELL}
        >
          {isDoor ? (
            <button
              type="button"
              onClick={onOpenFocus}
              aria-label={phase ? `Open indexing for ${name} (${phase})` : `Index ${name}`}
              className="flex items-center rounded-md px-2 -mx-2 transition-colors hover:bg-plate/60"
            >
              <StatusCell door {...(phase !== undefined ? { phase } : {})} {...(formFactor ? { formFactor: true } : {})} />
            </button>
          ) : (
            <StatusCell
              {...(phase !== undefined ? { phase } : {})}
              {...(formFactor ? { formFactor: true } : {})}
              {...(noExposures ? { noExposures: true } : {})}
            />
          )}
        </div>
      </div>
    </div>
  );
}
