import { forwardRef } from "react";
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
  /** Multi-select highlight set, forwarded to the gallery's `selectedIds`.
   *  Optional — the Corpus page drives per-exposure verdicts off the frame
   *  cursor, not this set; kept for the contact-sheet assembly story. */
  selectedExposureIds?: ReadonlySet<number>;
  /** The roving cursor's active frame within this row → a double-border cue on
   *  that thumb, distinct from the solid cull-selection border. */
  cursoredExposureId?: number;
  onSelectExposure?: (id: number) => void;
  kept: number;
  total: number;
  dropped?: number;
  /** The sample's exposures are LOADED and empty (SA-ZEROEXP) — there is nothing
   *  to index. Renders a terminal "No exposures" status. Supplied by the page only
   *  once the exposure fetch has resolved empty (never inferred from `total === 0`,
   *  which is also true mid-load). Default false. */
  noExposures?: boolean;
  /** Slot coordinate shown as a chip in the identity cell (e.g. "slot 5").
   *  When undefined, no slot chip renders. */
  slotIndex?: number;
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
  /** When set, double-clicking a thumbnail fires this with the exposure id,
   *  opening the loupe at that frame. Forwarded to ThumbnailGallery.onActivate. */
  onActivateExposure?: (id: number) => void;
  /** When provided and `dropped > 0`, a Restore button appears in the Kept
   *  cell — the pointer target for the Backspace=restore key. Mirrors the
   *  CullBar restore verb (null-status batchSet on the sample's exposures). */
  onRestore?: () => void;
  /** When provided, a checkbox column is rendered as the left-most cell. */
  checked?: boolean;
  indeterminate?: boolean;
  onCheck?: () => void;
  /** The roving keyboard cursor (↑/↓) is on this row → an inset accent outline
   *  so arrow-key navigation is visible (item 4). The cursor's active frame is
   *  separately highlighted by `selectedExposureId`. */
  cursored?: boolean;
  /** PLACEMENT-ONLY. Appended last to the root. */
  className?: string;
  // Roving-cursor DOM props — spread from cursor.rowProps(id) after destructuring ref.
  // SheetTable/SampleFold pass none of these → unchanged behavior.
  tabIndex?: 0 | -1;
  /** Root-level click for cursor parking (from rowProps). */
  onClick?: (e: import("react").MouseEvent<HTMLElement>) => void;
  onDoubleClick?: (e: import("react").MouseEvent<HTMLElement>) => void;
  "aria-current"?: "true" | undefined;
  "data-cursored"?: "true" | "false";
  role?: "row";
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
export const SampleTableRow = forwardRef<HTMLDivElement, SampleTableRowProps>(function SampleTableRow({
  name,
  sampleId,
  screened = false,
  exposures,
  selectedExposureId,
  selectedExposureIds,
  cursoredExposureId,
  onSelectExposure,
  kept,
  total,
  dropped,
  noExposures = false,
  slotIndex,
  tags,
  onAddTag,
  onRemoveTag,
  phase,
  formFactor = false,
  onOpenLoupe,
  onActivateExposure,
  onRestore,
  checked,
  indeterminate,
  onCheck,
  cursored = false,
  tabIndex,
  onClick,
  onDoubleClick,
  "aria-current": ariaCurrent,
  "data-cursored": dataCursored,
  role: roleProp,
  className,
}: SampleTableRowProps, ref): JSX.Element {
  const restTint = screened ? "" : " bg-paper-sunk";
  const hasCheckbox = onCheck !== undefined;
  // Cursor cue (item 4): the roving ↑/↓ row reads as one opaque yellowish band
  // (bg-row-cursor, the pages2 mockup's accent-into-paper-sunk mix). An `outline`
  // was clipped — the frozen sticky cells (z-10, opaque) painted over it on the
  // left and the horizontal scroller cut it on the right — so the whole row + its
  // sticky cells take the OPAQUE band instead (opaque so the frozen column still
  // hides content beneath). Distinct in hue from the neutral grey rest/hover, and
  // overrides hover so the cursor stays put under the mouse.
  const rowBg = cursored ? " bg-row-cursor" : ` hover:bg-paper-sunk${restTint}`;
  // The root also carries `scroll-mb-14` (56px): when the page scrolls the
  // cursored row into view (scrollIntoView), that bottom scroll-margin keeps a
  // downward-navigated row ABOVE the fixed Dock (~47px) instead of behind it.

  // Sticky identity cells (SheetTable owns the scroller; rows own the frozen
  // cells). The opaque background must mirror the row's own surface — bg-plate
  // for screened rows (transparent over the Card plate), bg-paper-sunk for the
  // unscreened resting tint — and follow the row hover via group-hover, so the
  // frozen column is indistinguishable from the row at wide viewports.
  const stickyBg = screened ? "bg-plate" : "bg-paper-sunk";
  // The frozen cells follow the cursor band so the left column doesn't stay grey
  // while the body turns terracotta (the "broken border" the outline produced).
  const sticky = cursored
    ? "sticky z-10 bg-row-cursor"
    : `sticky z-10 ${stickyBg} group-hover/row:bg-paper-sunk`;
  const sampleLeft = hasCheckbox ? CHECKBOX_TRACK_PX : 0;

  return (
    <div
      ref={ref}
      data-testid="sample-table-row"
      role={roleProp ?? "row"}
      data-screened={screened ? "true" : "false"}
      data-cursored={dataCursored ?? (cursored ? "true" : "false")}
      {...(tabIndex !== undefined ? { tabIndex } : {})}
      {...(onClick ? { onClick } : {})}
      {...(onDoubleClick ? { onDoubleClick } : {})}
      {...(ariaCurrent !== undefined ? { "aria-current": ariaCurrent } : {})}
      // The roving cursor focuses this div programmatically; the bg-row-cursor
      // band (above) IS the focus indicator, so suppress the UA outline that would
      // otherwise paint a clipped blue box over the sticky cells (the very clipping
      // the band was chosen to avoid).
      className={`group/row scroll-mb-14 border-b border-hair focus:outline-none${rowBg}${className ? ` ${className}` : ""}`}
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
            {...(slotIndex !== undefined ? { slotIndex } : {})}
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
            {...(cursoredExposureId != null ? { cursoredId: cursoredExposureId } : {})}
            {...(onSelectExposure ? { onSelect: onSelectExposure } : {})}
            {...(onActivateExposure ? { onActivate: onActivateExposure } : {})}
          />
        </div>
        <div role="cell" className={CELL}>
          <KeptCell
            kept={kept}
            total={total}
            {...(dropped != null ? { dropped } : {})}
            {...(onRestore ? { onRestore } : {})}
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
          <StatusCell
              {...(phase !== undefined ? { phase } : {})}
              {...(formFactor ? { formFactor: true } : {})}
              {...(noExposures ? { noExposures: true } : {})}
            />
        </div>
      </div>
    </div>
  );
});
