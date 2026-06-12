import { useLayoutEffect, useRef } from "react";
import { SpecCell } from "./SpecCell";
import { ThumbnailGallery } from "./ThumbnailGallery";
import type { GalleryExposure } from "./ThumbnailGallery";
import { KeptCell } from "./KeptCell";
import { StatusCell } from "./StatusCell";
import { TagList } from "../ui/TagList";
import type { Tag } from "../ui/tag";
import { Checkbox } from "../ui";
import { useRovingGridContext } from "../../lib/grid/useRovingGrid";

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
  /** This row's grid row index (header = 0, first data row = 1, …). The page
   *  passes `i+1`. Only consulted when rendered inside a roving SheetTable grid;
   *  absent / outside the grid context → the row renders unchanged (`role="cell"`,
   *  natural tab order). */
  rowIndex?: number;
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

/** Fixed Samples-grid column indices for the six cells (checkbox is col 0 only
 *  when a checkbox column is present; otherwise Sample is the left-most cell
 *  at col 1 and col 0 is unused for this row). */
const COL_CHECKBOX = 0;
const COL_SAMPLE = 1;
const COL_EXPOSURES = 2;
const COL_KEPT = 3;
const COL_TAGS = 4;
const COL_STATUS = 5;

/**
 * Hook a multi-widget gridcell (Exposures gallery / Tags list) into the roving
 * model: in NAVIGATION mode every focusable descendant is made inert
 * (`tabIndex = -1`) so the cell is a single roving stop; in INTERACTION mode the
 * overrides are removed so the cell's own controls become naturally tabbable
 * again. The gridcell div itself carries the roving tabindex in nav mode.
 *
 * Imperative (a layout effect querying descendants) rather than threading a prop
 * through Thumbnail/TagPill/Chip — those controls do NOT self-set tabIndex
 * (verified), so React never fights the imperative value. `interactionActive` is
 * derived from the roving context, so SSE re-renders re-run the effect
 * idempotently.
 */
function useMultiWidgetInertness(
  ref: React.RefObject<HTMLElement>,
  interactionActive: boolean,
  roving: boolean,
  // Re-run when the cell's content can change (exposure/tag counts).
  ...deps: unknown[]
): void {
  useLayoutEffect(() => {
    if (!roving) return;
    const el = ref.current;
    if (el === null) return;
    const focusables = el.querySelectorAll<HTMLElement>(
      "button, input, a, select, textarea, [tabindex]",
    );
    focusables.forEach((f, i) => {
      if (interactionActive) {
        // Restore natural tabbing: first widget tabbable, the rest reachable
        // via Tab (tabIndex 0 keeps DOM order; -1 would drop them). Removing the
        // attribute returns each to its element-default tab order.
        if (i === 0) f.tabIndex = 0;
        else f.removeAttribute("tabindex");
      } else {
        f.tabIndex = -1;
      }
    });
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [ref, interactionActive, roving, ...deps]);
}

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
  rowIndex,
  className,
}: SampleTableRowProps): JSX.Element {
  const restTint = screened ? "" : " bg-paper-sunk";
  const hasCheckbox = onCheck !== undefined;
  // Sticky identity cells (SheetTable owns the scroller; rows own the frozen
  // cells). The opaque background must mirror the row's own surface — bg-plate
  // for screened rows (transparent over the Card plate), bg-paper-sunk for the
  // unscreened resting tint — and follow the row hover via group-hover, so the
  // frozen column is indistinguishable from the row at wide viewports.
  const stickyBg = screened ? "bg-plate" : "bg-paper-sunk";
  const sticky = `sticky z-10 ${stickyBg} group-hover/row:bg-paper-sunk`;
  const sampleLeft = hasCheckbox ? CHECKBOX_TRACK_PX : 0;

  // ── Roving grid wiring (SA-ROVING). Context is INERT outside a roving
  // SheetTable, so `tabIndexFor` returns undefined and this row renders exactly
  // as before (role="cell", natural tab order). ───────────────────────────────
  const grid = useRovingGridContext();
  const r = rowIndex ?? -1;
  // "On" only when both the page supplies a rowIndex AND the grid context is
  // live (tabIndexFor returns a number, never undefined).
  const roving = rowIndex !== undefined && grid.tabIndexFor(r, COL_SAMPLE) !== undefined;
  const cellRole = roving ? "gridcell" : "cell";

  // Roving accessors (undefined when not roving → spread to nothing).
  const tab = (col: number): number | undefined =>
    roving ? grid.tabIndexFor(r, col) : undefined;
  const cellRef =
    (col: number) =>
    (el: HTMLElement | null): void => {
      if (roving) grid.registerCellEl(r, col, el);
    };
  // Pointer parity: mousedown makes this the active roving cell WITHOUT forcing
  // focus (`{ focus: false }`). The native click already focuses the clicked
  // element; programmatically re-focusing here would yank focus to the cell div.
  // NOT onFocus — an onFocus latch turned every stray focus during a Tab
  // transition into a permanent active-cell migration, trapping Tab inside the
  // grid (BUG D). The active cell is discovered by the container arrow handler,
  // never by a focus event.
  const activatePointer = (col: number): (() => void) | undefined =>
    roving ? () => grid.requestActivate(r, col, { focus: false }) : undefined;
  // A multi-widget gridcell enters interaction mode on Enter/F2 (the container
  // also handles this off activeCoord; this keydown handles the case where focus
  // is on the gridcell div directly).
  const cellKeyDown =
    (col: number) =>
    (e: React.KeyboardEvent): void => {
      if (!roving) return;
      if ((e.key === "Enter" || e.key === "F2") && (col === COL_EXPOSURES || col === COL_TAGS)) {
        e.preventDefault();
        grid.enterInteraction(r, col);
      }
    };

  const exposuresInteraction =
    roving && grid.interactionCoord?.row === r && grid.interactionCoord?.col === COL_EXPOSURES;
  const tagsInteraction =
    roving && grid.interactionCoord?.row === r && grid.interactionCoord?.col === COL_TAGS;

  // Multi-widget cell refs for the imperative inertness effect.
  const exposuresCellRef = useRef<HTMLElement | null>(null);
  const tagsCellRef = useRef<HTMLElement | null>(null);
  // Depend on a CONTENT signature, not the count: a same-length swap (different
  // ids, same length) mounts fresh focusable <button>s and the effect must
  // re-run to re-inert them, else a transient extra tab stop appears in nav
  // mode. Tag has no `label` field — its identity is key(+value).
  const exposuresSig = exposures.map((e) => e.id).join(",");
  const tagsSig = tags.map((t) => `${t.key}=${t.value ?? ""}`).join(",");
  useMultiWidgetInertness(exposuresCellRef, exposuresInteraction, roving, exposuresSig);
  useMultiWidgetInertness(tagsCellRef, tagsInteraction, roving, tagsSig);

  // Compose two refs (the registry callback + the local ref for inertness).
  const exposuresRef = (el: HTMLElement | null): void => {
    exposuresCellRef.current = el;
    cellRef(COL_EXPOSURES)(el);
  };
  const tagsRef = (el: HTMLElement | null): void => {
    tagsCellRef.current = el;
    cellRef(COL_TAGS)(el);
  };

  // The Status door / Kept / multi-widget cells take the roving tabIndex on the
  // gridcell div (Kept/Exposures/Tags) or the inner widget (Status door, Sample,
  // checkbox). Spread helper keeps the non-roving render byte-identical.
  const cellRoving = (col: number) =>
    roving
      ? {
          tabIndex: tab(col),
          ref: cellRef(col) as React.Ref<HTMLDivElement>,
          onMouseDown: activatePointer(col),
        }
      : {};

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
            role={cellRole}
            data-sticky="true"
            className={`${sticky} left-0 flex items-center justify-center px-1 h-[92px]`}
            {...(roving ? { onMouseDown: activatePointer(COL_CHECKBOX) } : {})}
          >
            <Checkbox
              checked={checked ?? false}
              {...(indeterminate ? { indeterminate: true } : {})}
              onChange={onCheck}
              aria-label="Select sample"
              {...(roving ? { tabIndex: tab(COL_CHECKBOX) ?? -1, ref: cellRef(COL_CHECKBOX) } : {})}
            />
          </div>
        )}
        <div
          role={cellRole}
          data-sticky="true"
          className={`${CELL} ${sticky} border-r border-hair-strong`}
          style={{ left: sampleLeft }}
          {...(roving ? { onMouseDown: activatePointer(COL_SAMPLE) } : {})}
        >
          <SpecCell
            name={name}
            sampleId={sampleId}
            screened={screened}
            {...(onOpenLoupe ? { onOpenLoupe } : {})}
            {...(roving ? { tabIndex: tab(COL_SAMPLE) ?? -1, ref: cellRef(COL_SAMPLE) } : {})}
          />
        </div>
        <div
          role={cellRole}
          className={CELL}
          {...(roving
            ? {
                tabIndex: tab(COL_EXPOSURES),
                ref: exposuresRef as React.Ref<HTMLDivElement>,
                onMouseDown: activatePointer(COL_EXPOSURES),
                onKeyDown: cellKeyDown(COL_EXPOSURES),
              }
            : {})}
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
        <div role={cellRole} className={CELL} {...cellRoving(COL_KEPT)}>
          <KeptCell
            kept={kept}
            total={total}
            {...(dropped != null ? { dropped } : {})}
          />
        </div>
        <div
          role={cellRole}
          className={CELL}
          {...(roving
            ? {
                tabIndex: tab(COL_TAGS),
                ref: tagsRef as React.Ref<HTMLDivElement>,
                onMouseDown: activatePointer(COL_TAGS),
                onKeyDown: cellKeyDown(COL_TAGS),
              }
            : {})}
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
          role={cellRole}
          className={CELL}
          {...(roving && !onOpenFocus ? { onMouseDown: activatePointer(COL_STATUS) } : {})}
          {...(roving && !onOpenFocus ? { tabIndex: tab(COL_STATUS), ref: cellRef(COL_STATUS) as React.Ref<HTMLDivElement> } : {})}
        >
          {onOpenFocus ? (
            <button
              type="button"
              onClick={onOpenFocus}
              aria-label={phase ? `Open indexing for ${name} (${phase})` : `Index ${name}`}
              className="flex items-center rounded-md px-2 -mx-2 transition-colors hover:bg-plate/60"
              {...(roving
                ? { tabIndex: tab(COL_STATUS) ?? -1, ref: cellRef(COL_STATUS) as React.Ref<HTMLButtonElement>, onMouseDown: activatePointer(COL_STATUS) }
                : {})}
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
