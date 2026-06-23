import { DetectorImage } from "../detector/DetectorImage";
import { RejectOverlay } from "../ui/RejectOverlay";
import { Dot } from "../ui/Dot";

export interface ThumbnailProps {
  /** Detector thumbnail source URL; null → DetectorImage placeholder. */
  src: string | null;
  /** Frame number label, e.g. "65". */
  frameNo?: string | number;
  representative?: boolean;
  rejected?: boolean;
  /** Screened-in (exposure status "accepted"). Tri-state semantics (SA-SCREENED):
   *  kept and rejected are EXPLICIT verdicts; an unscreened (null-status) frame
   *  is neither — so omit both flags rather than smearing it into either side.
   *  If both are passed, rejected wins and kept is ignored. */
  kept?: boolean;
  selected?: boolean;
  /** The roving keyboard cursor's active frame. Distinct from `selected` (the
   *  cull pick): renders a thinner DOUBLE border (accent line + inner light line)
   *  so the cursor stays legible even when it sits on an already-selected thumb,
   *  whose solid accent border would otherwise look identical. */
  cursored?: boolean;
  /** `"xs"` → 30px (dense Focus exposure strip); `"sm"` → 62px (contact-sheet strip); `"lg"` → 70px (loupe strip). */
  size?: "xs" | "sm" | "lg";
  onClick?: () => void;
  onDoubleClick?: () => void;
  title?: string;
  /** PLACEMENT ONLY. */
  className?: string;
}

const SIZE_PX: Record<"xs" | "sm" | "lg", number> = {
  xs: 30,
  sm: 62,
  lg: 70,
};

function buildDataState(props: {
  selected?: boolean;
  cursored?: boolean;
  rejected?: boolean;
  representative?: boolean;
  kept?: boolean;
}): string {
  const tokens: string[] = [];
  if (props.selected) tokens.push("selected");
  if (props.cursored) tokens.push("cursored");
  if (props.rejected) tokens.push("rejected");
  if (props.representative) tokens.push("representative");
  if (props.kept) tokens.push("kept");
  if (tokens.length === 0) tokens.push("normal");
  return tokens.join(" ");
}

export function Thumbnail({
  src,
  frameNo,
  representative = false,
  rejected = false,
  kept = false,
  selected = false,
  cursored = false,
  size = "sm",
  onClick,
  onDoubleClick,
  title,
  className,
}: ThumbnailProps): JSX.Element {
  const px = SIZE_PX[size];
  // A frame is never kept AND dropped — rejected wins so the two channels
  // can't contradict each other if a caller passes both.
  const showKept = kept && !rejected;
  const dataState = buildDataState({ selected, cursored, rejected, representative, kept: showKept });

  // Accessible name: "Frame N" (or "Exposure" when unnumbered), with state
  // suffixes for the visual markers (representative dot / reject overlay /
  // kept dot) so a screen-reader hears what a sighted user sees. Selection is
  // conveyed via aria-pressed (a toggle role), not folded into the label.
  let ariaLabel = frameNo != null ? `Frame ${frameNo}` : "Exposure";
  if (representative) ariaLabel += ", representative";
  if (rejected) ariaLabel += ", dropped";
  if (showKept) ariaLabel += ", kept";

  return (
    <button
      data-testid="thumbnail"
      data-state={dataState}
      data-size={size}
      // UI-RINGCLIP: inset the keyboard focus ring so the gallery's
      // overflow-x-auto (which clips the y-axis) can't crop its top/bottom.
      data-focus-ring="inset"
      aria-label={ariaLabel}
      {...(onClick ? { "aria-pressed": selected } : {})}
      onClick={onClick}
      {...(onDoubleClick ? { onDoubleClick } : {})}
      // SA-THUMBKEY: keyboard parity for the double-click→loupe activate. Plain
      // Enter stays the button's native click (select-toggle); Shift+Enter opens
      // the loupe AT this frame. preventDefault stops the default Enter-click so
      // the same keystroke can't ALSO toggle the selection.
      {...(onDoubleClick
        ? {
            onKeyDown: (e: React.KeyboardEvent) => {
              if (e.key === "Enter" && e.shiftKey) {
                e.preventDefault();
                onDoubleClick();
              }
            },
          }
        : {})}
      title={title}
      style={{ width: px, height: px }}
      className={`group relative inline-block flex-shrink-0 overflow-hidden rounded-sm bg-frame-edge border border-frame-edge p-0 cursor-pointer${className ? ` ${className}` : ""}`}
    >
      {/* Detector image — dimmed when rejected */}
      <div
        data-dimmed={rejected ? "true" : "false"}
        className={rejected ? "w-full h-full opacity-40" : "w-full h-full"}
      >
        <DetectorImage src={src} size="thumb" className="w-full h-full" />
      </div>

      {/* Frame-number label */}
      {frameNo != null && (
        <>
          {/* SA-RESCORE3 F10: a soft bottom scrim so the mono frame number stays
              legible over BRIGHT/busy detector content. The mockup assumed
              uniform-dark frames (no scrim); real frames can be bright enough to
              wash out the light caption. Gradient up from the frame-edge token
              (the window backing) — token-based, no raw literal. */}
          <span
            aria-hidden="true"
            data-role="thumb-fno-scrim"
            className="pointer-events-none absolute inset-x-0 bottom-0 h-1/3 bg-gradient-to-t from-frame-edge/85 to-transparent"
          />
          <span
            data-role="thumb-fno"
            className={`absolute left-1 bottom-0.5 font-mono text-frame-tag text-xs${rejected ? " opacity-50" : " opacity-80"}`}
          >
            {frameNo}
          </span>
        </>
      )}

      {/* Representative marker. SA-RESCORE3 F8: the meaning lived only in the
          button's aria-label (SR-only); a sighted mouse user saw an unexplained
          dot. A `title` surfaces it on hover too. */}
      {representative && (
        <span
          data-role="thumb-rep"
          title="Representative exposure"
          className="absolute top-0.5 right-0.5 flex"
        >
          <Dot tone="accent" size="md" bordered aria-hidden="true" />
        </span>
      )}

      {/* Kept (screened-in) marker — top-LEFT so a representative-and-kept
          frame shows both corners (bottom-left belongs to the frame number). */}
      {showKept && (
        <span
          data-role="thumb-kept"
          title="Kept (screened in)"
          className="absolute top-0.5 left-0.5 flex"
        >
          <Dot tone="success" size="md" bordered aria-hidden="true" />
        </span>
      )}

      {/* Reject overlay */}
      {rejected && <RejectOverlay />}

      {/* Selection / hover frame — an overlay painted ON TOP of the detector
          image. An inset ring (or any box-shadow) sits BEHIND child content, so
          the `object-fit: contain` canvas covered it wherever the image bled to
          the frame edge (top/bottom on a square frame holding a portrait image).
          An `inset-0` bordered span paints the full perimeter above the canvas,
          inside `overflow-hidden` — so it is neither occluded by the image nor
          clipped by the gallery's `overflow-x-auto`. */}
      <span
        aria-hidden="true"
        className={`pointer-events-none absolute inset-0 rounded-sm${
          selected
            ? " border-[3px] border-accent"
            : cursored
              ? " border-2 border-accent"
              : " border-2 border-transparent group-hover:border-hair-strong"
        }`}
      />
      {/* Cursor-on-selection = the SECOND, inner line of a double border. Only a
          SELECTED thumb needs disambiguating (its solid accent border already
          looks like the cursor's), so the inner light line is drawn only when the
          cursor sits on a selected thumb. A cursored-but-unselected thumb keeps
          the single 2px accent line above. */}
      {cursored && selected && (
        <span
          aria-hidden="true"
          className="pointer-events-none absolute inset-[3px] rounded-sm border border-paper"
        />
      )}
    </button>
  );
}
