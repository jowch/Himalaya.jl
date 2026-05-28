import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";
import { RejectXMark } from "./RejectXMark";
import { ThumbnailGallery } from "./ThumbnailGallery";

interface Props {
  /** The exposure shown in the big frame. */
  exposure: Exposure;
  /** All exposures for the sample — drives the strip. */
  exposures: Exposure[];
  onSelectExposure: (id: number) => void;
}

/**
 * LoupeFrame — the left column of the loupe: the full-size detector image
 * for the active exposure plus a thumbnail strip of every exposure.
 *
 * Iterates `exposures` by id and renders every one regardless of `kind`
 * (`file` / `averaged` / `background_subtracted`) — the loupe must not
 * assume one file per exposure (master plan §11; the deferred
 * derived-exposure feature must remain possible).
 */
export function LoupeFrame({
  exposure,
  exposures,
  onSelectExposure,
}: Props): JSX.Element {
  const isRejected = exposure.status === "rejected";
  const caption = exposure.filename ?? `Exposure #${exposure.id}`;

  return (
    <div data-testid="loupe-frame" className="flex flex-col gap-3">
      <div
        className={[
          "relative mx-auto aspect-square w-full max-w-[500px]",
          // T-1/T-8: the detector frame is the dark window set into the paper
          // (mockup `.big-frame` = frame-edge), not light paper.
          "overflow-hidden rounded border border-frame-edge bg-frame-edge",
        ].join(" ")}
      >
        <DetectorImage
          exposureId={exposure.id}
          imagePath={exposure.image_path}
          imageVersion={exposure.image_version}
          size="full"
          className={isRejected ? "h-full w-full opacity-40" : "h-full w-full"}
        />
        {/* M-10: the hand-skewed grease-pencil ✕ over a dropped big frame —
            the same mark as the contact sheet (mockup `.big-x`), scaled up by
            the SVG's full-bleed viewBox. */}
        {isRejected && <RejectXMark testId="loupe-reject-xmark" />}
        {isRejected && (
          <span
            data-testid="loupe-dropped-badge"
            className="absolute left-3 top-3 rounded bg-print-accent px-2 py-0.5
                       text-[10px] font-bold uppercase tracking-wide text-paper"
          >
            Dropped
          </span>
        )}
        {/* T-8: light caption over the dark frame (was dark `text-ink-faint`). */}
        <span className="absolute bottom-2 left-3 font-mono text-[11px] text-frame-tag">
          {caption}
        </span>
      </div>
      <div className="h-[88px]">
        <ThumbnailGallery
          exposures={exposures}
          selectedId={exposure.id}
          onSelect={onSelectExposure}
          className="justify-center"
        />
      </div>
    </div>
  );
}
