import type { Exposure } from "../api";
import { DetectorImage } from "./DetectorImage";
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
          "overflow-hidden rounded border border-border bg-bg",
        ].join(" ")}
      >
        <DetectorImage
          exposureId={exposure.id}
          imagePath={exposure.image_path}
          imageVersion={exposure.image_version}
          size="full"
          className={isRejected ? "h-full w-full opacity-40" : "h-full w-full"}
        />
        {isRejected && (
          <span
            data-testid="loupe-dropped-badge"
            className="absolute left-3 top-3 rounded bg-accent px-2 py-0.5
                       text-[10px] font-bold uppercase tracking-wide text-bg"
          >
            Dropped
          </span>
        )}
        <span className="absolute bottom-2 left-3 font-mono text-[11px] text-ink-faint">
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
