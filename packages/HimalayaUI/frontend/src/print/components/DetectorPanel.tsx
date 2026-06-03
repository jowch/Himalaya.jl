// src/print/components/DetectorPanel.tsx
import type { ReactNode } from "react";
import { Card } from "../ui";
import { DetectorImage, DetectorRings, buildRingPlacements } from "../detector";
import type { DetectorLutVariant } from "../detector/detectorLut";
import type { DetectorCalibration, RingInput } from "../detector/detectorGeometry";
import { PanelHeader } from "./PanelHeader";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface DetectorPanelProps {
  /** Current exposure image URL; null → DetectorImage placeholder. */
  src: string | null;
  /** Section label. Default "Detector image". */
  label?: ReactNode;
  /** Header right slot — the exposure switcher (ThumbnailGallery). */
  tools?: ReactNode;
  /** Detector colormap; "neutral" (default) lets the ring overlay own colour. */
  lutVariant?: DetectorLutVariant;
  /** Phase-coloured Debye–Scherrer rings to overlay — the assigned reflections,
   *  each `{ q, color, ghost? }` (or a bare q). Empty/omitted → image only. */
  rings?: RingInput[];
  /** Beam-center + q→radius calibration. `null`/omitted → presentational
   *  centered fallback (radii spread by q-range) until the backend supplies it. */
  calibration?: DetectorCalibration | null;
  /** Normalized (0–1) beam center, overriding the placement's. Use when the beam
   *  center is MEASURED but full optical calibration (for true radii) isn't yet —
   *  the real intermediate state. Omitted → the placement's center (calibration
   *  or centered fallback). */
  beamCenter?: { x: number; y: number };
  /** Displayed image aspect ratio (width / height). Sets the frame's box ratio so
   *  the contained image fills it with no letterbox (rings register), and corrects
   *  the ring y-radius so they stay round. Default 1 (square). */
  imageAspect?: number;
  /** q hovered elsewhere (trace peak / comb tooth) → the matching ring lights —
   *  the cross-panel "triple lights up" link. */
  hoveredQ?: number;
  /** Fired on ring enter (q) / leave (undefined) → drives the same link outward. */
  onHoverQ?: (q?: number) => void;
  /** Caption under the frame. */
  hint?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function DetectorPanel({
  src,
  label = "Detector image",
  tools,
  lutVariant,
  rings,
  calibration,
  beamCenter,
  imageAspect,
  hoveredQ,
  onHoverQ,
  hint,
  className,
}: DetectorPanelProps): JSX.Element {
  // Rings are geometrically real (beam-centered). null calibration → centered
  // presentational fallback radii; a measured `beamCenter` overrides the center
  // (real-center + fallback-radii hybrid until full optical calibration lands).
  const placed =
    rings && rings.length > 0
      ? buildRingPlacements(rings, calibration ?? null)
      : null;
  const center = beamCenter ?? placed?.beamCenter;
  const aspect = imageAspect ?? 1;
  return (
    <Card
      as="section"
      data-testid="detector-panel"
      className={cx("flex flex-col p-4", className)}
    >
      <PanelHeader label={label}>{tools}</PanelHeader>
      <div
        data-testid="detector-frame"
        // aspectRatio is geometry (not appearance) → inline style is guard-clean
        // and lets the frame match the real (often portrait) detector image so
        // the contained image fills edge-to-edge and the rings register. Bounded
        // to the mockup's det-size (320px) so the frame — esp. the empty
        // placeholder — never blows up to fill an unbounded container; floored so
        // it stays legible in a narrow column.
        style={{ aspectRatio: aspect, maxWidth: 320, minWidth: 160 }}
        className="relative bg-frame-edge border border-frame-edge rounded overflow-hidden mx-auto w-full"
      >
        <DetectorImage
          src={src}
          size="full"
          className="h-full w-full"
          {...(lutVariant !== undefined ? { lutVariant } : {})}
        />
        {placed && center && (
          <DetectorRings
            beamCenter={center}
            rings={placed.rings}
            imageAspect={aspect}
            {...(hoveredQ !== undefined ? { hoveredQ } : {})}
            {...(onHoverQ !== undefined ? { onHoverQ } : {})}
          />
        )}
      </div>
      {hint != null && (
        <div className="text-caption text-ink-faint mt-2.5">{hint}</div>
      )}
    </Card>
  );
}
