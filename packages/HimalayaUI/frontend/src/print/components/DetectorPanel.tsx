// src/print/components/DetectorPanel.tsx
import type { ReactNode } from "react";
import { Card } from "../ui";
import { DetectorImage } from "../detector";
import type { DetectorLutVariant } from "../detector/detectorLut";
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
  /** Detector colormap; "neutral" lets the ring overlay own colour. */
  lutVariant?: DetectorLutVariant;
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
  hint,
  className,
}: DetectorPanelProps): JSX.Element {
  return (
    <Card
      as="section"
      data-testid="detector-panel"
      className={cx("flex flex-col p-4", className)}
    >
      <PanelHeader label={label}>{tools}</PanelHeader>
      <div
        data-testid="detector-frame"
        className="bg-frame-edge border border-frame-edge rounded overflow-hidden aspect-square mx-auto w-full"
      >
        <DetectorImage
          src={src}
          size="full"
          {...(lutVariant !== undefined ? { lutVariant } : {})}
        />
      </div>
      {hint != null && (
        <div className="text-caption text-ink-faint mt-2.5">{hint}</div>
      )}
    </Card>
  );
}
