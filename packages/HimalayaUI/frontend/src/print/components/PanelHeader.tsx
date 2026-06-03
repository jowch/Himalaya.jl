// src/print/components/PanelHeader.tsx
import type { ReactNode } from "react";
import { Kicker } from "../ui";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export interface PanelHeaderProps {
  /** The uppercase section label (.panel-h), e.g. "Detector image". */
  label: ReactNode;
  /** Optional right-side tools slot (exposure strip / comb-view toggle). */
  children?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function PanelHeader({ label, children, className }: PanelHeaderProps): JSX.Element {
  return (
    <div
      data-testid="panel-header"
      className={cx("flex items-center justify-between gap-2.5 mb-3", className)}
    >
      <Kicker tone="faint" className="flex-shrink-0">{label}</Kicker>
      {children != null && (
        // min-w-0 lets a scrollable tools child (e.g. an overflow-x-auto
        // ThumbnailGallery) shrink + scroll IN PLACE rather than widen the header.
        <div data-testid="panel-header-tools" className="min-w-0">
          {children}
        </div>
      )}
    </div>
  );
}
