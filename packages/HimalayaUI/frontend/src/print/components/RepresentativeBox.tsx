import type { ReactNode } from "react";
import { Button, Card, Dot } from "../ui";

export interface RepresentativeBoxProps {
  /** true = this exposure IS the representative. */
  isRepresentative: boolean;
  /** Set this exposure as representative. */
  onSetRepresentative?: () => void;
  /** Body explanation override. */
  children?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function RepresentativeBox({
  isRepresentative,
  onSetRepresentative,
  children,
  className,
}: RepresentativeBoxProps): JSX.Element {
  return (
    <Card
      selected={isRepresentative}
      padding="sm"
      {...(className ? { className } : {})}
    >
      {isRepresentative && (
        <div className="flex items-center gap-1.5 text-print-accent font-bold text-meta">
          <Dot tone="accent" size="sm" aria-hidden />
          Representative for indexing
        </div>
      )}

      <p className="text-caption text-ink-soft mt-1">
        {children ??
          "One exposure per sample carries forward to the Index stage. Pick the cleanest, strongest frame."}
      </p>

      <Button
        variant="outline"
        className="mt-2"
        {...(onSetRepresentative ? { onClick: onSetRepresentative } : {})}
      >
        Set as representative
      </Button>
    </Card>
  );
}
