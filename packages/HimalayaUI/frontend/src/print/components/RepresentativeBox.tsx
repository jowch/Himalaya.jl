import { useLayoutEffect, useRef } from "react";
import type { ReactNode } from "react";
import { Button, Card, Dot } from "../ui";

export interface RepresentativeBoxProps {
  /** true = this exposure IS the representative. */
  isRepresentative: boolean;
  /** true = the SAMPLE's representative exposure (whichever frame that is)
   *  is currently dropped. The backend still carries it to the Index stage
   *  (status is never consulted there), so the box must say so honestly. */
  representativeDropped?: boolean;
  /** Set this exposure as representative. */
  onSetRepresentative?: () => void;
  /** Body explanation override. */
  children?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function RepresentativeBox({
  isRepresentative,
  representativeDropped = false,
  onSetRepresentative,
  children,
  className,
}: RepresentativeBoxProps): JSX.Element {
  // Focus retention: activating the set button flips isRepresentative
  // optimistically, unmounting the still-focused button (controls-don't-lie)
  // — without this, the DOM dumps focus on document.body and a keyboard
  // user loses their Tab position at the moment of success. The flag is set
  // ONLY by the button's own activation, so a flip arriving from elsewhere
  // (another user via SSE) never yanks focus. Card doesn't forward refs
  // (plain function component, React 18), so the focus target is an inner
  // wrapper div rather than the Card itself.
  const boxRef = useRef<HTMLDivElement>(null);
  const focusWasInside = useRef(false);
  const prevIsRep = useRef(isRepresentative);
  useLayoutEffect(() => {
    const was = prevIsRep.current;
    prevIsRep.current = isRepresentative;
    if (!was && isRepresentative && focusWasInside.current) {
      boxRef.current?.focus(); // before paint, so focus never flashes on body
    }
    focusWasInside.current = false; // consume — the flag lives one transition
  }, [isRepresentative]);

  const handleSetClick = onSetRepresentative
    ? () => {
        focusWasInside.current = true;
        onSetRepresentative();
      }
    : undefined;

  return (
    <Card
      selected={isRepresentative}
      padding="sm"
      {...(className ? { className } : {})}
    >
      <div ref={boxRef} tabIndex={-1} data-testid="representative-box">
        {isRepresentative && (
          <div className="flex items-center gap-1.5 text-print-accent font-bold text-meta">
            <Dot tone="accent" size="sm" aria-hidden />
            Representative for indexing
          </div>
        )}

        <p className="text-caption text-ink-soft mt-1">
          {children ??
            "One frame per sample carries forward to the Index stage. Pick the cleanest, strongest frame."}
        </p>

        {representativeDropped && (
          // Severity-as-word (the Toast pattern): bold leading "Warning." in
          // default ink tones — small warning-hued text fails AA on plate.
          <p data-testid="rep-dropped-warning" className="text-caption mt-1.5">
            <strong className="font-bold text-ink">Warning.</strong> The
            representative frame is dropped. It still carries to the Index
            stage. Restore it or set another frame as representative.
          </p>
        )}

        {/* Controls-don't-lie: on the representative itself the button would
            be a no-op, so it is OMITTED — the accent header already states it. */}
        {!isRepresentative && (
          <Button
            variant="outline"
            className="mt-2"
            {...(handleSetClick ? { onClick: handleSetClick } : {})}
          >
            Set as representative
          </Button>
        )}
      </div>
    </Card>
  );
}
