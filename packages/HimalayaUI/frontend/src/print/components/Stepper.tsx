import type { ReactNode } from "react";
import { IconButton } from "../ui";

export interface StepperProps {
  /** Primary line — the sample name/title. */
  label: ReactNode;
  /** Subtitle line, e.g. "sample 4 of 9". Optional. */
  position?: ReactNode;
  onPrev?: () => void;
  onNext?: () => void;
  /** Disable prev (at first sample). */
  prevDisabled?: boolean;
  /** Disable next (at last sample). */
  nextDisabled?: boolean;
  /** PLACEMENT-ONLY (margin/flex position). */
  className?: string;
}

export function Stepper({
  label,
  position,
  onPrev,
  onNext,
  prevDisabled,
  nextDisabled,
  className,
}: StepperProps): JSX.Element {
  return (
    <div className={`flex items-center gap-0.5${className ? ` ${className}` : ""}`}>
      <IconButton
        label="previous sample"
        disabled={prevDisabled}
        {...(onPrev ? { onClick: onPrev } : {})}
      >
        ‹
      </IconButton>

      <span className="min-w-48 text-center leading-none px-1.5">
        <span className="block text-meta font-semibold text-ink">
          {label}
        </span>
        {position && (
          <span className="block text-caption uppercase tracking-wide text-ink-faint mt-0.5">
            {position}
          </span>
        )}
      </span>

      <IconButton
        label="next sample"
        disabled={nextDisabled}
        {...(onNext ? { onClick: onNext } : {})}
      >
        ›
      </IconButton>
    </div>
  );
}
