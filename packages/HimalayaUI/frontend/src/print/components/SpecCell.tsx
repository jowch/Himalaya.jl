import { CheckCircle } from "../ui/CheckCircle";

export interface SpecCellProps {
  name: string;
  sampleId: string;
  screened?: boolean;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function SpecCell({
  name,
  sampleId,
  screened,
  className,
}: SpecCellProps): JSX.Element {
  return (
    <div
      data-testid="spec-cell"
      className={`flex items-start gap-2.5${className ? ` ${className}` : ""}`}
    >
      <CheckCircle
        checked={!!screened}
        label={screened ? "Screened" : "Not screened"}
        className="mt-0.5"
      />
      <div className="min-w-0">
        <span data-role="spec-name" className="text-body font-semibold block leading-tight line-clamp-2">
          {name}
        </span>
        <span
          data-role="spec-id"
          className="text-data text-ink-faint block mt-0.5"
        >
          {sampleId}
        </span>
      </div>
    </div>
  );
}
