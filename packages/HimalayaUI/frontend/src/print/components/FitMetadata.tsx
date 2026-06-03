function cx(...p: Array<string | false | null | undefined>): string {
  return p.filter(Boolean).join(" ");
}

export interface FitMetadataProps {
  landed: number;
  total: number;
  paramName: string;
  paramValue: string;
  unit?: string;
  snapped?: boolean;
  className?: string;
}

export function FitMetadata({
  landed,
  total,
  paramName,
  paramValue,
  unit = "Å",
  snapped = false,
  className,
}: FitMetadataProps): JSX.Element {
  return (
    <div data-testid="fit-metadata" className={cx("text-meta text-ink-soft", className)}>
      Lands on <b className="text-ink font-bold">{landed}</b> of {total} observed peaks
      {" · "}
      <b className="text-ink font-bold">{paramName}</b> ={" "}
      <b className="text-ink font-bold">
        {paramValue} {unit}
      </b>
      {snapped && (
        <>
          {" · "}
          <span data-testid="fit-snapped" className="text-print-accent">
            snapped
          </span>
        </>
      )}
    </div>
  );
}
