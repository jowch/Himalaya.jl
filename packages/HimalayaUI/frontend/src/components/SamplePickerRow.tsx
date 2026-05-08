import { useState } from "react";
import type { PickerSampleRow as PickerRow } from "../api";

interface Props {
  row: PickerRow;
  /** Whole-row checked state (the sample is in the picks set). */
  checked: boolean;
  onCheckedChange: (next: boolean) => void;
  /** Override exposure id, or undefined for default. */
  overrideExposureId: number | undefined;
  onOverrideChange: (exposureId: number) => void;
  /** When true, the row is locked because this exposure is already a draft member. */
  alreadyAdded?: boolean;
}

function sampleLabel(s: PickerRow["sample"]): string {
  return s.name ?? s.label ?? `Sample #${s.id}`;
}

export function SamplePickerRow({
  row, checked, onCheckedChange, overrideExposureId, onOverrideChange,
  alreadyAdded = false,
}: Props): JSX.Element {
  const [expanded, setExpanded] = useState(false);
  const disabled = row.indexing_exposure_id === null;
  const resolvedExposureId = overrideExposureId ?? row.indexing_exposure_id;

  const dataAttrs: Record<string, string> = {};
  dataAttrs["data-sample-id"] = String(row.sample.id);
  if (disabled) dataAttrs["data-disabled"] = "true";
  if (alreadyAdded) dataAttrs["data-locked"] = "true";
  if (resolvedExposureId !== null && resolvedExposureId !== undefined) {
    dataAttrs["data-exposure-id"] = String(resolvedExposureId);
  }

  return (
    <div
      data-testid="sample-picker-row"
      {...dataAttrs}
      className={
        "flex items-start gap-3 px-3 py-2 rounded " +
        (disabled || alreadyAdded
          ? "opacity-60 cursor-not-allowed"
          : "cursor-pointer hover:bg-bg-elevated") +
        (checked && !disabled && !alreadyAdded ? " ring-1 ring-accent/30" : "")
      }
    >
      {!disabled && (
        <input
          type="checkbox"
          data-testid="sample-picker-row-checkbox"
          checked={alreadyAdded ? true : checked}
          disabled={alreadyAdded}
          onChange={alreadyAdded ? undefined : (e) => onCheckedChange(e.target.checked)}
          className="mt-1 shrink-0"
        />
      )}
      <div className="flex-1 min-w-0">
        <div className="flex items-baseline gap-2">
          <span className="font-medium text-fg truncate">
            {sampleLabel(row.sample)}
          </span>
          {alreadyAdded && (
            <span className="text-xs text-fg-dim italic shrink-0">already added</span>
          )}
          <span className="text-xs text-fg-dim shrink-0">
            {row.all_exposures.length} {row.all_exposures.length === 1 ? "exposure" : "exposures"}
          </span>
          {!disabled && (
            <button
              type="button"
              data-testid="sample-picker-row-caret"
              aria-expanded={expanded}
              aria-label="Show exposures"
              onClick={(e) => { e.stopPropagation(); setExpanded((v) => !v); }}
              className="ml-auto text-fg-muted hover:text-fg text-xs px-1"
            >
              {expanded ? "▾" : "▸"}
            </button>
          )}
        </div>
        {row.sample.notes && (
          <div className="text-xs text-fg-muted line-clamp-2" title={row.sample.notes}>
            {row.sample.notes}
          </div>
        )}
        {expanded && (
          <ul role="radiogroup" aria-label="Override exposure"
              className="mt-2 ml-4 border-l border-border pl-2 flex flex-col">
            {row.all_exposures.map((e) => {
              const radioId = `sprow-radio-${e.id}`;
              return (
                <li key={e.id}>
                  <label
                    htmlFor={radioId}
                    className="flex items-center gap-2 px-2 py-1 rounded cursor-pointer hover:bg-bg-elevated"
                  >
                    <input
                      id={radioId}
                      type="radio"
                      checked={resolvedExposureId === e.id}
                      onChange={() => onOverrideChange(e.id)}
                      className="shrink-0"
                    />
                    <span className="text-sm text-fg truncate">
                      {e.filename ?? `#${e.id}`}
                    </span>
                  </label>
                </li>
              );
            })}
          </ul>
        )}
      </div>
    </div>
  );
}
