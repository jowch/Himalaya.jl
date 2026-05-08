/**
 * ExposureListRow (Plan §Phase 5, Task 5.1).
 *
 * Single row representing one exposure in the context of its parent sample.
 * Used by the picker modal (`ComparisonPicker`) and reserved for any future
 * Inspect surfaces that need a list-style exposure presentation. The Inspect
 * page itself currently uses `ThumbnailGallery` (image filmstrip) — it does
 * NOT consume this component yet because its current layout is image-driven,
 * not list-driven. The component is built standalone so the picker can use
 * it without forcing a presentation rewrite on Inspect.
 *
 * Action surface:
 *   - `onClick(exposureId)` — fires on row click. Suppressed when `locked`.
 *   - `checked` + `onCheckedChange` — when supplied, renders a checkbox.
 *   - `actionElement` — escape-hatch slot for any custom trailing content.
 *
 * Visuals are intentionally minimal — the picker scopes layout / hover / focus
 * styling at the list level. Notes truncation falls out of `line-clamp-2`
 * with the full text on the `title` attribute for native hover surfacing.
 */
import type { ReactNode, MouseEvent } from "react";
import type { Exposure, Sample } from "../api";

interface Props {
  exposure: Exposure;
  sample: Sample;
  /** Fires on row click. Receives the exposure id. Suppressed when `locked`. */
  onClick?: (exposureId: number) => void;
  /** Checkbox value. When supplied (defined), the row renders a checkbox slot. */
  checked?: boolean;
  /** Receives the new checkbox value. Required when `checked` is supplied. */
  onCheckedChange?: (next: boolean) => void;
  /** Lock the row — disables checkbox + suppresses onClick. */
  locked?: boolean;
  /**
   * Hint shown next to a locked checkbox (e.g. "already added"). Typed
   * `string | undefined` (not `string?:`) so callers can pass an undefined
   * literal under `exactOptionalPropertyTypes: true` — see the picker's
   * conditional `r.alreadyAdded ? "already added" : undefined`.
   */
  lockedReason?: string | undefined;
  /** Input control type. Defaults to "checkbox". */
  control?: "checkbox" | "radio";
  /** Escape-hatch trailing slot for callers that don't want checkbox/onClick. */
  actionElement?: ReactNode;
  /** Optional class merge for the row container. */
  className?: string;
}

function exposureLabel(exposure: Exposure): string {
  if (exposure.filename) return exposure.filename.replace(/\.dat$/i, "");
  return `#${exposure.id}`;
}

function sampleLabel(sample: Sample): string {
  // Spec: prefer sample.name, fall back to sample.label.
  return sample.name ?? sample.label ?? `Sample #${sample.id}`;
}

export function ExposureListRow({
  exposure,
  sample,
  onClick,
  checked,
  onCheckedChange,
  locked = false,
  lockedReason,
  control = "checkbox",
  actionElement,
  className,
}: Props): JSX.Element {
  const hasCheckbox = checked !== undefined;
  const handleRowClick = (e: MouseEvent<HTMLDivElement>): void => {
    if (locked) return;
    if (e.target instanceof HTMLInputElement) return; // checkbox handles itself
    if (onClick) onClick(exposure.id);
  };

  return (
    <div
      data-testid="exposure-list-row"
      data-exposure-id={String(exposure.id)}
      data-locked={locked ? "true" : undefined}
      onClick={handleRowClick}
      className={
        "flex items-start gap-3 px-3 py-2 rounded " +
        (locked
          ? "opacity-60 cursor-not-allowed"
          : "cursor-pointer hover:bg-bg-elevated") +
        (className ? ` ${className}` : "")
      }
    >
      {hasCheckbox && (
        <input
          type={control}
          data-testid="exposure-list-row-checkbox"
          checked={checked}
          disabled={locked}
          onChange={(e) => onCheckedChange?.(e.target.checked)}
          // stopPropagation so a click on the checkbox doesn't also bubble
          // up and try to invoke onClick (some callers wire both surfaces).
          onClick={(e) => e.stopPropagation()}
          className="mt-1 shrink-0"
        />
      )}
      <div className="flex-1 min-w-0">
        <div className="flex items-baseline gap-2">
          <span className="font-medium text-fg truncate">
            {exposureLabel(exposure)}
          </span>
          <span className="text-fg-muted text-sm truncate">
            {sampleLabel(sample)}
          </span>
          {locked && lockedReason && (
            <span className="text-xs text-fg-dim italic">
              {lockedReason}
            </span>
          )}
        </div>
        {sample.notes && (
          <div
            data-testid="exposure-list-row-notes"
            // line-clamp-2 truncates visually; title surfaces full text on hover.
            className="text-xs text-fg-muted line-clamp-2"
            title={sample.notes}
          >
            {sample.notes}
          </div>
        )}
      </div>
      {actionElement && <div className="shrink-0">{actionElement}</div>}
    </div>
  );
}
