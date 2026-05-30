import { Card } from "./ui";

interface AutogroupCardProps {
  sampleCount: number;
  orderingVariable: string | null;
  onAdjust: () => void;
}

/**
 * AutogroupCard — the confident-expert "Auto-grouped" note in the series
 * builder read rail (R8 / B-I). Mockup `series-builder.html` `.autogroup`:
 * a spark glyph + "Auto-grouped" title, a sentence naming how many samples
 * Himalaya read as one series and the ordering variable, and Confirm / Adjust
 * link buttons. Plate surface, hair border.
 */
export function AutogroupCard({
  sampleCount, orderingVariable, onAdjust,
}: AutogroupCardProps): JSX.Element {
  return (
    <Card
      data-testid="autogroup-card"
      className="p-3.5"
    >
      <div className="mb-1.5 flex items-center gap-1.5">
        <svg className="h-[15px] w-[15px] shrink-0" viewBox="0 0 16 16" fill="none" aria-hidden="true">
          <path
            d="M8 1.4l1.7 4.2 4.5.3-3.5 2.9 1.2 4.4L8 10.9 4.1 13.2l1.2-4.4L1.8 5.9l4.5-.3z"
            fill="var(--color-print-accent)"
          />
        </svg>
        <span className="text-xs font-bold text-ink">Auto-grouped</span>
      </div>
      <p className="text-xs leading-relaxed text-ink-soft">
        Himalaya read <b className="font-semibold text-ink">{sampleCount} samples</b> as one
        series from{" "}
        {orderingVariable !== null ? (
          <>
            their names, ordered by{" "}
            <b className="font-semibold text-ink">{orderingVariable}</b>.
          </>
        ) : (
          <>their names.</>
        )}
      </p>
      <div className="mt-2 flex gap-3">
        <button
          type="button"
          data-testid="autogroup-adjust"
          onClick={onAdjust}
          className="text-xs font-semibold text-print-accent hover:underline"
        >
          Adjust grouping
        </button>
      </div>
    </Card>
  );
}
