import type { ReactNode } from "react";
import { Card, SegmentedControl, Swatch } from "../ui";
import { PlateHeader } from "./PlateHeader";
import { ToolBar } from "./ToolBar";
import { WaterfallChart } from "../waterfall";
import type { WaterfallRow } from "../waterfall/waterfallModel";

function cx(...parts: Array<string | false | null | undefined>): string {
  return parts.filter(Boolean).join(" ");
}

export type SeriesScale = "log" | "lin";

export interface SeriesPlateProps {
  kicker?: ReactNode;
  title: ReactNode;
  subtitle?: ReactNode;
  /** Member rows low→high (rendered bottom-up by WaterfallChart). */
  rows: WaterfallRow[];
  /** Global waterfall trace-offset; scales inter-trace separation (1 = unchanged). */
  offsetScale?: number;
  scale: SeriesScale;
  onScaleChange: (next: SeriesScale) => void;
  /** Controlled hot row + cursor q, lifted to the page (synced with MemberList). */
  hoveredKey?: string;
  onHoverRow?: (key?: string) => void;
  hoveredQ?: number;
  onHoverQ?: (q?: number) => void;
  /** Foot legend phase dots (e.g. ["Pn3m","Lamellar"]). */
  legendPhases?: string[];
  /** Foot hint (e.g. "peaks are light anchors — hover a trace to read its indexing"). */
  footHint?: ReactNode;
  /** Mono right-aligned foot note (e.g. "offset ×1.0 · log I"). */
  footNote?: ReactNode;
  /** PLACEMENT-ONLY. */
  className?: string;
}

export function SeriesPlate({
  kicker, title, subtitle, rows, offsetScale, scale, onScaleChange,
  hoveredKey, onHoverRow, hoveredQ, onHoverQ,
  legendPhases, footHint, footNote, className,
}: SeriesPlateProps): JSX.Element {
  const hasFoot = (legendPhases && legendPhases.length > 0) || footHint != null || footNote != null;
  return (
    <Card as="section" elevated data-testid="series-plate" className={cx("px-6 pt-5 pb-[18px]", className)}>
      <PlateHeader kicker={kicker} title={title} subtitle={subtitle} as="h1">
        <ToolBar>
          <SegmentedControl
            options={[{ value: "log", label: "log q" }, { value: "lin", label: "linear q" }]}
            value={scale}
            onChange={onScaleChange}
            aria-label="q scale"
          />
        </ToolBar>
      </PlateHeader>
      <WaterfallChart
        rows={rows}
        xType={scale === "log" ? "log" : "linear"}
        className="mt-2"
        {...(offsetScale !== undefined ? { offsetScale } : {})}
        {...(hoveredKey !== undefined ? { hoveredKey } : {})}
        {...(onHoverRow ? { onHoverRow } : {})}
        {...(hoveredQ !== undefined ? { hoveredQ } : {})}
        {...(onHoverQ ? { onHoverQ } : {})}
      />
      {hasFoot && (
        <div data-testid="series-plate-foot" className="mt-3 pt-[11px] border-t border-hair flex items-center justify-between text-meta text-ink-faint">
          <div className="flex items-center gap-[14px]">
            {legendPhases?.map((p) => (
              <span key={p} className="inline-flex items-center gap-[5px]">
                <Swatch phase={p} shape="circle" />
                <span>{p}</span>
              </span>
            ))}
            {footHint != null && <span>{footHint}</span>}
          </div>
          {footNote != null && <span className="text-data">{footNote}</span>}
        </div>
      )}
    </Card>
  );
}
