import { useState } from "react";
import type { MouseEvent, ReactNode } from "react";
import type { CollectionTimeline, ReviewProgress } from "../api";
import { clockTime, dayAndTime, describeRibbon, duration, reviewSummary, ribbonAt, ribbonLayout, toReview } from "../lib/collectionRibbon";
import type { RibbonPoint } from "../lib/collectionRibbon";
import { Tooltip } from "./Tooltip";

const VIEW_W = 1000;
const VIEW_H = 10;
const MIN_SPAN = 3;
const MIN_MARK = 6;

interface CollectionRibbonProps {
  timeline: CollectionTimeline | null | undefined;
  review: ReviewProgress | undefined;
  className?: string;
}

function plural(n: number, word: string): string {
  return `${n} ${word}${n === 1 ? "" : "s"}`;
}

function pointCopy(tl: CollectionTimeline, p: RibbonPoint): ReactNode {
  if (p.kind === "gap") return `No images for ${duration(p.to - p.from)}`;
  const when = p.to - p.from >= 60
    ? `${dayAndTime(tl.start, p.from)} to ${clockTime(tl.start, p.to)}`
    : dayAndTime(tl.start, p.to);
  return (
    <span className="flex flex-col">
      <span>{when}</span>
      <span>{p.reviewed ? "Reviewed" : "Not reviewed yet"}</span>
      {p.beamOff > 0 && <span className="text-warning">{plural(p.beamOff, "image")} with the beam off</span>}
    </span>
  );
}

function Caption({ left, review }: { left: string; review: ReviewProgress | undefined }): JSX.Element {
  return (
    <div className="mt-1 flex justify-between gap-2 whitespace-nowrap text-xs text-ink-faint">
      <span className="truncate">{left}</span>
      {review && <span>{toReview(review)}</span>}
    </div>
  );
}

export function CollectionRibbon({ timeline, review, className }: CollectionRibbonProps): JSX.Element {
  const [point, setPoint] = useState<RibbonPoint | null>(null);
  const caption = timeline ? describeRibbon(timeline) : null;

  if (!timeline || !caption) {
    return (
      <div className={className} data-testid="collection-ribbon" data-empty="true">
        <svg className="block h-3 w-full" viewBox={`0 0 ${VIEW_W} ${VIEW_H + 2}`} preserveAspectRatio="none" aria-hidden="true">
          <rect className="fill-paper-sunk" x={0} y={1} width={VIEW_W} height={VIEW_H} />
        </svg>
        <Caption left="No image times yet" review={review} />
      </div>
    );
  }

  const collected = `Collected from ${caption.started} over ${caption.length}.`;
  const { spans, end } = ribbonLayout(timeline);
  const total = Math.max(end, 1);
  const x = (t: number): number => (t / total) * VIEW_W;

  const onMove = (e: MouseEvent<SVGSVGElement>): void => {
    const r = e.currentTarget.getBoundingClientRect();
    const t = ((e.clientX - r.left) / Math.max(r.width, 1)) * total;
    setPoint(ribbonAt(timeline, t));
  };

  return (
    <div className={className} data-testid="collection-ribbon">
      <Tooltip content={point ? pointCopy(timeline, point) : collected}>
        <svg
          className="block h-3 w-full"
          viewBox={`0 0 ${VIEW_W} ${VIEW_H + 2}`}
          preserveAspectRatio="none"
          shapeRendering="crispEdges"
          role="img"
          aria-label={review ? `${reviewSummary(review)}. ${collected}` : collected}
          onMouseMove={onMove}
          onMouseLeave={() => setPoint(null)}
        >
          <rect className="fill-paper-sunk" x={0} y={1} width={VIEW_W} height={VIEW_H} />
          {spans.map((s, i) => (
            <rect
              key={i}
              data-reviewed={s.reviewed}
              className={s.reviewed ? "fill-success" : "fill-hair-strong"}
              x={x(s.from)}
              y={1}
              width={Math.max(x(s.to) - x(s.from), MIN_SPAN)}
              height={VIEW_H}
            />
          ))}
          {timeline.beam_off.map((t, i) => (
            <rect key={`off-${i}`} data-beam-off className="fill-warning" x={Math.max(x(t) - MIN_MARK, 0)} y={1} width={MIN_MARK} height={VIEW_H} />
          ))}
        </svg>
      </Tooltip>
      <Caption left={`${caption.started} · ${caption.length}`} review={review} />
    </div>
  );
}
