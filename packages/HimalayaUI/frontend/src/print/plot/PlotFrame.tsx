import { useEffect, useRef, useState, type ReactNode } from "react";

export interface Margins {
  top: number;
  right: number;
  bottom: number;
  left: number;
}

export interface PlotDims {
  width: number;
  height: number;
  plotWidth: number;
  plotHeight: number;
  margins: Margins;
}

export interface PlotFrameProps {
  height: number;
  margins: Margins;
  /** Controlled width; when omitted the frame measures its container. */
  width?: number;
  /** Width used before the first measurement (and in non-DOM tests). */
  defaultWidth?: number;
  className?: string;
  "data-testid"?: string;
  /** Container-relative pixel gestures. q-translation happens in the caller. */
  onWheelPx?: (deltaY: number, px: number, py: number) => void;
  onClickPx?: (px: number, py: number, altKey: boolean) => void;
  onDoubleClickPx?: () => void;
  /** Render the plot body, given the resolved pixel dimensions. */
  render: (dims: PlotDims) => ReactNode;
}

export function PlotFrame({
  height,
  margins,
  width,
  defaultWidth = 640,
  className,
  "data-testid": testid,
  onWheelPx,
  onClickPx,
  onDoubleClickPx,
  render,
}: PlotFrameProps): JSX.Element {
  const containerRef = useRef<HTMLDivElement | null>(null);
  const [measured, setMeasured] = useState<number | null>(null);

  // Measure container width unless an explicit width is supplied.
  useEffect(() => {
    if (width !== undefined) return;
    const el = containerRef.current;
    if (!el || typeof ResizeObserver === "undefined") return;
    const ro = new ResizeObserver((entries) => {
      const w = entries[0]?.contentRect.width ?? 0;
      if (w > 0) setMeasured(w);
    });
    ro.observe(el);
    return () => ro.disconnect();
  }, [width]);

  // React's onWheel is passive and cannot preventDefault(); bind non-passive.
  useEffect(() => {
    const el = containerRef.current;
    if (!el || !onWheelPx) return;
    function handle(ev: WheelEvent): void {
      ev.preventDefault();
      const rect = el!.getBoundingClientRect();
      onWheelPx!(ev.deltaY, ev.clientX - rect.left, ev.clientY - rect.top);
    }
    el.addEventListener("wheel", handle, { passive: false });
    return () => el.removeEventListener("wheel", handle);
  }, [onWheelPx]);

  const w = width ?? measured ?? defaultWidth;
  const plotWidth = Math.max(0, w - margins.left - margins.right);
  const plotHeight = Math.max(0, height - margins.top - margins.bottom);
  const dims: PlotDims = { width: w, height, plotWidth, plotHeight, margins };

  function handleClick(ev: React.MouseEvent): void {
    if (!onClickPx) return;
    const rect = containerRef.current!.getBoundingClientRect();
    onClickPx(ev.clientX - rect.left, ev.clientY - rect.top, ev.altKey);
  }

  return (
    <div
      ref={containerRef}
      className={className}
      style={{ width: width !== undefined ? width : "100%" }}
    >
      <svg
        width={w}
        height={height}
        viewBox={`0 0 ${w} ${height}`}
        data-testid={testid}
        onClick={handleClick}
        {...(onDoubleClickPx ? { onDoubleClick: onDoubleClickPx } : {})}
      >
        <g transform={`translate(${margins.left},${margins.top})`}>
          {render(dims)}
        </g>
      </svg>
    </div>
  );
}
