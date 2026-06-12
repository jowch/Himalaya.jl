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
  /** Accessible role for the svg (e.g. "img"), so a figure svg announces as a
   *  named graphic in the a11y tree instead of a nameless img (WCAG 1.1.1). */
  svgRole?: string;
  /** Accessible name paired with `svgRole`. */
  svgLabel?: string;
  /** Container-relative pixel gestures. q-translation happens in the caller. */
  onWheelPx?: (deltaY: number, px: number, py: number) => void;
  onClickPx?: (px: number, py: number, altKey: boolean) => void;
  onDoubleClickPx?: () => void;
  /** Container-relative pointer position (same rect math as wheel/click). */
  onPointerMovePx?: (px: number, py: number) => void;
  onPointerLeave?: () => void;
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
  svgRole,
  svgLabel,
  onWheelPx,
  onClickPx,
  onDoubleClickPx,
  onPointerMovePx,
  onPointerLeave,
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
    const el = containerRef.current;
    if (!onClickPx || !el) return;
    const rect = el.getBoundingClientRect();
    onClickPx(ev.clientX - rect.left, ev.clientY - rect.top, ev.altKey);
  }

  function handlePointerMove(ev: React.PointerEvent): void {
    const el = containerRef.current;
    if (!onPointerMovePx || !el) return;
    const rect = el.getBoundingClientRect();
    onPointerMovePx(ev.clientX - rect.left, ev.clientY - rect.top);
  }

  return (
    <div
      ref={containerRef}
      className={className}
      style={{ width: width !== undefined ? width : "100%" }}
      onPointerMove={onPointerMovePx ? handlePointerMove : undefined}
      onPointerLeave={onPointerLeave}
    >
      <svg
        width={w}
        height={height}
        viewBox={`0 0 ${w} ${height}`}
        data-testid={testid}
        {...(svgRole ? { role: svgRole } : {})}
        {...(svgLabel ? { "aria-label": svgLabel } : {})}
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
