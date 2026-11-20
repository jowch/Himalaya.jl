import { useEffect, useRef, useState } from "react";
import type { Peak, Trace } from "../../api";
import { traceGeometry, type CombLine } from "../plot";

interface Props {
  trace: Trace | undefined;
  peaks: readonly Peak[];
  combs: readonly CombLine[];
  height?: number;
}

function cssVar(el: Element, name: string): string {
  return getComputedStyle(el).getPropertyValue(name).trim() || "currentColor";
}

export function TraceCanvas({ trace, peaks, combs, height = 360 }: Props): JSX.Element {
  const wrapRef = useRef<HTMLDivElement | null>(null);
  const canvasRef = useRef<HTMLCanvasElement | null>(null);
  const [width, setWidth] = useState(0);

  useEffect(() => {
    const el = wrapRef.current;
    if (!el) return;
    const ro = new ResizeObserver((entries) => setWidth(Math.floor(entries[0]!.contentRect.width)));
    ro.observe(el);
    return () => ro.disconnect();
  }, []);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || width === 0) return;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr;
    canvas.height = height * dpr;
    const ctx = canvas.getContext("2d");
    if (!ctx) return;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, width, height);
    const g = traceGeometry(trace, peaks, combs, width, height);
    if (!g) return;
    const ink = cssVar(canvas, "--color-ink");
    const soft = cssVar(canvas, "--color-ink-faint");
    const hair = cssVar(canvas, "--color-hair");
    const accent = cssVar(canvas, "--color-accent");
    ctx.font = `12px ${cssVar(canvas, "--font-sans")}`;

    ctx.strokeStyle = hair;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(g.plot.left, g.plot.bottom + 0.5);
    ctx.lineTo(g.plot.right, g.plot.bottom + 0.5);
    ctx.stroke();
    ctx.fillStyle = soft;
    ctx.textAlign = "center";
    ctx.textBaseline = "top";
    for (const t of g.xTicks) ctx.fillText(t.label, t.x, g.plot.bottom + 6);
    ctx.textAlign = "right";
    ctx.fillText("q (Å⁻¹)", g.plot.right, g.plot.bottom + 22);

    for (const c of g.combs) {
      ctx.strokeStyle = c.color;
      ctx.globalAlpha = c.emphasis ? 0.9 : 0.45;
      ctx.setLineDash(c.emphasis ? [] : [3, 3]);
      for (const x of c.xs) {
        ctx.beginPath();
        ctx.moveTo(x + 0.5, g.plot.top);
        ctx.lineTo(x + 0.5, g.plot.bottom);
        ctx.stroke();
      }
    }
    ctx.setLineDash([]);
    ctx.globalAlpha = 1;

    ctx.strokeStyle = ink;
    ctx.lineWidth = 1.25;
    ctx.beginPath();
    g.path.forEach(([x, y], i) => (i === 0 ? ctx.moveTo(x, y) : ctx.lineTo(x, y)));
    ctx.stroke();

    for (const p of g.peaks) {
      ctx.beginPath();
      ctx.arc(p.x, p.y - 7, 3, 0, Math.PI * 2);
      if (p.excluded) {
        ctx.strokeStyle = soft;
        ctx.lineWidth = 1;
        ctx.stroke();
      } else {
        ctx.fillStyle = accent;
        ctx.fill();
      }
    }
  }, [trace, peaks, combs, width, height]);

  return (
    <div ref={wrapRef} className="w-full" style={{ height }}>
      <canvas
        ref={canvasRef}
        role="img"
        aria-label={trace ? `Trace with ${peaks.filter((p) => !p.excluded).length} peaks` : "No trace"}
        data-testid="workspace-trace"
        style={{ width, height }}
      />
    </div>
  );
}
