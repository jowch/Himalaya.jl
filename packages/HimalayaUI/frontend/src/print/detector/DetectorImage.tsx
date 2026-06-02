import { useCallback, useEffect, useRef, useState, type CSSProperties } from "react";
import { decideOrient } from "../../lib/detectorOrient";
import { buildPrintDetectorLut } from "./detectorLut";

interface Props {
  /** Image URL to fetch + paint; null -> frame-window placeholder. Caller builds
   *  it (Storybook: a fixture asset URL; composite: the API image URL). */
  src: string | null;
  /** `thumb` locks portrait + scroll-gates the fetch; `full` rotates + eager. */
  size: "thumb" | "full";
  className?: string;
}

interface Layout { orient: "portrait" | "landscape"; caps: { maxW: number; maxH: number } | null }

export function DetectorImage({ src, size, className }: Props): JSX.Element {
  const wrapperRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [layout, setLayout] = useState<Layout>({ orient: "portrait", caps: null });
  const [hasIntersected, setHasIntersected] = useState<boolean>(
    () => size === "full" || typeof window === "undefined" ||
      typeof window.IntersectionObserver !== "function",
  );

  const evaluateOrient = useCallback(() => {
    const wrapper = wrapperRef.current, canvas = canvasRef.current;
    if (!wrapper || !canvas || !canvas.width || !canvas.height) return;
    if (size === "thumb") { setLayout({ orient: "portrait", caps: null }); return; }
    const cw = wrapper.clientWidth, ch = wrapper.clientHeight;
    if (cw === 0 || ch === 0) return;
    setLayout(decideOrient({
      containerW: cw, containerH: ch, imageW: canvas.width, imageH: canvas.height,
      viewportW: typeof window !== "undefined" ? window.innerWidth : 0,
    }));
  }, [size]);

  const renderImage = useCallback(async () => {
    const canvas = canvasRef.current;
    if (!canvas || !src || !hasIntersected) return;
    const res = await fetch(src);
    if (!res.ok) return;
    const bitmap = await createImageBitmap(await res.blob());
    const { width, height } = bitmap;
    const off = new OffscreenCanvas(width, height);
    const offCtx = off.getContext("2d")!;
    offCtx.drawImage(bitmap, 0, 0);
    bitmap.close();
    const imageData = offCtx.getImageData(0, 0, width, height);

    // Warm, NON-inverting Print LUT: index the 256-entry colormap by the
    // grayscale intensity (R channel of the log-normalized PNG).
    const lut = buildPrintDetectorLut();
    const data = imageData.data;
    for (let i = 0; i < data.length; i += 4) {
      const t = data[i];           // 0..255 intensity
      data[i]     = lut[t * 3];
      data[i + 1] = lut[t * 3 + 1];
      data[i + 2] = lut[t * 3 + 2];
      data[i + 3] = 255;
    }
    canvas.width = width;
    canvas.height = height;
    canvas.getContext("2d")?.putImageData(imageData, 0, 0);
    evaluateOrient();
  }, [src, hasIntersected, evaluateOrient]);

  useEffect(() => { renderImage(); }, [renderImage]);

  useEffect(() => {
    if (hasIntersected) return;
    const wrapper = wrapperRef.current;
    if (!wrapper) return;
    if (typeof window === "undefined" || typeof window.IntersectionObserver !== "function") {
      setHasIntersected(true); return;
    }
    const io = new window.IntersectionObserver((entries) => {
      for (const e of entries) if (e.isIntersecting) { setHasIntersected(true); io.disconnect(); break; }
    }, { rootMargin: "200px 0px" });
    io.observe(wrapper);
    return () => io.disconnect();
  }, [hasIntersected]);

  useEffect(() => {
    const wrapper = wrapperRef.current;
    if (!wrapper) return;
    const ro = new ResizeObserver(() => evaluateOrient());
    ro.observe(wrapper);
    return () => ro.disconnect();
  }, [evaluateOrient]);

  if (!src) {
    return (
      <div data-testid="detector-image-placeholder" data-variant="frame-window"
        className={`flex items-center justify-center bg-frame-edge text-frame-tag font-mono text-sm tracking-wide ${className ?? ""}`}>
        No image
      </div>
    );
  }

  const canvasStyle: CSSProperties = {
    imageRendering: "pixelated",
    ...(layout.orient === "landscape" && layout.caps
      ? { maxWidth: `${layout.caps.maxW}px`, maxHeight: `${layout.caps.maxH}px`,
          transform: "rotate(90deg)", transformOrigin: "center" }
      : { maxWidth: "100%", maxHeight: "100%" }),
  };

  return (
    <div ref={wrapperRef} data-orient={layout.orient}
         className={`flex items-center justify-center w-full h-full overflow-hidden ${className ?? ""}`}>
      <canvas ref={canvasRef} role="img" aria-label="Detector image" style={canvasStyle} />
    </div>
  );
}
