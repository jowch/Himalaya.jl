import { useCallback, useEffect, useRef, useState, type CSSProperties } from "react";
import { decideOrient } from "../../lib/detectorOrient";
import { buildPrintDetectorLut, type DetectorLutVariant } from "./detectorLut";

interface Props {
  /** Image URL to fetch + paint; null -> frame-window placeholder. Caller builds
   *  it (Storybook: a fixture asset URL; composite: the API image URL). */
  src: string | null;
  /** `thumb` locks portrait + scroll-gates the fetch; `full` rotates + eager. */
  size: "thumb" | "full";
  /** Colormap: "neutral" (default) keeps the image a warm-gray so the phase-ring
   *  overlay owns all the colour; "warm" is the saturated beauty ramp. */
  lutVariant?: DetectorLutVariant;
  className?: string;
  /** Raw detector pixel size from X-Image-Width/Height headers (notify-only). */
  onRawSize?: (w: number, h: number) => void;
  /** Display orientation, fired only on a true transition (notify-only). */
  onOrient?: (o: "portrait" | "landscape") => void;
}

interface Layout { orient: "portrait" | "landscape"; caps: { maxW: number; maxH: number } | null }

export function DetectorImage({ src, size, lutVariant = "neutral", className, onRawSize, onOrient }: Props): JSX.Element {
  const wrapperRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [layout, setLayout] = useState<Layout>({ orient: "portrait", caps: null });
  const [hasIntersected, setHasIntersected] = useState<boolean>(
    () => size === "full" || typeof window === "undefined" ||
      typeof window.IntersectionObserver !== "function",
  );
  const onRawSizeRef = useRef(onRawSize); onRawSizeRef.current = onRawSize;
  const onOrientRef = useRef(onOrient); onOrientRef.current = onOrient;
  const rawSizeRef = useRef<{ w: number; h: number } | null>(null);
  const prevOrientRef = useRef<"portrait" | "landscape">("portrait");

  const evaluateOrient = useCallback(() => {
    const wrapper = wrapperRef.current, canvas = canvasRef.current;
    if (!wrapper || !canvas || !canvas.width || !canvas.height) return;
    let next: Layout;
    if (size === "thumb") {
      next = { orient: "portrait", caps: null };
    } else {
      const cw = wrapper.clientWidth, ch = wrapper.clientHeight;
      if (cw === 0 || ch === 0) return;
      next = decideOrient({
        containerW: cw, containerH: ch, imageW: canvas.width, imageH: canvas.height,
        viewportW: typeof window !== "undefined" ? window.innerWidth : 0,
      });
    }
    setLayout(next);
    if (prevOrientRef.current !== next.orient) {
      prevOrientRef.current = next.orient;
      onOrientRef.current?.(next.orient);
    }
  }, [size]);

  const renderImage = useCallback(async () => {
    const canvas = canvasRef.current;
    if (!canvas || !src || !hasIntersected) return;
    const res = await fetch(src);
    if (!res.ok) return;
    const rw = Number(res.headers?.get?.("X-Image-Width"));
    const rh = Number(res.headers?.get?.("X-Image-Height"));
    if (Number.isFinite(rw) && Number.isFinite(rh) && rw > 0 && rh > 0 &&
        (rawSizeRef.current?.w !== rw || rawSizeRef.current?.h !== rh)) {
      rawSizeRef.current = { w: rw, h: rh };
      onRawSizeRef.current?.(rw, rh);
    }
    const bitmap = await createImageBitmap(await res.blob());
    const { width, height } = bitmap;
    const off = new OffscreenCanvas(width, height);
    const offCtx = off.getContext("2d")!;
    offCtx.drawImage(bitmap, 0, 0);
    bitmap.close();
    const imageData = offCtx.getImageData(0, 0, width, height);

    // NON-inverting Print LUT: index the 256-entry colormap by the grayscale
    // intensity (R channel of the log-normalized PNG). Neutral warm-gray by
    // default so the phase-ring overlay carries the colour.
    const lut = buildPrintDetectorLut(lutVariant);
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
  }, [src, hasIntersected, lutVariant, evaluateOrient]);

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
      // Fill the frame (scale up AND down), aspect-preserved. maxWidth:100% alone
      // only ever scaled DOWN, so a sub-frame image (e.g. a low-res capture) drew
      // at its native pixel size and floated tiny in the box. The frame is the
      // canonical display size; the image fits it.
      : { width: "100%", height: "100%", objectFit: "contain" }),
  };

  return (
    <div ref={wrapperRef} data-orient={layout.orient}
         className={`flex items-center justify-center w-full h-full overflow-hidden ${className ?? ""}`}>
      <canvas ref={canvasRef} role="img" aria-label="Detector image" style={canvasStyle} />
    </div>
  );
}
