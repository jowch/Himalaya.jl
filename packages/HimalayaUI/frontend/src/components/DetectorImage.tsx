import { useCallback, useEffect, useRef, useState, type CSSProperties } from "react";

interface Props {
  exposureId: number;
  imagePath: string | null;
  /**
   * Cache-busting token from the backend (Exposure.image_version).
   * Appended to the URL as `?v=<token>` so the browser can cache the PNG
   * aggressively while still picking up real changes (TIFF mtime moved or
   * IMAGE_PROCESSING_VERSION bumped).
   */
  imageVersion: string;
  size: "thumb" | "full";
  className?: string;
}

function getCssColor(varName: string): [number, number, number] {
  const raw = getComputedStyle(document.documentElement)
    .getPropertyValue(varName)
    .trim();
  const c = document.createElement("canvas");
  c.width = c.height = 1;
  const ctx = c.getContext("2d");
  if (!ctx) return [0, 0, 0];
  ctx.fillStyle = raw;
  ctx.fillRect(0, 0, 1, 1);
  const d = ctx.getImageData(0, 0, 1, 1).data;
  return [d[0], d[1], d[2]];
}

interface Layout {
  orient: "portrait" | "landscape";
  caps: { maxW: number; maxH: number } | null;
}

// Rotate when the container is meaningfully wider than the image's natural
// aspect — keeps near-square images upright unless space really wants landscape.
const ROTATE_THRESHOLD = 1.25;

// Auto-rotate is gated to viewports ≥ this width. Below the WorkspaceGrid
// breakpoint the layout stacks into a single column where the image card is
// wider than tall by default — but at that breakpoint the user has chosen a
// portrait-friendly layout, so we keep the exposure upright and let the image
// fit by width. Matching the WorkspaceGrid breakpoint keeps the two decisions
// synchronized: one CSS source of truth for "small layout."
const ROTATE_MIN_VIEWPORT = 1400;

export function DetectorImage({
  exposureId,
  imagePath,
  imageVersion,
  size,
  className,
}: Props): JSX.Element {
  const wrapperRef = useRef<HTMLDivElement>(null);
  const canvasRef = useRef<HTMLCanvasElement>(null);
  const [layout, setLayout] = useState<Layout>({
    orient: "portrait",
    caps: null,
  });

  const evaluateOrient = useCallback(() => {
    const wrapper = wrapperRef.current;
    const canvas = canvasRef.current;
    if (!wrapper || !canvas || !canvas.width || !canvas.height) return;
    const cw = wrapper.clientWidth;
    const ch = wrapper.clientHeight;
    if (cw === 0 || ch === 0) return;
    // Below the small-screen breakpoint, force portrait — the stacked layout
    // gives the image card a wide-but-shallow slot, so rotating to landscape
    // would make the diffraction pattern read sideways for no real gain.
    const viewportW = typeof window !== "undefined" ? window.innerWidth : 0;
    if (viewportW < ROTATE_MIN_VIEWPORT) {
      setLayout({ orient: "portrait", caps: null });
      return;
    }
    const containerAspect = cw / ch;
    const imageAspect = canvas.width / canvas.height;
    if (containerAspect > imageAspect * ROTATE_THRESHOLD) {
      // Pre-rotation max-width must be capped by container HEIGHT (becomes
      // visual height after rotation), and max-height by container WIDTH.
      setLayout({ orient: "landscape", caps: { maxW: ch, maxH: cw } });
    } else {
      setLayout({ orient: "portrait", caps: null });
    }
  }, []);

  const renderImage = useCallback(async () => {
    const canvas = canvasRef.current;
    if (!canvas || !imagePath) return;

    // The `?v=<imageVersion>` token makes the URL unique per (TIFF mtime,
    // IMAGE_PROCESSING_VERSION) pair so the browser can cache the PNG
    // aggressively. No `cache: "no-store"` — the URL itself is the cache key.
    const params = new URLSearchParams();
    if (size === "thumb") params.set("thumb", "1");
    if (imageVersion) params.set("v", imageVersion);
    const qs = params.toString();
    const url = `/api/exposures/${exposureId}/image${qs ? `?${qs}` : ""}`;
    const res = await fetch(url);
    if (!res.ok) return;

    const blob = await res.blob();
    const bitmap = await createImageBitmap(blob);

    // Draw grayscale to offscreen canvas to read pixel data
    const { width, height } = bitmap;
    const off = new OffscreenCanvas(width, height);
    const offCtx = off.getContext("2d")!;
    offCtx.drawImage(bitmap, 0, 0);
    bitmap.close();
    const imageData = offCtx.getImageData(0, 0, width, height);

    // Build LUT: bg color (intensity 0) → fg color (intensity 255)
    const [br, bg, bb] = getCssColor("--color-bg");
    const [fr, fg, fb] = getCssColor("--color-fg");
    const data = imageData.data;
    for (let i = 0; i < data.length; i += 4) {
      const t = data[i] / 255;
      data[i]     = Math.round(br + t * (fr - br));
      data[i + 1] = Math.round(bg + t * (fg - bg));
      data[i + 2] = Math.round(bb + t * (fb - bb));
      data[i + 3] = 255;
    }

    canvas.width  = width;
    canvas.height = height;
    canvas.getContext("2d")?.putImageData(imageData, 0, 0);

    // New intrinsic dims may flip the orient decision.
    evaluateOrient();
  }, [exposureId, imagePath, size, evaluateOrient]);

  useEffect(() => {
    renderImage();
  }, [renderImage]);

  // Re-apply colormap when theme changes (AppShell toggles class on <html>)
  useEffect(() => {
    const observer = new MutationObserver(() => renderImage());
    observer.observe(document.documentElement, {
      attributes: true,
      attributeFilter: ["class"],
    });
    return () => observer.disconnect();
  }, [renderImage]);

  // Watch wrapper size — rotate when container becomes much wider than image.
  useEffect(() => {
    const wrapper = wrapperRef.current;
    if (!wrapper) return;
    const ro = new ResizeObserver(() => evaluateOrient());
    ro.observe(wrapper);
    return () => ro.disconnect();
  }, [evaluateOrient]);

  if (!imagePath) {
    return (
      <div
        data-testid="detector-image-placeholder"
        className={`flex items-center justify-center text-fg-muted text-xs ${className ?? ""}`}
      >
        No image
      </div>
    );
  }

  const canvasStyle: CSSProperties = {
    imageRendering: "pixelated",
    ...(layout.orient === "landscape" && layout.caps
      ? {
          maxWidth: `${layout.caps.maxW}px`,
          maxHeight: `${layout.caps.maxH}px`,
          transform: "rotate(90deg)",
          transformOrigin: "center",
        }
      : {
          maxWidth: "100%",
          maxHeight: "100%",
        }),
  };

  return (
    <div
      ref={wrapperRef}
      data-orient={layout.orient}
      className={`flex items-center justify-center w-full h-full overflow-hidden ${className ?? ""}`}
    >
      <canvas
        ref={canvasRef}
        role="img"
        aria-label="Detector image"
        style={canvasStyle}
      />
    </div>
  );
}
