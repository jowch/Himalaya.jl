import { useCallback, useEffect, useRef, useState, type CSSProperties } from "react";
import { decideOrient } from "../lib/detectorOrient";

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

// The orient decision (thresholds, viewport gate, caps swap) lives in the
// shared `decideOrient` helper so DetectorRingOverlay (#180) can rotate in
// lockstep with this canvas.

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
  /**
   * R2-M14 (#207): contact-sheet thumbs only fetch their PNG once on-screen.
   * The corpus has 139 samples × ~4 exposures ≈ 556 thumbs; an eager fetch on
   * every row mount blew past Chromium's 6-per-host HTTP/1.1 limit on the
   * live surface (`ERR_INSUFFICIENT_RESOURCES`, rows 2-7+ stuck at "Loading
   * frames…"). The IntersectionObserver gate keeps the request volume to
   * what's actually visible during a scroll.
   *
   * The gate is skipped (auto-visible) for `size === "full"` — the loupe big
   * frame is always rendered above-the-fold, and gating the loupe would make
   * its bones visible to JSDOM (which has no IntersectionObserver), regressing
   * existing unit tests.
   *
   * IntersectionObserver is absent in JSDOM, so we degrade-to-eager when the
   * constructor is undefined — the existing thumb tests keep painting.
   */
  const [hasIntersected, setHasIntersected] = useState<boolean>(
    () =>
      size === "full" ||
      typeof window === "undefined" ||
      typeof window.IntersectionObserver !== "function",
  );

  const evaluateOrient = useCallback(() => {
    const wrapper = wrapperRef.current;
    const canvas = canvasRef.current;
    // Preserve prior layout on a transient zero-size — do NOT setLayout here
    // (a later observe corrects it). decideOrient's own zero-guard would
    // instead clobber to {portrait,null}, so the skip must stay in the caller.
    if (!wrapper || !canvas || !canvas.width || !canvas.height) return;
    // U-3 (#256): contact-sheet thumbs lock portrait so the sheet reads as a
    // uniform grid of windows — the eye compares the *content* of each window,
    // not its orientation. Only the loupe/focus `size="full"` frame rotates, to
    // maximise the data area on a wide canvas. The shared `decideOrient` helper
    // (also driving DetectorRingOverlay, #180) is intentionally NOT consulted
    // here; the overlay only renders on the `full` frame, so locking thumbs has
    // no cross-effect on ring rotation.
    if (size === "thumb") {
      setLayout({ orient: "portrait", caps: null });
      return;
    }
    const cw = wrapper.clientWidth;
    const ch = wrapper.clientHeight;
    if (cw === 0 || ch === 0) return;
    setLayout(decideOrient({
      containerW: cw,
      containerH: ch,
      imageW: canvas.width,
      imageH: canvas.height,
      viewportW: typeof window !== "undefined" ? window.innerWidth : 0,
    }));
  }, [size]);

  const renderImage = useCallback(async () => {
    const canvas = canvasRef.current;
    if (!canvas || !imagePath) return;
    // R2-M14: thumb scroll-gating. The component still mounts (so its
    // wrapper occupies layout space and the observer can detect it), but
    // we skip the actual fetch+decode work until the wrapper enters the
    // viewport.
    if (!hasIntersected) return;

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

    // U-1 (#255): warm the detector LUT. The endpoints are the detector-window
    // tokens, not the page bg/fg: intensity 0 (background / window backing) maps
    // to `frame-edge` (warm near-black, the dark window set into the paper) and
    // intensity 255 (bright signal) maps to `frame-signal` (warm off-white). This
    // is NON-inverting — the brightest pixels are the lightest output — so Bragg
    // rings read as "ink on dark paper" inside the warm window and bridge to the
    // surrounding `frame-edge` border, instead of the prior paper->ink ramp that
    // rendered the image cold (light) against the dark frame. The backend PNG
    // stays grayscale; only this display ramp moved (the chroma is applied here).
    const [er, eg, eb] = getCssColor("--color-frame-edge");
    const [sr, sg, sb] = getCssColor("--color-frame-signal");
    const data = imageData.data;
    for (let i = 0; i < data.length; i += 4) {
      const t = data[i] / 255;
      data[i]     = Math.round(er + t * (sr - er));
      data[i + 1] = Math.round(eg + t * (sg - eg));
      data[i + 2] = Math.round(eb + t * (sb - eb));
      data[i + 3] = 255;
    }

    canvas.width  = width;
    canvas.height = height;
    canvas.getContext("2d")?.putImageData(imageData, 0, 0);

    // New intrinsic dims may flip the orient decision.
    evaluateOrient();
  }, [exposureId, imagePath, size, hasIntersected, evaluateOrient]);

  useEffect(() => {
    renderImage();
  }, [renderImage]);

  // R2-M14: latch hasIntersected when the wrapper enters the viewport.
  // One-shot — once a thumb is decoded, we keep it (so scrolling back doesn't
  // re-fetch). A small rootMargin warms thumbs just outside the viewport so
  // the scroll feels live, not staggered.
  useEffect(() => {
    if (hasIntersected) return;
    const wrapper = wrapperRef.current;
    if (!wrapper) return;
    if (typeof window === "undefined" ||
        typeof window.IntersectionObserver !== "function") {
      setHasIntersected(true);
      return;
    }
    const io = new window.IntersectionObserver(
      (entries) => {
        for (const entry of entries) {
          if (entry.isIntersecting) {
            setHasIntersected(true);
            io.disconnect();
            break;
          }
        }
      },
      { rootMargin: "200px 0px" },
    );
    io.observe(wrapper);
    return () => io.disconnect();
  }, [hasIntersected]);

  // (R0c #223) The `<html>` theme-class MutationObserver was removed: R0a
  // retired the dark↔light theme toggle, so nothing ever mutates the class
  // and the observer never fired — dead scaffolding.

  // Watch wrapper size — rotate when container becomes much wider than image.
  useEffect(() => {
    const wrapper = wrapperRef.current;
    if (!wrapper) return;
    const ro = new ResizeObserver(() => evaluateOrient());
    ro.observe(wrapper);
    return () => ro.disconnect();
  }, [evaluateOrient]);

  if (!imagePath) {
    // U-2 / R3-S06 (#255): the missing-image placeholder is a `frame-edge`
    // window (matching the live detector treatment) with a `frame-tag` mono
    // caption — not light `text-fg-muted` text on paper. Empty states are where
    // The Print is most exposed, so the dark window persists when data is absent.
    return (
      <div
        data-testid="detector-image-placeholder"
        data-variant="frame-window"
        className={`flex items-center justify-center bg-frame-edge text-frame-tag font-mono text-[11px] tracking-wide ${className ?? ""}`}
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
