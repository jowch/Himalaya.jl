import { useCallback, useEffect, useRef, useState, type CSSProperties } from "react";
import { decideOrient } from "../../lib/detectorOrient";
import type { DetectorLutVariant } from "./detectorLut";

interface Props {
  /** Image URL to fetch + paint; null -> frame-window placeholder. Caller builds
   *  it (Storybook: a fixture asset URL; composite: the API image URL). */
  src: string | null;
  /** `thumb` locks portrait + scroll-gates the fetch; `full` rotates + eager. */
  size: "thumb" | "full";
  /** Colormap: "neutral" (default) keeps the image a warm-gray so the phase-ring
   *  overlay owns all the colour; "warm" is the saturated beauty ramp. Applied as
   *  a GPU SVG filter (DetectorLutFilters), so switching it never refetches. */
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
  const [layout, setLayout] = useState<Layout>({ orient: "portrait", caps: null });
  const [hasIntersected, setHasIntersected] = useState<boolean>(
    () => size === "full" || typeof window === "undefined" ||
      typeof window.IntersectionObserver !== "function",
  );
  // Object URL for the fetched blob. We fetch in JS (not a bare <img src>) ONLY
  // to read the X-Image-Width/Height calibration headers — a plain <img> exposes
  // just the capped naturalWidth, not the raw detector dims the q-rings need.
  // `objectUrlRef` mirrors the committed URL so the lifecycle (revoke the prior
  // one when a new blob lands, revoke the last on unmount) needs no state in the
  // fetch effect's deps.
  const [objectUrl, setObjectUrl] = useState<string | null>(null);
  const objectUrlRef = useRef<string | null>(null);
  // True while a fetch for the current `src` is in flight. The <img> keeps
  // showing the PRIOR blob until the new one lands, so without this the frame
  // silently displays a stale image on a sample/exposure switch. Drives the
  // dim + spinner "new image coming" cue below.
  const [loading, setLoading] = useState<boolean>(false);
  const onRawSizeRef = useRef(onRawSize); onRawSizeRef.current = onRawSize;
  const onOrientRef = useRef(onOrient); onOrientRef.current = onOrient;
  // Raw detector dims (from headers). Aspect is preserved by the uniform 1536
  // server cap, so these drive the portrait/landscape decision too.
  const rawSizeRef = useRef<{ w: number; h: number } | null>(null);
  const prevOrientRef = useRef<"portrait" | "landscape">("portrait");

  const evaluateOrient = useCallback(() => {
    const wrapper = wrapperRef.current;
    const rs = rawSizeRef.current;
    if (!wrapper || !rs || !rs.w || !rs.h) return;
    let next: Layout;
    if (size === "thumb") {
      next = { orient: "portrait", caps: null };
    } else {
      const cw = wrapper.clientWidth, ch = wrapper.clientHeight;
      if (cw === 0 || ch === 0) return;
      next = decideOrient({
        containerW: cw, containerH: ch, imageW: rs.w, imageH: rs.h,
        viewportW: typeof window !== "undefined" ? window.innerWidth : 0,
      });
    }
    setLayout(next);
    if (prevOrientRef.current !== next.orient) {
      prevOrientRef.current = next.orient;
      onOrientRef.current?.(next.orient);
    }
  }, [size]);

  // Fetch the bytes once: read calibration headers, then hand the blob to the
  // <img> as an object URL. Deliberately NOT keyed on `lutVariant` — the colormap
  // is a CSS/SVG filter applied at render, so a variant toggle re-styles without
  // re-downloading (the canvas path used to refetch on every toggle).
  useEffect(() => {
    if (!src || !hasIntersected) return;
    let cancelled = false;
    setLoading(true);
    void (async () => {
      try {
        const res = await fetch(src);
        if (cancelled || !res.ok) return;
        const rw = Number(res.headers?.get?.("X-Image-Width"));
        const rh = Number(res.headers?.get?.("X-Image-Height"));
        if (Number.isFinite(rw) && Number.isFinite(rh) && rw > 0 && rh > 0 &&
            (rawSizeRef.current?.w !== rw || rawSizeRef.current?.h !== rh)) {
          rawSizeRef.current = { w: rw, h: rh };
          onRawSizeRef.current?.(rw, rh);
        }
        const blob = await res.blob();
        // `cancelled` is re-checked after every await, so a URL is only created
        // for a still-live effect — never leaked past unmount / a src change.
        if (cancelled) return;
        const url = URL.createObjectURL(blob);
        // Swap: revoke the URL this <img> was showing before pointing it at the new
        // blob. One owner of the lifecycle (this ref), so no double-revoke.
        if (objectUrlRef.current) URL.revokeObjectURL(objectUrlRef.current);
        objectUrlRef.current = url;
        setObjectUrl(url);
        evaluateOrient();
      } finally {
        // Clear the cue on success OR error, but NOT on cancel — a cancelled run
        // means `src` changed and the next run already set loading=true, so the
        // spinner should stay up for that in-flight fetch instead of flickering.
        if (!cancelled) setLoading(false);
      }
    })();
    return () => { cancelled = true; };
  }, [src, hasIntersected, evaluateOrient]);

  // Revoke the last committed URL on unmount (the swap above handles replacement
  // while mounted). Frees the blob — createObjectURL pins it until revoked.
  useEffect(() => () => {
    if (objectUrlRef.current) { URL.revokeObjectURL(objectUrlRef.current); objectUrlRef.current = null; }
  }, []);

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

  const imgStyle: CSSProperties = {
    // GPU colormap (DetectorLutFilters must be mounted app-wide / in Storybook).
    filter: `url(#detector-lut-${lutVariant})`,
    imageRendering: "pixelated",
    ...(layout.orient === "landscape" && layout.caps
      ? { maxWidth: `${layout.caps.maxW}px`, maxHeight: `${layout.caps.maxH}px`,
          transform: "rotate(90deg)", transformOrigin: "center" }
      // Fill the frame (scale up AND down), aspect-preserved. maxWidth:100% alone
      // only ever scaled DOWN, so a sub-frame image (e.g. a low-res capture) drew
      // at its native pixel size and floated tiny in the box. The frame is the
      // canonical display size; the image fits it.
      : { width: "100%", height: "100%", objectFit: "contain" }),
    // Dim the outgoing (stale) image while its replacement loads, so the frame
    // reads as "new image coming" rather than silently showing the old one.
    opacity: loading ? 0.4 : 1,
    transition: "opacity 150ms ease",
  };

  return (
    <div ref={wrapperRef} data-orient={layout.orient}
         className={`relative flex items-center justify-center w-full h-full overflow-hidden ${className ?? ""}`}>
      <img
        {...(objectUrl ? { src: objectUrl } : {})}
        alt="Detector image"
        style={imgStyle}
      />
      {loading && (
        <div data-testid="detector-image-loading"
             className="absolute inset-0 flex items-center justify-center pointer-events-none">
          <div className={`${size === "thumb" ? "h-4 w-4" : "h-6 w-6"} rounded-full border-2 border-frame-tag/30 border-t-frame-tag animate-spin`} />
        </div>
      )}
    </div>
  );
}
