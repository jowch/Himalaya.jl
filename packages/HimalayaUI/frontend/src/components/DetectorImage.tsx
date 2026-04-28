import { useCallback, useEffect, useRef } from "react";

interface Props {
  exposureId: number;
  imagePath: string | null;
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

export function DetectorImage({
  exposureId,
  imagePath,
  size,
  className,
}: Props): JSX.Element {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  const renderImage = useCallback(async () => {
    const canvas = canvasRef.current;
    if (!canvas || !imagePath) return;

    const url = `/api/exposures/${exposureId}/image${size === "thumb" ? "?thumb=1" : ""}`;
    const res = await fetch(url, { cache: "no-store" });
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
  }, [exposureId, imagePath, size]);

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

  return (
    <canvas
      ref={canvasRef}
      role="img"
      aria-label="Detector image"
      className={`object-contain ${className ?? ""}`}
      style={{ imageRendering: "pixelated" }}
    />
  );
}
