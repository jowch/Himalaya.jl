import type { Exposure } from "../api";

type Counted = Pick<Exposure, "frames_total" | "frames_rejected">;

function imagesKept(e: Counted): { kept: number; total: number } | null {
  if (e.frames_total == null) return null;
  return { kept: e.frames_total - (e.frames_rejected ?? 0), total: e.frames_total };
}

export function exposureKindLabel(e: Counted & Pick<Exposure, "kind">): string {
  const images = imagesKept(e);
  if (images) {
    if (images.total === 1 && images.kept === 1) return "1 image";
    return images.kept === images.total ? `${images.total} averaged` : `${images.kept} of ${images.total} averaged`;
  }
  if (e.kind === "averaged") return "Averaged";
  if (e.kind === "background_subtracted") return "Subtracted";
  return "Single";
}

export function exposureTimeLabel(e: Counted & Pick<Exposure, "exposure_time">): string {
  if (e.exposure_time == null) return "—";
  const perImage = `${e.exposure_time.toFixed(1)} s`;
  const images = imagesKept(e);
  return images ? `${images.kept} × ${perImage}` : perImage;
}
