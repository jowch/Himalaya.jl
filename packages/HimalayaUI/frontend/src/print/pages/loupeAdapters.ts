import type { Exposure, SampleTag } from "../../api";
import type { GalleryExposure } from "../components/ThumbnailGallery";
import type { Tag } from "../ui";
import { representativeExposure } from "../../lib/exposure";

/** Loupe-local tag view-model: the ui `Tag` ({key, value?}) plus the backend
 *  identity (`id`) and provenance (`source`). The `id` is what makes two
 *  byte-identical pills (manual dose=10 + scoping dose=10) individually
 *  addressable for removal/edit — `ui/tag.ts` is deliberately NOT widened. */
export interface LoupeTag extends Tag {
  id: number;
  source: string;
}

/**
 * Detector image URL. Mirrors the legacy DetectorImage builder exactly so the
 * browser cache key matches (cache-coherence): `?thumb=1` for the strip,
 * `?v=<image_version>` only when present. null image_path → null (placeholder).
 */
export function buildExposureImageUrl(
  exposure: Exposure,
  opts?: { thumb?: boolean },
): string | null {
  if (exposure.image_path === null) return null;
  const params = new URLSearchParams();
  if (opts?.thumb) params.set("thumb", "1");
  if (exposure.image_version) params.set("v", exposure.image_version);
  const qs = params.toString();
  return `/api/exposures/${exposure.id}/image${qs ? `?${qs}` : ""}`;
}

/** Map the per-sample exposure list to the filmstrip view-model. Status is
 *  binary: a frame is either dropped (rejected) or kept (the default), so only
 *  the `rejected` channel is carried — kept needs no marker. */
export function toGalleryExposures(exposures: Exposure[]): GalleryExposure[] {
  const representative = representativeExposure(exposures)?.id;
  return exposures.map((e, i) => ({
    id: e.id,
    src: buildExposureImageUrl(e, { thumb: true }),
    frameNo: i + 1,
    representative: e.id === representative,
    rejected: e.status === "rejected",
  }));
}

/** SampleTag[] → LoupeTag[]; omit empty value (exactOptionalPropertyTypes).
 *  Each tag carries its backend `id` and `source` so two byte-identical pills
 *  (e.g. manual dose=10 + scoping dose=10) are individually addressable for
 *  removal — see LO-TAGDUP. */
export function toLoupeTags(tags: SampleTag[]): LoupeTag[] {
  return tags.map((t): LoupeTag =>
    t.value
      ? { id: t.id, key: t.key, value: t.value, source: t.source }
      : { id: t.id, key: t.key, source: t.source },
  );
}
