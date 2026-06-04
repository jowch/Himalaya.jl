import type { Exposure, SampleTag } from "../../api";
import type { GalleryExposure } from "../components/ThumbnailGallery";
import type { MetaEntry, Tag } from "../ui";

/** Default exposure when the loupe opens: representative → first accepted → first. */
export function defaultExposureId(exposures: Exposure[]): number | undefined {
  const representative = exposures.find((e) => e.selected);
  if (representative) return representative.id;
  const firstAccepted = exposures.find((e) => e.status === "accepted");
  if (firstAccepted) return firstAccepted.id;
  return exposures[0]?.id;
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

/** Map the per-sample exposure list to the filmstrip view-model. */
export function toGalleryExposures(exposures: Exposure[]): GalleryExposure[] {
  return exposures.map((e, i) => ({
    id: e.id,
    src: buildExposureImageUrl(e, { thumb: true }),
    frameNo: i + 1,
    representative: e.selected,
    rejected: e.status === "rejected",
  }));
}

/**
 * "This exposure" metadata rows. `integration`/`collected` are placeholders
 * until the backend lands those fields (#256). The legacy signal-meter row is
 * intentionally dropped — the greenfield LoupeSidePanel has no signal block.
 */
export function toMetaEntries(active: Exposure, exposures: Exposure[]): MetaEntry[] {
  const idx = exposures.findIndex((e) => e.id === active.id);
  const position = idx >= 0 ? `${idx + 1} of ${exposures.length}` : "—";
  return [
    { key: "frame", value: position },
    { key: "integration", value: "—" },
    { key: "collected", value: "—" },
  ];
}

/** SampleTag[] → greenfield Tag[]; omit empty value (exactOptionalPropertyTypes). */
export function toLoupeTags(tags: SampleTag[]): Tag[] {
  return tags.map((t) => (t.value ? { key: t.key, value: t.value } : { key: t.key }));
}

/** Resolve a greenfield Tag back to its SampleTag id for removal (key+value match). */
export function findSampleTagId(tags: SampleTag[], tag: Tag): number | undefined {
  const want = tag.value ?? "";
  return tags.find((t) => t.key === tag.key && (t.value ?? "") === want)?.id;
}
