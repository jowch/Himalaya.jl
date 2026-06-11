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

/** Map the per-sample exposure list to the filmstrip view-model. The SA-SCREENED
 *  tri-state holds: kept means EXPLICITLY accepted (the same truth as the loupe
 *  caption and the contact-sheet Kept count); a null status is unscreened —
 *  neither kept nor rejected. */
export function toGalleryExposures(exposures: Exposure[]): GalleryExposure[] {
  return exposures.map((e, i) => ({
    id: e.id,
    src: buildExposureImageUrl(e, { thumb: true }),
    frameNo: i + 1,
    representative: e.selected,
    rejected: e.status === "rejected",
    kept: e.status === "accepted",
  }));
}

/**
 * The rejection reason of a dropped frame, read from the `rejection_reason`
 * exposure-tag convention (key `"rejection_reason"`, `source = "manual"` — the
 * legacy Inspect page wrote these via POST /api/exposures/:id/tags; see
 * docs/superpowers/specs/2026-04-25-inspect-page-design.md). DISPLAY-ONLY here:
 * PATCH /api/exposures/:id/status accepts no reason field, so the loupe's drop
 * verb cannot set one without composing a second mutation — out of scope for
 * this surface. Multiple reason tags join with "; ". `undefined` when absent.
 */
export function rejectionReason(exposure: Exposure): string | undefined {
  const reasons = exposure.tags
    .filter((t) => t.key === "rejection_reason" && t.value.trim() !== "")
    .map((t) => t.value);
  return reasons.length > 0 ? reasons.join("; ") : undefined;
}

/**
 * "This exposure" metadata rows. The mockup's `integration`/`collected` rows
 * are OMITTED, not stubbed: the exposures API carries neither field (the #256
 * stub rows showed an eternal "—" on every frame — controls-don't-lie says a
 * row whose data source doesn't exist yet isn't shown). Reinstate them when
 * the backend actually serves the values. The legacy signal-meter row is
 * intentionally dropped too — the greenfield LoupeSidePanel has no signal
 * block. A dropped frame carrying a rejection_reason tag gets a quiet
 * `reason` row (show the reasoning; display-only, see `rejectionReason`).
 */
export function toMetaEntries(active: Exposure, exposures: Exposure[]): MetaEntry[] {
  const idx = exposures.findIndex((e) => e.id === active.id);
  const position = idx >= 0 ? `${idx + 1} of ${exposures.length}` : "—";
  const entries: MetaEntry[] = [{ key: "frame", value: position }];
  if (active.status === "rejected") {
    const reason = rejectionReason(active);
    if (reason !== undefined) entries.push({ key: "reason", value: reason });
  }
  return entries;
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
