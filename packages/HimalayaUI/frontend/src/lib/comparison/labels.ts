/**
 * Shared per-member display-label resolver for the Compare review and edit
 * pages. Returns a Map keyed by `member.id` whose values follow the
 * documented fallback chain (issue #52, extended for #69):
 *
 *   1. `member.label_override` (user-set, wins outright)
 *   2. `${sample.display_name || sample.name} · ${exposure.filename}`
 *   3. `exposure.filename` alone
 *   4. `sample.display_name || sample.name` alone
 *   5. `Exposure #N` legacy fallback (still loading / cache cold)
 *   6. `(deleted exposure)` when `member.exposure_id === null`
 *
 * Both pages must agree on this chain — see ComparePage.tsx (review) and
 * ComparePageEdit.tsx (edit). Drift between the two surfaces the same
 * regression class that needed `comparisonHash8()` to fix #62.
 */
import type { SeriesMember, Exposure, Sample } from "../../api";

export function resolveDisplayLabels(
  members: readonly SeriesMember[],
  exposures: ReadonlyMap<number, Exposure>,
  samples: ReadonlyMap<number, Sample>,
): Map<number, string> {
  const map = new Map<number, string>();
  for (const m of members) {
    if (m.label_override) {
      map.set(m.id, m.label_override);
      continue;
    }
    if (m.exposure_id === null) {
      map.set(m.id, "(deleted exposure)");
      continue;
    }
    const exposure = exposures.get(m.exposure_id);
    const sample = exposure ? samples.get(exposure.sample_id) : undefined;
    // `||` (not `??`) so an empty-string display_name/name falls through to the
    // next fallback rather than rendering as a leading separator like " · run.dat".
    const sampleName = sample ? (sample.display_name || sample.name || null) : null;
    const filename = exposure?.filename || null;
    if (sampleName && filename) {
      map.set(m.id, `${sampleName} · ${filename}`);
    } else if (filename) {
      map.set(m.id, filename);
    } else if (sampleName) {
      map.set(m.id, sampleName);
    } else {
      map.set(m.id, `Exposure #${m.exposure_id}`);
    }
  }
  return map;
}
