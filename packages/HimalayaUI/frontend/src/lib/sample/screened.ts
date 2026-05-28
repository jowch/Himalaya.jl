import type { CorpusSample, Exposure } from "../../api";

/**
 * Whether a sample has been "screened" — its frames have been looked at and
 * triaged. M-2 affordance (the filled/hollow dot + unscreened-row tint).
 *
 * COORDINATION (#162): the authoritative `screened` flag is owned by #162's
 * backend, which is not yet wired. Until then this derives a sensible signal:
 *
 *   1. If the sample carries an explicit `screened` boolean (the #162 field,
 *      surfaced here via the optional `CorpusSample.screened`), trust it.
 *   2. Otherwise derive from the exposures: a sample reads as screened once
 *      every exposure has been triaged (a non-null status — accepted or
 *      rejected). A sample with no exposures is not screened.
 *
 * This keeps the visual affordance honest against today's data and flips to
 * the real flag automatically when #162 lands (step 1 wins).
 */
export function isSampleScreened(
  sample: Pick<CorpusSample, "screened"> & Partial<CorpusSample>,
  exposures: readonly Exposure[] | undefined,
): boolean {
  if (typeof sample.screened === "boolean") return sample.screened;
  if (!exposures || exposures.length === 0) return false;
  return exposures.every((e) => e.status !== null);
}
