import type { CorpusSample, Exposure } from "../../api";
import type { Tag } from "../ui";
import type { GalleryExposure } from "../components/ThumbnailGallery";
import { sampleDisplayName } from "../../lib/sample/displayName";
import { isSampleScreened } from "../../lib/sample/screened";
import { toGalleryExposures, toLoupeTags } from "./loupeAdapters";

export interface SampleRowModel {
  name: string;
  sampleId: string;
  screened: boolean;
  exposures: GalleryExposure[];
  kept: number;
  total: number;
  dropped: number;
  tags: Tag[];
  phase: string | null | undefined;
}

/** CorpusSample + its (possibly unloaded) exposures → SampleTableRow props.
 *  undefined exposures = "not yet fetched" → empty derivation (the page shows
 *  the boneyard skeleton meanwhile). */
export function toSampleRowModel(
  sample: CorpusSample,
  exposures: Exposure[] | undefined,
): SampleRowModel {
  const exps = exposures ?? [];
  const total = exps.length;
  // SA-SCREENED: "Kept" means EXPLICITLY accepted — the same truth the loupe
  // caption and the hero metric tell. The old `!== "rejected"` derivation
  // counted unscreened (null) frames as kept, which is the lie the Keep verb
  // fixes. Symmetrically, dropped counts only actual rejects: an unscreened
  // frame is neither, so kept + dropped can be < total on an untriaged corpus.
  const kept = exps.filter((e) => e.status === "accepted").length;
  const dropped = exps.filter((e) => e.status === "rejected").length;
  return {
    name: sampleDisplayName(sample),
    sampleId: `#${sample.id}`,
    screened: isSampleScreened(sample, exposures),
    exposures: toGalleryExposures(exps),
    kept,
    total,
    dropped,
    tags: toLoupeTags(sample.tags),
    phase: sample.phase,
  };
}
