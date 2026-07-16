import type { CorpusSample, Exposure } from "../../api";
import type { Tag } from "../ui";
import type { GalleryExposure } from "../components/ThumbnailGallery";
import { sampleDisplayName } from "../../lib/sample/displayName";
import { toGalleryExposures, toLoupeTags } from "./loupeAdapters";

export interface SampleRowModel {
  name: string;
  sampleId: string;
  exposures: GalleryExposure[];
  kept: number;
  total: number;
  dropped: number;
  tags: Tag[];
  phase: string | null | undefined;
  /** The representative exposure is declared form_factor → the status cell reads
   *  "Form factor" instead of "Not indexed". */
  formFactor: boolean;
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
  // Binary status: "Kept" means simply not-dropped (status !== rejected).
  // kept + dropped === total — there is no untriaged middle state.
  const dropped = exps.filter((e) => e.status === "rejected").length;
  const kept = total - dropped;
  return {
    name: sampleDisplayName(sample),
    sampleId: `#${sample.id}`,
    exposures: toGalleryExposures(exps),
    kept,
    total,
    dropped,
    tags: toLoupeTags(sample.tags),
    phase: sample.phase,
    formFactor: sample.assignment_state === "form_factor",
  };
}
