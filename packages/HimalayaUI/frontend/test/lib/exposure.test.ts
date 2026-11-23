import { describe, it, expect } from "vitest";
import cases from "../fixtures/analyzed-exposure.json";
import { analysisOf, analyzedExposure, keptExposures, selectedFrames } from "../../src/lib/exposure";
import type { Exposure } from "../../src/api";

type Frame = Pick<Exposure, "id" | "status"> & { kind?: string; selected?: boolean; sources?: number[] };
type Case = { case: string; frames: Frame[]; analyzed: number | null };

function exposures(c: Case): Exposure[] {
  return c.frames.map((f) => ({
    sample_id: 1, filename: `f${f.id}`, kind: "file", image_path: null, ...f,
    sources: (f.sources ?? []).map((id) => ({ source_exposure_id: id, role: "signal" })),
  }) as unknown as Exposure);
}

describe.each(cases as Case[])("analyzed image: $case", (c) => {
  it("is what the shared rule picks", () => {
    expect(analyzedExposure(exposures(c))?.id ?? null).toBe(c.analyzed);
  });
});

describe("analysisOf", () => {
  it("lists the selected kept images in id order, never an average or a dropped image", () => {
    const list = exposures({ case: "", analyzed: null, frames: [
      { id: 3, status: null }, { id: 1, status: null }, { id: 2, status: "rejected" },
      { id: 4, status: null, selected: false }, { id: 9, status: null, kind: "combined", sources: [1, 3] },
    ] });
    expect(selectedFrames(list).map((e) => e.id)).toEqual([1, 3]);
    const a = analysisOf(list);
    expect(a.selected.map((e) => e.id)).toEqual([1, 3]);
    expect(a.analyzed?.id).toBe(9);
  });
});

describe("keptExposures", () => {
  it("drops rejected frames and keeps accepted and unmarked ones", () => {
    const out = keptExposures([{ status: "accepted" }, { status: "rejected" }, { status: null }]);
    expect(out.map((e) => e.status)).toEqual(["accepted", null]);
  });

  it("treats a not-yet-loaded list as empty", () => {
    expect(keptExposures(undefined)).toEqual([]);
  });
});
